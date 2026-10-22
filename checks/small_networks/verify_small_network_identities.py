#!/usr/bin/env python3
"""Exact finite fingerprints for the analytic small-network proof.

This verifies algebraic examples and the finite three-vertex topology reduction.
It is not a numerical proof over latency parameters or demand vectors.
No third-party dependencies; checks remain active under python -O -S.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations
from pathlib import Path
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def weakly_connected(edges):
    seen = {0}
    while True:
        new = seen | {v for u, v in edges if u in seen} | {u for u, v in edges if v in seen}
        if new == seen:
            return len(seen) == 3
        seen = new


def acyclic(edges):
    return any(all(order.index(u) < order.index(v) for u, v in edges)
               for order in permutations(range(3)))


def canonical(edges):
    return min(tuple(sorted((perm[u], perm[v]) for u, v in edges))
               for perm in permutations(range(3)))


def topology_fingerprints():
    arcs = [(u, v) for u in range(3) for v in range(3) if u != v]
    all_graphs = list(combinations_with_replacement(arcs, 4))
    connected = [edges for edges in all_graphs if weakly_connected(edges)]
    dags = [edges for edges in connected if acyclic(edges)]
    classes = {}
    for edges in dags:
        simple_pairs = {tuple(sorted(e)) for e in edges}
        typ = "tree_of_parallel_banks" if len(simple_pairs) == 2 else "transitive_triangle_with_one_duplicated_side"
        check(len(simple_pairs) in (2, 3), "Connected three-vertex graph has wrong simple-pair count")
        if len(simple_pairs) == 3:
            counts = sorted(sum(tuple(sorted(e)) == pair for e in edges) for pair in simple_pairs)
            check(counts == [1, 1, 2], "Unexpected triangular multiplicities")
        classes.setdefault(canonical(edges), typ)
    return {
        "scope": "loopless directed multigraphs on three labelled vertices, four edges; no demand or cost enumeration",
        "labelled_multigraphs": len(all_graphs),
        "weakly_connected": len(connected),
        "connected_directed_acyclic": len(dags),
        "isomorphism_classes_of_connected_dags": len(classes),
        "classes": [{"edges": list(edges), "type": typ} for edges, typ in sorted(classes.items())],
    }


def algebraic_fingerprints():
    out = []
    for p in (2, 3, 4):
        # Two physical parallel edges have loads 1 and 3 and common delay 10.
        # Their individual derivatives are both one. Reconstruct a,b from
        # these constraints, independently of the curvature expression.
        xs = [F(1), F(3)]
        bs = [F(1, p) / x ** (p-1) for x in xs]
        aas = [F(10) - b*x**p for b, x in zip(bs, xs)]
        ds = [p*b*x**(p-1) for b, x in zip(bs, xs)]
        seconds = [p*(p-1)*b*x**(p-2) for b, x in zip(bs, xs)]
        throughput_prime = sum(1/d for d in ds)
        throughput_second = -sum(s/d**3 for s, d in zip(seconds, ds))
        latency_prime = 1/throughput_prime
        latency_second = -throughput_second / throughput_prime**3
        excess = latency_second - (p-1)*latency_prime/sum(xs)
        cauchy_residual = sum(xs)*sum(1/(x*d*d) for x, d in zip(xs, ds)) - throughput_prime**2
        formula_excess = F(p-1)/(sum(xs)*throughput_prime**3)*cauchy_residual
        check(excess == formula_excess == F(p-1, 24), "Parallel curvature mismatch")
        check(all(a > 0 for a in aas), "Nonpositive offset")

        # Old triangle: demands R=2,A=0,B=2; x=1 at t=1.
        # f,h,g loads are 1,1,3; all first derivatives are one.
        # a_h=a_g=10, a_f=20+3/p enforce f=h+g exactly.
        loads = [F(1), F(1), F(3)]
        coeff_b = [F(1,p), F(1,p), F(1,p*3**(p-1))]
        coeff_a = [F(20)+F(3,p), F(10), F(10)]
        costs = [a+b*z**p for a,b,z in zip(coeff_a, coeff_b, loads)]
        slopes = [p*b*z**(p-1) for b,z in zip(coeff_b, loads)]
        curvatures = [p*(p-1)*b*z**(p-2) for b,z in zip(coeff_b, loads)]
        check(costs[0] == costs[1]+costs[2], "Triangle route costs differ")
        check(slopes == [F(1)]*3, "Triangle slopes differ")
        R, A, B = F(2), F(0), F(2)
        L, M = R+A, B-A
        J = sum(slopes)
        h_prime = (L*slopes[0]-M*slopes[2])/J
        check(h_prime == 0, "Triangle target is not stationary")
        implicit_h_second = (L*L*curvatures[0]-M*M*curvatures[2])/J
        formula_h_second = F(p-1)*L*slopes[0]*loads[1]*(L+M)/(J*loads[0]*loads[2])
        check(implicit_h_second == formula_h_second == F(8*(p-1),9), "Stationary leg curvature mismatch")
        out.append({
            "p": p,
            "parallel": {
                "a": [str(x) for x in aas], "b": [str(x) for x in bs],
                "loads": [str(x) for x in xs], "latency": "10",
                "L_prime": str(latency_prime), "L_second": str(latency_second),
                "curvature_excess": str(excess), "cauchy_residual": str(cauchy_residual),
            },
            "triangle": {
                "edge_order": ["f_direct", "h_connector", "g_shared"],
                "a": [str(x) for x in coeff_a], "b": [str(x) for x in coeff_b],
                "loads": [str(x) for x in loads], "demand_coefficients_R_A_B": ["2", "0", "2"],
                "H_prime": str(h_prime), "H_second": str(implicit_h_second),
            },
        })
    return out


def main():
    result = {
        "status": "all finite exact checks passed",
        "proof_boundary": "The parameter-uniform theorem is proved analytically in global_minimality_proof.md; these are independent finite fingerprints.",
        "topology": topology_fingerprints(),
        "algebra": algebraic_fingerprints(),
    }
    output = Path(__file__).resolve().with_name("small_network_verification.json")
    output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "output": str(output),
                      "connected_dag_classes": result["topology"]["isomorphism_classes_of_connected_dags"],
                      "exact_exponents_checked": [x["p"] for x in result["algebra"]]}, indent=2))


if __name__ == "__main__":
    main()
