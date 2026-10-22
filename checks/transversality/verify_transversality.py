#!/usr/bin/env python3
"""Exact checks from physical quartic coefficients; standard library only.

The general results are proved in results.md. This script tests the displayed
constants, non-symmetric rays, derivative anisotropy and finite moment ranks.
It reuses the main bundle's exact power-series arithmetic, but reconstructs
the physical scaled costs and performs new implicit series substitutions.
Checks deliberately do not use assert, so they also run under python -O -S.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import importlib.util
import json

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("exact_arithmetic", ROOT / "code" / "verify_exact.py")
exact = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exact)
mul, compose, power, check = exact.mul, exact.compose, exact.power, exact.check


def determinant(rows):
    a = [[Q(x) for x in row] for row in rows]
    n = len(a)
    check(all(len(row) == n for row in a), "Determinant requires square matrix")
    d = Q(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return Q(0)
        if pivot != j:
            a[pivot], a[j] = a[j], a[pivot]
            d = -d
        d *= a[j][j]
        p = a[j][j]
        for i in range(j + 1, n):
            c = a[i][j] / p
            for k in range(j + 1, n):
                a[i][k] -= c * a[j][k]
    return d


def implicit_on_path(f, g, a, b, v, w, n):
    """Baseline bank loads a,b; path increments v*delta,w*delta."""
    x = [Q(0)] * (n + 1)
    jac = -f[1] / a - g[1] / b

    def bank_values(x, k):
        pf = [-z / a for z in x[:k + 1]]
        pg = [z / b for z in x[:k + 1]]
        pf[1] += v / a
        pg[1] += w / b
        return compose(f, pf, k), compose(g, pg, k)

    def residual(x, k):
        ff, gg = bank_values(x, k)
        hh = power(x, 4, k)
        return [ff[j] - gg[j] - hh[j] for j in range(k + 1)]

    for k in range(1, n + 1):
        x[k] = -residual(x, k)[k] / jac
    check(all(z == 0 for z in residual(x, n)), "Implicit connector residual")

    def social_cost(x):
        ff, gg = bank_values(x, n)
        ff[0], gg[0] = Q(3), Q(2)
        ecf, ecg = mul([a, v], ff, n), mul([b, w], gg, n)
        return [ecf[j] + ecg[j] for j in range(n + 1)]

    E, E0 = social_cost(x), social_cost([Q(0)] * (n + 1))
    return x, [E[j] - E0[j] for j in range(n + 1)]


def asymmetric_case(m, a, b):
    n = 2 * m - 1
    N = Q(2 ** (m - 1))
    edges = exact.physical_links(m)
    for e in edges:
        factor = a / N if e["group"] == "F" else b / N
        e["b"] /= factor ** 4
        e["baseline_load"] *= factor
    qf = exact.throughput_series(edges, "F", n)
    qg = exact.throughput_series(edges, "G", n)
    f = exact.inverse_series(qf, a, n)
    g = exact.inverse_series(qg, b, n)
    x, correction = implicit_on_path(f, g, a, b, a, b, n)
    d = f[m] - g[m]
    L = 1 / qf[1] + 1 / qg[1]
    check(d > 0 and L > 0, "Asymmetric contact signs")
    check(all(z == 0 for z in x[:m]) and x[m] == d / L,
          "Asymmetric flow contact")
    expected = -m * d * d / L
    check(all(z == 0 for z in correction[:2*m-1]) and correction[-1] == expected,
          "Asymmetric aggregate contact")
    return {"m": m, "demand_a": str(a), "demand_b": str(b),
            "gap_coefficient": str(d), "connector_coefficient": str(x[m]),
            "cost_coefficient": str(correction[-1])}


def codimension_case(m):
    N = Q(2 ** (m - 1))
    edges = exact.physical_links(m)
    groups = {s: [j for j in range(m + 1) if edges[j]["group"] == s]
              for s in ["F", "G"]}
    anchors = {s: rows[-1] for s, rows in groups.items()}
    free = [j for j in range(m + 1) if j not in anchors.values()]
    matrix = []
    for k in range(1, m):
        row = []
        for j in free:
            s = edges[j]["group"]
            sign = 1 if s == "F" else -1
            row.append(Q(sign, 1) * ((j + 1)**k - (anchors[s]+1)**k) / N)
        matrix.append(row)
    det = determinant(matrix)
    check(det != 0, "Pinned physical moment Jacobian loses rank")
    for k in range(m):
        moment = sum(Q((-1)**j * comb(m,j) * (j+1)**k, N) for j in range(m+1))
        check(moment == 0, "Reference moment does not vanish")
    return {"m": m, "dimension_of_fixed_total_chart": m-1,
            "rank": m-1, "restricted_moment_determinant": str(det),
            "matrix": [[str(z) for z in row] for row in matrix]}


def transverse_case(m, v, w):
    N = Q(2**(m-1))
    n = 3
    edges = exact.physical_links(m)
    qf = exact.throughput_series(edges, "F", n)
    qg = exact.throughput_series(edges, "G", n)
    f = exact.inverse_series(qf, N, n)
    g = exact.inverse_series(qg, N, n)
    x, correction = implicit_on_path(f, g, N, N, v, w, n)
    f1 = f[1]/N
    f2 = 2*f[2]/N**2
    chi = f1/2 + N*f2/4
    stated = Q(4*(2*(m+2)**2+3*(m+1)*(m+4)), 2*N*(m+2)**3)
    check(chi == stated and chi > 0, "Anisotropy coefficient formula")
    check(x[1] == (v-w)/2, "Transverse connector derivative")
    check(correction[0] == correction[1] == 0, "First-order cost discontinuity")
    check(correction[2] == -chi*(v-w)**2, "Transverse cost quadratic coefficient")
    return {"m": m, "path_v": str(v), "path_w": str(w),
            "chi": str(chi), "connector_active_derivative": str(x[1]),
            "cost_active_minus_inactive_second_derivative": str(correction[2]*2),
            "active_side": "right" if v>w else "left"}


def main():
    output = {
        "method": "exact Fraction arithmetic from physical quartic edges",
        "asymmetric_rays": [asymmetric_case(m, Q(a), Q(b))
                            for m,a,b in [(3,2,7),(3,Q(1,3),Q(11,2)),(3,3,5),(5,2,7)]],
        "physical_codimension": [codimension_case(m) for m in [3,5,7]],
        "transverse_paths": [transverse_case(m,Q(v),Q(w))
                             for m,v,w in [(3,1,0),(3,2,-3),(3,-1,4),(5,1,0),(7,1,0)]],
    }
    A,d,B = Q(5,2),Q(336,625),Q(8,5)
    threshold_cube = B/d
    kink_coefficient_cube = (A*9/2)**3*d**2*B**4
    check(threshold_cube == Q(125,42), "Unequal-demand threshold constant")
    check(kink_coefficient_cube == Q(5267275776,1953125), "Unequal-demand kink constant")
    output["demand_perturbation_m3_p4"] = {
        "absolute_threshold_prefactor_cubed": str(threshold_cube),
        "kink_prefactor_cubed": str(kink_coefficient_cube),
        "kink_prefactor_approximate": float(kink_coefficient_cube)**(1/3),
    }
    destination = Path(__file__).with_name("exact_checks.json")
    destination.write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps({"status":"pass", "asymmetric_rays":4,
                      "physical_moment_ranks":3, "transverse_paths":5,
                      "output":str(destination)}, indent=2))


if __name__ == "__main__":
    main()
