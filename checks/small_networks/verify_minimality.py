#!/usr/bin/env python3
"""Exact finite checks accompanying the section "Five edges are minimal".

Standard library only. All arithmetic is exact (fractions.Fraction).
These are finite checks of the combinatorial and algebraic ingredients of the
proof; the parameter-uniform theorem is proved analytically in the text.

Usage: python3 verify_minimality.py [--output PATH]
"""
import argparse
import itertools
import json
import sys
from fractions import Fraction as Fr

FAILURES = []


def check(cond, msg):
    if not cond:
        FAILURES.append(msg)
    return bool(cond)


# ---------------------------------------------------------------- linear algebra
def rank(rows):
    m = [[Fr(x) for x in r] for r in rows if any(x != 0 for x in r)]
    if not m:
        return 0
    ncol = len(m[0])
    r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c] / m[r][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        r += 1
        if r == len(m):
            break
    return r


# ---------------------------------------------------------------- graphs
def weakly_connected(n, edges):
    seen = {0}
    changed = True
    while changed:
        changed = False
        for u, v in edges:
            if (u in seen) != (v in seen):
                seen |= {u, v}
                changed = True
    return len(seen) == n


def incidence_rank(n, edges):
    rows = []
    for x in range(n):
        rows.append([(1 if v == x else 0) - (1 if u == x else 0) for u, v in edges])
    return rank(rows)


def underlying_is_tree(n, edges):
    pairs = {frozenset(e) for e in edges}
    return len(pairs) == n - 1  # connected + n-1 pairs => tree


def has_directed_cycle(n, edges):
    return not any(all(order.index(u) < order.index(v) for u, v in edges)
                   for order in itertools.permutations(range(n)))


def canonical(n, edges):
    return min(tuple(sorted((p[u], p[v]) for u, v in edges))
               for p in itertools.permutations(range(n)))


def simple_paths(edges, o, d):
    """All simple directed o-d paths as tuples of edge indices."""
    out = []

    def dfs(x, visited, path):
        if x == d:
            out.append(tuple(path))
            return
        for i, (u, v) in enumerate(edges):
            if u == x and v not in visited:
                dfs(v, visited | {v}, path + [i])
    dfs(o, {o}, [])
    return out


def vec(path, m):
    v = [0] * m
    for i in path:
        v[i] += 1
    return v


def topology_checks():
    res = {}
    for n in (2, 3, 4, 5):
        arcs = [(u, v) for u in range(n) for v in range(n) if u != v]
        counts = {"connected": 0, "tree": 0, "nontree": 0,
                  "nontree_kernel_le_1": 0, "nontree_kernel_2": 0}
        classes_nontree_e4 = {}
        for e in range(1, 5):
            for edges in itertools.combinations_with_replacement(arcs, e):
                if not weakly_connected(n, edges):
                    continue
                counts["connected"] += 1
                kdim = e - incidence_rank(n, edges)
                check(kdim == e - n + 1, f"kernel dim formula n={n} {edges}")
                if underlying_is_tree(n, edges):
                    counts["tree"] += 1
                    continue
                counts["nontree"] += 1
                check(n != 2, "two-vertex graph with non-tree underlying graph")
                if kdim <= 1:
                    counts["nontree_kernel_le_1"] += 1
                    check(n >= 4 or e <= 3, f"kernel<=1 outside predicted range {n},{edges}")
                else:
                    counts["nontree_kernel_2"] += 1
                    check(n == 3 and e == 4, f"kernel 2 outside n=3,e=4: {n},{edges}")
                    cyc = has_directed_cycle(n, edges)
                    if not cyc:
                        pairs = [frozenset(x) for x in edges]
                        mult = sorted(pairs.count(p) for p in set(pairs))
                        check(mult == [1, 1, 2], "acyclic non-tree is not triangle+dup")
                    key = canonical(n, edges)
                    classes_nontree_e4[key] = "acyclic_transitive_triangle_with_duplicated_side" if not cyc else "has_directed_cycle"
        entry = dict(counts)
        if n == 3:
            entry["nontree_e4_isomorphism_classes"] = {
                "acyclic": sum(1 for t in classes_nontree_e4.values() if t.startswith("acyclic")),
                "with_directed_cycle": sum(1 for t in classes_nontree_e4.values() if t == "has_directed_cycle"),
            }
            check(entry["nontree_e4_isomorphism_classes"]["acyclic"] == 3,
                  "expected three acyclic classes (duplicate on f, h or g)")
        res[f"n={n}"] = entry
    # directed triangle: every ordered pair has at most one simple path
    cyc = [(0, 1), (1, 2), (2, 0)]
    ok = all(len(simple_paths(cyc, o, d)) == 1 for o in range(3) for d in range(3) if o != d)
    check(ok, "directed triangle has a pair with several routes")
    res["directed_triangle_unique_routes"] = ok
    return res


# ---------------------------------------------------------------- triangle with duplicated side
S, U, V = 0, 1, 2
SIDE_ARC = {"f": (S, V), "h": (S, U), "g": (U, V)}


def dup_triangle(dup):
    """Edges and side labels; the duplicate is the last edge."""
    sides = ["f", "h", "g", dup]
    edges = [SIDE_ARC[x] for x in sides]
    return edges, sides


def is_pi(v, sides, dup):
    nz = [i for i, x in enumerate(v) if x != 0]
    return len(nz) == 2 and all(sides[i] == dup for i in nz) and sorted(v[i] for i in nz) == [-1, 1]


def is_tau(v, sides):
    w = [x for x in v]
    if any(abs(x) > 1 for x in w):
        return False
    plus = [sides[i] for i, x in enumerate(w) if x == 1]
    minus = [sides[i] for i, x in enumerate(w) if x == -1]
    return (sorted(plus) == ["g", "h"] and minus == ["f"]) or (sorted(minus) == ["g", "h"] and plus == ["f"])


def classification_checks():
    out = {}
    pairs = {"sv": (S, V), "su": (S, U), "uv": (U, V)}
    for dup in ("f", "h", "g"):
        edges, sides = dup_triangle(dup)
        m = len(edges)
        check(m - incidence_rank(3, edges) == 2, "kernel of dup triangle not 2")
        # only sv, su, uv have routes
        for o in range(3):
            for d in range(3):
                if o != d:
                    has = len(simple_paths(edges, o, d)) > 0
                    check(has == ((o, d) in pairs.values()), "unexpected route existence")
        routes = {k: simple_paths(edges, *pairs[k]) for k in pairs}
        n_families = 0
        dims = {0: 0, 1: 0, 2: 0}
        rank1_types = {"pi": 0, "tau": 0}
        for r in range(1, 4):
            for comms in itertools.combinations(sorted(pairs), r):
                choices = [[c for j in range(1, len(routes[k]) + 1)
                            for c in itertools.combinations(routes[k], j)] for k in comms]
                for fam in itertools.product(*choices):
                    n_families += 1
                    gens = []
                    for rs in fam:
                        for P, Q in itertools.combinations(rs, 2):
                            g = [a - b for a, b in zip(vec(P, m), vec(Q, m))]
                            check(is_pi(g, sides, dup) or is_tau(g, sides), "generator not pi/tau")
                            gens.append(g)
                    k = rank(gens)
                    check(k <= 2, "rank > 2")
                    dims[k] += 1
                    if k == 1:
                        nz = [g for g in gens if any(g)]
                        if any(is_pi(g, sides, dup) for g in nz):
                            check(all(is_pi(g, sides, dup) for g in nz), "mixed pi/tau with rank 1")
                            rank1_types["pi"] += 1
                        else:
                            taus = {tuple(g) if g[0] == -1 else tuple(-x for x in g) for g in nz}
                            check(len(taus) == 1, "two distinct triangle differences with rank 1")
                            rank1_types["tau"] += 1
        # complete regimes from all cost vectors in {1..4}^4
        n_cost = 0
        regimes_by_dim = {0: 0, 1: 0, 2: 0}
        for cost in itertools.product(range(1, 5), repeat=m):
            for r in range(1, 4):
                for comms in itertools.combinations(sorted(pairs), r):
                    n_cost += 1
                    fam = []
                    for k in comms:
                        pc = [sum(cost[i] for i in P) for P in routes[k]]
                        mn = min(pc)
                        fam.append({P for P, c in zip(routes[k], pc) if c == mn})
                    # side-structure formula, Lemma sides(a)
                    msid = {x: min(cost[i] for i in range(m) if sides[i] == x) for x in "fhg"}
                    Mset = {x: {i for i in range(m) if sides[i] == x and cost[i] == msid[x]} for x in "fhg"}
                    delta = msid["h"] + msid["g"] - msid["f"]
                    for k, rs in zip(comms, fam):
                        if k == "su":
                            pred = {(i,) for i in Mset["h"]}
                        elif k == "uv":
                            pred = {(i,) for i in Mset["g"]}
                        else:
                            direct = {(i,) for i in Mset["f"]}
                            indirect = {(i, j) for i in Mset["h"] for j in Mset["g"]}
                            pred = direct if delta > 0 else indirect if delta < 0 else direct | indirect
                        check(rs == pred, f"side-structure formula fails {dup} {cost} {k}")
                    gens = [[a - b for a, b in zip(vec(P, m), vec(Q, m))]
                            for rs in fam for P, Q in itertools.combinations(sorted(rs), 2)]
                    k = rank(gens)
                    active = any(sides[i] == dup for rs in fam for P in rs for i in P)
                    cond_a = active and len(Mset[dup]) == 2
                    cond_b = ("sv" in comms) and delta == 0
                    check(k == int(cond_a) + int(cond_b), f"rank formula fails {dup} {cost} {comms}")
                    if k == 1 and cond_b:
                        check(len(Mset[dup]) == 1, "type T with two cheapest duplicated edges")
                    regimes_by_dim[k] += 1
        out[f"duplicate_on_{dup}"] = {
            "route_counts": {k: len(v) for k, v in routes.items()},
            "route_families_enumerated": n_families,
            "families_by_rank": {str(k): v for k, v in dims.items()},
            "rank_one_types": rank1_types,
            "complete_regimes_checked": n_cost,
            "complete_regimes_by_rank": {str(k): v for k, v in regimes_by_dim.items()},
        }
    return out


# ---------------------------------------------------------------- algebra
class LCG:
    def __init__(self, seed):
        self.s = seed

    def nxt(self):
        self.s = (1103515245 * self.s + 12345) % (2 ** 31)
        return self.s

    def rat(self, lo=1, hi=9):
        num = lo + self.nxt() % (hi - lo + 1)
        den = 1 + self.nxt() % 5
        return Fr(num, den)


def curvature_checks():
    rng = LCG(20261022)
    out = []
    for p in (2, 3, 4):
        trials = 0
        for m in (2, 3):
            for _ in range(20):
                xs = [rng.rat() for _ in range(m)]
                bs = [rng.rat() for _ in range(m)]
                lam = 1 + max(b * x ** p for b, x in zip(bs, xs)) + rng.rat()
                aas = [lam - b * x ** p for b, x in zip(bs, xs)]
                check(all(a > 0 for a in aas), "nonpositive offset")
                ds = [p * b * x ** (p - 1) for b, x in zip(bs, xs)]
                c2 = [p * (p - 1) * b * x ** (p - 2) for b, x in zip(bs, xs)]
                T = sum(xs)
                Np = sum(1 / d for d in ds)
                Npp = -sum(s / d ** 3 for s, d in zip(c2, ds))
                lp, lpp = 1 / Np, -Npp / Np ** 3
                eta = Np
                lhs = lpp - (p - 1) * lp / T
                rhs = Fr(p - 1) / (T * eta ** 3) * (T * sum(1 / (x * d * d) for x, d in zip(xs, ds)) - eta ** 2)
                check(lhs == rhs, "curvature identity")
                check(all(x * d == p * (lam - a) for x, d, a in zip(xs, ds, aas)), "x d = p(l-a)")
                eq = len(set(aas)) == 1
                check((rhs == 0) == eq and rhs >= 0, "curvature sign/equality")
                trials += 1
            # equal offsets
            ws = [rng.rat() for _ in range(m)]
            bs = [1 / w ** p for w in ws]
            a = rng.rat()
            xs = ws
            lam = a + 1
            ds = [p * b * x ** (p - 1) for b, x in zip(bs, xs)]
            T = sum(xs)
            eta = sum(1 / d for d in ds)
            resid = T * sum(1 / (x * d * d) for x, d in zip(xs, ds)) - eta ** 2
            kappa_inv_root = sum(ws)  # sum b^{-1/p}
            check(resid == 0, "equal offsets residual")
            check(lam == a + T ** p / kappa_inv_root ** p, "equal offsets latency formula")
            trials += 1
        # dimension-lemma identity D q = p (c0 - a)/t0
        t0 = rng.rat()
        q = [rng.rat() for _ in range(4)] + [Fr(0)]
        aa = [rng.rat() for _ in range(5)]
        bb = [rng.rat() for _ in range(5)]
        z0 = [t0 * x for x in q]
        c0 = [a + b * z ** p for a, b, z in zip(aa, bb, z0)]
        Dq = [p * b * z ** (p - 1) * x for b, z, x in zip(bb, z0, q)]
        check(Dq == [p * (c - a) / t0 for c, a in zip(c0, aa)], "Dq identity")
        # type-T velocities and stationary-leg curvature
        for _ in range(20):
            R, A, B = rng.rat(), rng.rat(), rng.rat()
            ds = [rng.rat() for _ in range(3)]  # d_f, d_h, d_g
            J = sum(ds)
            xp = (R * ds[0] - A * ds[1] - B * ds[2]) / J
            check(R - xp == ((R + A) * ds[1] + (R + B) * ds[2]) / J, "y_f'")
            check(A + xp == ((R + A) * ds[0] - (B - A) * ds[2]) / J, "y_h'")
            check(B + xp == ((R + B) * ds[0] - (A - B) * ds[1]) / J, "y_g'")
        for _ in range(20):
            t, H = rng.rat(), rng.rat()
            L = rng.rat() + H / t  # ensures y_f = L t - H > 0
            M = rng.rat()
            yf, yg = L * t - H, M * t + H
            bf = rng.rat()
            df = p * bf * yf ** (p - 1)
            bg = L * df / (M * p * yg ** (p - 1))  # enforce L d_f = M d_g
            dg = p * bg * yg ** (p - 1)
            bh = rng.rat()
            dh = p * bh * H ** (p - 1)
            J = df + dh + dg
            check(L * df == M * dg, "stationarity constraint")
            f2 = p * (p - 1) * bf * yf ** (p - 2)
            g2 = p * (p - 1) * bg * yg ** (p - 2)
            Hpp = (L * L * f2 - M * M * g2) / J
            U = L * df
            check(Hpp == (p - 1) * U * H * (L + M) / (J * yf * yg), "H'' formula")
            check(Hpp > 0, "H'' positive")
            check(L * yg - M * yf == H * (L + M), "L y_g - M y_f identity")
        out.append({"p": p, "random_bank_trials": trials, "status": "ok"})
    return out


def heterogeneous_checks():
    # G(t) = c_h(t) + c_g(t) - c_f(t) with c_f=1+t^4, c_h=1/3+t^2, c_g=1/3+t^6/3
    G = {0: Fr(1, 3) + Fr(1, 3) - 1, 2: Fr(1), 6: Fr(1, 3), 4: Fr(-1)}
    target = {0: Fr(-1, 3), 2: Fr(1), 4: Fr(-1), 6: Fr(1, 3)}  # (t^2-1)^3/3
    check(G == target, "heterogeneous triangle gap is not (t^2-1)^3/3")
    signs = ["+" if G[k] > 0 else "-" for k in sorted(G)]
    changes = sum(1 for a, b in zip(signs, signs[1:]) if a != b)
    check(signs == ["-", "+", "-", "+"] and changes == 3, "sign pattern")
    return {"gap_coefficients_by_exponent": {str(k): str(v) for k, v in sorted(G.items())},
            "equals_(t^2-1)^3/3": G == target, "sign_pattern": "".join(signs),
            "sign_changes": changes}


# ---------------------------------------------------------------- reversal symmetry
def min_sets(edges, cost, comms):
    res = {}
    for (o, d) in comms:
        ps = simple_paths(edges, o, d)
        pc = [sum(cost[i] for i in P) for P in ps]
        mn = min(pc)
        res[(o, d)] = (sorted(frozenset(P) for P in ps),
                       sorted(frozenset(P) for P, c in zip(ps, pc) if c == mn))
    return res


def reversal_checks():
    edges, sides = dup_triangle("h")
    rev = [(v, u) for u, v in edges]
    comms = [(S, V), (S, U), (U, V)]
    rng = LCG(7)
    n_ok = 0
    for _ in range(25):
        cost = [rng.rat(1, 4) for _ in edges]
        a = min_sets(edges, cost, comms)
        b = min_sets(rev, cost, [(d, o) for o, d in comms])
        for (o, d) in comms:
            check(a[(o, d)][0] == b[(d, o)][0], "route edge sets differ under reversal")
            check(a[(o, d)][1] == b[(d, o)][1], "minimum sets differ under reversal")
        n_ok += 1
    # exchange: reversal then s<->v maps dup on h to dup on g, su commodity to uv
    swap = {S: V, V: S, U: U}
    ex = [(swap[u], swap[v]) for u, v in rev]
    ex_sides = []
    for (u, v) in ex:
        ex_sides.append([k for k, arc in SIDE_ARC.items() if arc == (u, v)][0])
    check(sorted(ex_sides) == ["f", "g", "g", "h"], "exchange does not move duplicate to g")
    check(ex_sides[0] == "f", "exchange does not fix f")
    su_routes = simple_paths(edges, S, U)
    ex_uv = simple_paths(ex, U, V)
    check(sorted(map(frozenset, su_routes)) == sorted(map(frozenset, ex_uv)),
          "exchange does not map s->u routes to u->v routes")
    return {"instances": n_ok, "exchange_side_images": ex_sides}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default="./minimality_checks.json")
    args = ap.parse_args()
    result = {
        "topology": topology_checks(),
        "route_space_classification": classification_checks(),
        "curvature_and_triangle_algebra": curvature_checks(),
        "heterogeneous_triangle": heterogeneous_checks(),
        "reversal_symmetry": reversal_checks(),
    }
    result["failures"] = FAILURES
    result["status"] = "PASS" if not FAILURES else "FAIL"
    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(result, fh, sort_keys=True, indent=2)
        fh.write("\n")
    if FAILURES:
        print(f"FAIL: {len(FAILURES)} check(s) failed; first: {FAILURES[0]}")
        sys.exit(1)
    print("PASS: topology, route-space classification, curvature (p=2,3,4), "
          "triangle algebra, heterogeneous gap, reversal symmetry")


if __name__ == "__main__":
    main()
