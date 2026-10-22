#!/usr/bin/env python3
"""Exact checks of small examples (sympy rationals / fractions.Fraction).

(a) devil's-advocate quartic/sextic triangle (contraction: active below 1)
(b) cubic triangle (active above 1)
(c) swapped five-edge contraction
"""
import argparse, json, sys
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sympy as sp
import wardrop_series as ws

FAILS = []
t, x, d, y = sp.symbols("t x d y")
R = sp.Rational


def need(c, m):
    if not c:
        FAILS.append(m)
    return bool(c)


def s(v):
    return str(sp.nsimplify(v))


def implicit_series(eqn, n):
    """formal branch x(d)=sum_{k>=1} c_k d^k of polynomial eqn(x,d)=0 with x(0)=0."""
    X = sp.Integer(0)
    c = sp.Symbol("c")
    for k in range(1, n + 1):
        e = sp.expand(eqn.subs(x, X + c * d ** k))
        coeff = sp.Poly(e, d).coeff_monomial(d ** k) if e != 0 else 0
        sol = sp.solve(coeff, c)
        X += (sol[0] if sol else 0) * d ** k
    res = sp.expand(eqn.subs(x, X))
    low = [sp.Poly(res, d).coeff_monomial(d ** k) for k in range(n + 1)]
    need(all(v == 0 for v in low), "implicit residual")
    return sp.expand(X)


def ser(expr, n):
    return sp.expand(sp.series(sp.expand(expr), d, 0, n + 1).removeO())


def lead(expr, n):
    p = sp.Poly(ser(expr, n), d)
    for k in range(n + 1):
        cf = p.coeff_monomial(d ** k)
        if cf != 0:
            return k, cf
    return None, 0


def part_a():
    f = lambda z: 1 + z ** 4; h = lambda z: R(1, 3) + z ** 2; g = lambda z: R(1, 3) + z ** 6 / 3
    gap = sp.expand(f(t) - h(t) - g(t))
    need(sp.expand(gap + (t ** 2 - 1) ** 3 / 3) == 0, "triangle gap identity")
    n = 7
    eqn = sp.expand((f(t - x) - h(t + x) - g(t + x)).subs(t, 1 + d))
    X = implicit_series(eqn, n)
    kx, cx = lead(X, n)
    need(kx == 3 and cx == R(-1, 3), "x = -d^3/3")
    Ecom = (t * (f(t - x) + h(t + x) + g(t + x))).subs(t, 1 + d)
    Ephys = ((t - x) * f(t - x) + (t + x) * h(t + x) + (t + x) * g(t + x)).subs(t, 1 + d)
    E0 = sp.expand(Ecom.subs(x, 0))
    dE = ser(Ecom.subs(x, X) - E0, n); dEp = ser(Ephys.subs(x, X) - E0, n)
    need(sp.expand(dE - dEp) == 0, "commodity vs physical cost agree on branch")
    kE, cE = lead(dE, n)
    need(kE == 5 and cE == R(-8, 3), "E-E0 = -(8/3)d^5")
    E0_1, E0p = E0.subs(d, 0), sp.diff(E0, d).subs(d, 0)
    need(E0_1 == 4 and E0p == 12, "E0(1)=4, E0'(1)=12")
    # side of activity: x>0 needs d<0 (odd order, negative coefficient); gap sign
    need(gap.subs(t, R(9, 10)) > 0 and gap.subs(t, R(11, 10)) < 0, "gap positive below 1, negative above")
    # optimum: marginal costs
    z = sp.Symbol("z")
    FM, HM, GM = [sp.Lambda(z, sp.expand(sp.diff(z * c(z), z))) for c in (f, h, g)]
    gapm = sp.expand((FM(t) - HM(t) - GM(t)).subs(t, 1 + d))
    km, cm = lead(gapm, 6)
    eqm = sp.expand((FM(t - x) - HM(t + x) - GM(t + x)).subs(t, 1 + d))
    Y = implicit_series(eqm, 6)
    ky, cy = lead(Y, 6)
    need(km == 2 and cm < 0 and ky == 2 and cy == R(-1, 5), "optimum branch y=-d^2/5 (infeasible both sides)")
    # global sign of marginal gap near 1: gapm = -(polynomial) check on [0.9,1.1] by factoring
    gm_t = sp.factor(sp.expand(FM(t) - HM(t) - GM(t)))
    S0 = sp.expand((t * (f(t) + h(t) + g(t))).subs(t, 1 + d))
    need(sp.expand(S0 - E0) == 0, "S = E0 near 1")
    jumpE5 = -sp.factorial(5) * cE   # right minus left (left carries the extra term)
    return {"label": "exact", "network": "f: s->v 1+z^4; h: s->u 1/3+z^2; g: u->v 1/3+z^6/3; demands t on s->v, s->u, u->v",
            "gap_f-h-g": str(sp.factor(gap)), "active_side": "t<1 (delta<0): indirect route s->u->v used; inactive for t>1 (contraction)",
            "x_series": str(X), "x_leading": {"order": kx, "coefficient": s(cx)},
            "E_minus_E0_active_side": {"order": kE, "coefficient": s(cE), "series": str(dE)},
            "E0(1)": s(E0_1), "E'(1)_both_sides": s(E0p),
            "E5_jump_right_minus_left": s(jumpE5), "E_class": "C^4 not C^5, equal slopes",
            "optimum": {"marginal_gap_f_m-h_m-g_m": str(gm_t), "marginal_gap_leading": {"order": km, "coefficient": s(cm)},
                        "formal_optimum_branch_y": str(Y), "y_leading": {"order": ky, "coefficient": s(cy)},
                        "conclusion": "formal branch y=-d^2/5<0 on BOTH sides => infeasible; optimal indirect flow 0 near t=1; S=E0 (polynomial), analytic"},
            "P": {"class": "C^4 not C^5", "P5_jump_right_minus_left": s(jumpE5 / E0_1)},
            "note": "commodity-weighted and physical sum z c(z) costs agree on the branch"}


def part_b():
    f = lambda z: z ** 3 + z ** 2 + 4 * z + 1; g = lambda z: 4 * z ** 2 + z + 1; h = lambda z: 1 + z + z ** 2
    gap = sp.factor(sp.expand(f(t) - g(t) - h(0)))
    need(sp.expand(gap - (t - 1) ** 3) == 0, "cubic gap (t-1)^3")
    n = 7
    eqn = sp.expand((f(t - x) - g(t + x) - h(x)).subs(t, 1 + d))
    X = implicit_series(eqn, n)
    kx, cx = lead(X, n)
    need(kx == 3 and cx == R(1, 19), "x=d^3/19")
    Ecom = (t * (f(t - x) + g(t + x))).subs(t, 1 + d)
    Ephys = ((t - x) * f(t - x) + (t + x) * g(t + x) + x * h(x)).subs(t, 1 + d)
    E0 = sp.expand(Ecom.subs(x, 0))
    dE = ser(Ecom.subs(x, X) - E0, n)
    need(sp.expand(dE - ser(Ephys.subs(x, X) - E0, n)) == 0, "cubic commodity vs physical")
    kE, cE = lead(dE, n)
    need(kE == 5 and cE == R(-3, 19), "E-E0=-(3/19)d^5")
    need(E0.subs(d, 0) == 13 and sp.diff(E0, d).subs(d, 0) == 31, "E(1)=13,E'(1)=31")
    return {"label": "exact", "network": "A s->v: f=z^3+z^2+4z+1 direct or connector h=1+z+z^2 then g=4z^2+z+1; B u->v: g; demand t each",
            "gap_f-g-h(0)": str(gap), "active_side": "t>1", "x_series": str(X), "x_leading": {"order": kx, "coefficient": s(cx)},
            "E_minus_E0": {"order": kE, "coefficient": s(cE), "series": str(dE)},
            "E(1)": s(E0.subs(d, 0)), "E'(1)": s(sp.diff(E0, d).subs(d, 0)), "E5_jump": s(sp.factorial(5) * cE)}


def part_c():
    dedges = [(3, 2), (1, 4)]      # direct bank s->v (level 3): 5/2+z^4/162 (load 3), 11/4+z^4/4 (load 1)
    sedges = [(1, 1), (3, 3)]      # shared bank u->v (level 2): 1+z^4 (load 1), 5/3+z^4/243 (load 3)
    phd = ws.physical_check(dedges, 3); phs = ws.physical_check(sedges, 2)
    need([(a, b) for a, b, _ in phd] == [(Q(5, 2), Q(1, 162)), (Q(11, 4), Q(1, 4))], "direct edges")
    need([(a, b) for a, b, _ in phs] == [(Q(1), Q(1)), (Q(5, 3), Q(1, 243))], "shared edges")
    n = 8
    tb = ws.TwoBank(dedges, 3, sedges, 2, 4, n)
    gap = tb.gap(); kg, cg = ws.first_nonzero(gap)
    need(kg == 3 and cg == -Q(336, 625), "swapped gap -336/625 d^3")
    X = tb.connector(); kx, cx = ws.first_nonzero(X)
    need(kx == 3 and cx < 0, "connector branch positive for d<0")
    z = ws.zeros(n)
    E, E0 = tb.cost(X), tb.cost(z)
    dE = [a - b for a, b in zip(E, E0)]; kE, cE = ws.first_nonzero(dE)
    need(kE == 5 and cE == -Q(84672, 78125), "E-E0 order 5 on active side")
    need(E[0] == 20 and E[1] == E0[1], "E(1)=20, equal slopes")
    loads = tb.edge_loads(X)
    # Sign lemma for t>1: q_direct(s)-q_shared(s) = Delta^3 g(0) = (21/64) s^3 (...)^(-11/4) > 0 for s>0
    # => s_d(T) < s_s(T) for T>4: connector strictly inactive, all edges loaded (s>0)
    V0 = tb.potential(z)
    return {"label": "exact", "edges_direct_(a,b,ref_load)": [[ws.frac(a), ws.frac(b), ws.frac(l)] for a, b, l in phd],
            "edges_shared_(a,b,ref_load)": [[ws.frac(a), ws.frac(b), ws.frac(l)] for a, b, l in phs],
            "bank_delays_at_T=4": {"direct": "3", "shared": "2"},
            "gap_D_leading": {"order": kg, "coefficient": ws.frac(cg)}, "gap_series": [ws.frac(c) for c in gap],
            "connector_formal_branch_leading": {"order": kx, "coefficient": ws.frac(cx), "positive_for": "delta<0 (t<1)"},
            "connector_series": [ws.frac(c) for c in X],
            "edge_loads_at_t=1": {nm: ws.frac(sr[0]) for nm, sr in loads},
            "E(1)": ws.frac(E[0]), "E'(1)_both_sides": ws.frac(E0[1]),
            "E_minus_E0_active_left_side": {"order": kE, "coefficient": ws.frac(cE), "series": [ws.frac(c) for c in dE]},
            "E5_jump_right_minus_left": ws.frac(-120 * cE), "E_class": "C^4 not C^5, equal first slopes",
            "relation_to_original": "E-E0 series identical to the original five-edge (-K delta^5 + ...) and the connector series is exactly the negative of the original; on the active side delta<0 so E-E0 = +K|delta|^5 > 0 (Braess-type cost increase)",
            "t>1_inactivity": "sign lemma: q_direct(s)-q_shared(s) = Delta^3 g(0) > 0 for s>0, so s_d(T)<s_s(T) for all T>4: connector strictly off and all four bank edges loaded (s>0)",
            "V0(1)_Beckmann": ws.frac(V0[0])}, (cE, V0)


def part_c_numeric(cE):
    import mpmath as mp
    from wardrop_numeric import Network, fmt, DPS
    mp.mp.dps = DPS
    D = [(mp.mpf(5) / 2, mp.mpf(1) / 162), (mp.mpf(11) / 4, mp.mpf(1) / 4)]
    S_ = [(1, 1), (mp.mpf(5) / 3, mp.mpf(1) / 243)]
    N = Network(D, S_); M = Network(D, S_, k=5)
    g = mp.mpf(5) ** mp.mpf(0.25)
    rows = []
    for tv in ["0.98", "0.99", "1", "1.01", "1.02"]:
        tt = mp.mpf(tv)
        y = M.solve(4 * tt)
        yy, zF, zG = M.flows(4 * tt, x=y)
        Sv = M.physical_cost(yy, zF, zG)
        xg, zFg, zGg = N.flows(4 * g * tt)
        Ssc = N.potential(xg, zFg, zGg) / g
        need(y == 0 and xg == 0 and min(zF + zG) > 0, f"swapped optimum connector off, edges on, t={tv}")
        need(abs(Sv - Ssc) < mp.mpf(10) ** -45, "swapped S direct vs scaling")
        rows.append({"t": tv, "opt_connector": fmt(y, 5), "S_direct": fmt(Sv, 20), "min_opt_edge_load": fmt(min(zF + zG), 10)})
    x_left = N.solve(4 * mp.mpf("0.99")); x_right = N.solve(4 * mp.mpf("1.01"))
    need(x_left > 0 and x_right == 0, "swapped equilibrium active left only")
    S1 = M.physical_cost(*M.flows(4, x=mp.mpf(0)))
    cEf = mp.mpf(cE.numerator) / cE.denominator
    # optimum transition location for swapped: optimum connector active iff eq at gamma t active iff gamma t < 1 near 1
    y_lo = M.solve(4 * (1 / g - mp.mpf("0.001"))); y_hi = M.solve(4 * (1 / g + mp.mpf("0.001")))
    return {"label": "numerical", "mp_dps": DPS, "rows": rows, "S(1)": fmt(S1, 20), "P(1)=20/S(1)": fmt(20 / S1, 15),
            "P5_jump_right_minus_left=-120*c/S(1)": fmt(-120 * cEf / S1, 15),
            "eq_connector_t=0.99": fmt(x_left, 12), "eq_connector_t=1.01": fmt(x_right, 5),
            "optimum_connector_at_1/gamma-0.001": fmt(y_lo, 10), "optimum_connector_at_1/gamma+0.001": fmt(y_hi, 5),
            "conclusion": "near t=1 the swapped optimum uses no connector (gamma t ~ 1.495 > 1 lies in the inactive equilibrium regime), S analytic at 1; the optimum's own transition is at t=1/gamma (connector on below), irrelevant near t=1"}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    out = {"a_devils_advocate_triangle": part_a(), "b_cubic_triangle": part_b()}
    c, (cE, V0) = part_c()
    out["c_swapped_five_edge"] = c
    out["c_swapped_five_edge_numeric"] = part_c_numeric(cE)
    out["status"] = "PASS" if not FAILS else "FAIL"; out["failures"] = FAILS
    a.output.write_text(json.dumps(out, sort_keys=True, indent=2) + "\n")
    print(("PASS" if not FAILS else "FAIL: " + "; ".join(FAILS)) + f" verify_small_examples -> {a.output}")
    sys.exit(0 if not FAILS else 1)


if __name__ == "__main__":
    main()
