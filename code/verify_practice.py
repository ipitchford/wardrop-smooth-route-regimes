#!/usr/bin/env python3
"""Practical-magnitude verifier for the five-edge quartic example.

Records labelled "numerical" use mpmath (dps stated); records labelled
"exact" use fractions.Fraction series from wardrop_series.py.
"""
import argparse, json, sys
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mpmath as mp
from wardrop_numeric import Network, Bank, FIVE_DIRECT, FIVE_SHARED, fmt, DPS
import wardrop_series as ws

FAILS = []


def need(cond, msg):
    if not cond:
        FAILS.append(msg)
    return bool(cond)


def num(rec):
    rec["label"] = "numerical"; rec["mp_dps"] = DPS
    return rec


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    out = {}
    N = Network(FIVE_DIRECT, FIVE_SHARED)
    E = lambda t, eps=0: N.eq_cost_formula(t * 4, eps=eps)
    E0 = lambda t: N.eq_cost_formula(4 * t, x=mp.mpf(0))
    K = mp.mpf(84672) / 78125
    # exact series (reference coefficients)
    F, G = ws.binomial_family_banks(3)
    tb = ws.TwoBank(F, 3, G, 2, 4, 8)
    xs = tb.connector(); z = ws.zeros(8)
    Es, E0s = tb.cost(xs), tb.cost(z)
    ws.check(xs[3] == Q(84, 125) and Es[5] - E0s[5] == -Q(84672, 78125) and Es[0] == 20 and Es[1] == Q(164, 5), "exact base")
    out["exact_reference"] = {"label": "exact", "x_series_delta^3..^5": [ws.frac(c) for c in xs[3:6]],
                              "E_minus_E0_series_delta^5..^6": [ws.frac(Es[k] - E0s[k]) for k in (5, 6)],
                              "E(1)": ws.frac(Es[0]), "E'(1)": ws.frac(Es[1]), "gap_a": ws.frac(tb.gap()[3]),
                              "beta=q_F'(0)": ws.frac(tb.beta_d)}
    # (a) table
    rows = []
    for ds in ["0.01", "0.05", "0.10", "0.20"]:
        d = mp.mpf(ds); t = 1 + d
        x = N.solve(4 * t); e, e0 = E(t), E0(t)
        # equilibrium residual check
        res = N.F.delay(4 * t - x) - N.G.delay(4 * t + x) - 1 - x ** 4
        need(abs(res) < mp.mpf(10) ** -45 and x > 0, f"eq residual {ds}")
        rows.append({"delta": ds, "x": fmt(x, 15), "x_percent_of_A_demand": fmt(100 * x / (4 * t), 12),
                     "rel_cost_change_(E-E0)/E0": fmt((e - e0) / e0, 12), "E_minus_E0": fmt(e - e0, 15),
                     "E": fmt(e, 20), "E0": fmt(e0, 20), "x_over_delta3": fmt(x / d ** 3, 12),
                     "(E-E0)_over_delta5": fmt((e - e0) / d ** 5, 12)})
        if ds == "0.01":
            need(abs(x / d ** 3 - mp.mpf(84) / 125) < 0.05, "x/d^3 near 84/125")
            need(abs((e - e0) / d ** 5 + K) < 0.15, "(E-E0)/d^5 near -K")
    out["a_activation_table"] = num({"rows": rows})
    # (b) hidden activation vs tolerance
    def largest(fn, tau):
        lo, hi = mp.mpf("1e-6"), mp.mpf(4)
        need(fn(hi) > tau, "bracket")
        for _ in range(60):
            mid = (lo + hi) / 2
            if fn(mid) < tau: lo = mid
            else: hi = mid
        return lo
    fx = lambda d: N.solve(4 * (1 + d)) / (4 * (1 + d))
    fe = lambda d: abs(E(1 + d) - E0(1 + d)) / E0(1 + d)
    grid = [mp.mpf(k) / 50 for k in range(1, 201)]
    vx = [fx(d) for d in grid]; ve = [fe(d) for d in grid]
    need(all(vx[i] < vx[i + 1] for i in range(199)) and all(ve[i] < ve[i + 1] for i in range(199)),
         "monotone ratios on grid delta in [0.02,4]")
    thr = []
    for tau in ["1e-4", "1e-5", "1e-6"]:
        T = mp.mpf(tau)
        dx, de = largest(fx, T), largest(fe, T)
        thr.append({"tau": tau, "largest_delta_x_over_4t_below_tau": fmt(dx, 10),
                    "largest_delta_relcost_below_tau": fmt(de, 10),
                    "leading_order_estimate_x": fmt((T * 4 * 125 / 84) ** (mp.mpf(1) / 3), 10),
                    "leading_order_estimate_cost": fmt((20 * T / K) ** (mp.mpf(1) / 5), 10)})
    out["b_hidden_activation_vs_tolerance"] = num({"monotonicity_checked_on_grid": "delta=0.02..4.00 step 0.02",
                                                   "note": "bisection, 60 halvings of [1e-6,4]", "rows": thr})
    # (c) S(1): direct physical optimisation (KKT with marginal costs a+5bz^4) vs scaling
    M = Network(FIVE_DIRECT, FIVE_SHARED, k=5)
    y = M.solve(4)
    yy, zF, zG = M.flows(4, x=y)
    S_direct = M.physical_cost(yy, zF, zG)
    mF, mG = M.F.delay(4 - y), M.G.delay(4 + y)
    kkt = mF - mG - (1 + 5 * y ** 4)
    need(abs(kkt) < mp.mpf(10) ** -45 and y > 0, "KKT optimum")
    gam = mp.mpf(5) ** mp.mpf(0.25)
    xg, zFg, zGg = N.flows(4 * gam)
    V = N.potential(xg, zFg, zGg); S_scal = V / gam
    need(abs(S_direct - S_scal) < mp.mpf(10) ** -45, "S direct vs scaling")
    need(abs(S_direct - mp.mpf("19.762875963175")) < 1e-11, "S(1) expected value")
    zmatch = max(abs(a - b / gam) for a, b in zip([y] + zF + zG, [xg] + zFg + zGg))
    need(zmatch < mp.mpf(10) ** -45, "z_opt = z_eq(gamma t)/gamma")
    # scipy SLSQP brute force cross-check (double precision)
    import numpy as np
    from scipy.optimize import minimize
    a = np.array([2, 8 / 3, 3 / 2, 7 / 4, 1.]); b = np.array([1, 1 / 243, 1 / 162, 1 / 4, 1.])
    def loads(f):  # A paths: F0, F2, h-G1, h-G3 ; B paths: G1, G3
        return np.array([f[0], f[1], f[2] + f[4], f[3] + f[5], f[2] + f[3]])
    obj = lambda f: float(np.sum(loads(f) * (a + b * loads(f) ** 4)))
    cons = [{"type": "eq", "fun": lambda f: f[:4].sum() - 4}, {"type": "eq", "fun": lambda f: f[4:].sum() - 4}]
    best = None
    for seed in range(8):
        r0 = np.random.default_rng(seed)
        f0 = np.r_[r0.dirichlet(np.ones(4)) * 4, r0.dirichlet(np.ones(2)) * 4]
        r = minimize(obj, f0, constraints=cons, bounds=[(0, None)] * 6, method="SLSQP", options={"ftol": 1e-15, "maxiter": 3000})
        if r.success and (best is None or r.fun < best.fun): best = r
    need(abs(best.fun - float(S_direct)) < 1e-8, "SLSQP agrees")
    Pt = 20 / S_direct
    out["c_social_optimum_t1"] = num({
        "method_direct": "KKT of min sum z c(z): bank marginal costs a+5bz^4 equalised, connector marginal 1+5y^4; convex so KKT sufficient",
        "S1_direct": fmt(S_direct, 25), "S1_scaling_V(gamma)/gamma": fmt(S_scal, 25),
        "abs_diff": fmt(abs(S_direct - S_scal), 3), "S1_scipy_SLSQP_double": f"{best.fun:.12f}",
        "P1=E(1)/S(1)": fmt(Pt, 15), "optimal_connector_load_t1": fmt(y, 15),
        "optimal_loads_F0_F2_G1_G3": [fmt(v, 12) for v in zF + zG],
        "KKT_residual": fmt(abs(kkt), 3), "max|z_opt - z_eq(gamma)/gamma|": fmt(zmatch, 3),
        "path_marginal_costs": {"A_direct": fmt(mF, 15), "A_via_connector": fmt(1 + 5 * y ** 4 + mG, 15), "B": fmt(mG, 15)}})
    need(abs(y - mp.mpf("0.0064")) < 0.0005, "opt connector ~0.0064")
    # (d) P^(5) jump
    J5 = -120 * K / S_direct
    out["d_P5_jump_t1"] = num({"formula": "P^(5)(1+)-P^(5)(1-) = -5! K / S(1) (S analytic at 1 since optimum connector active on both sides; E jump -5!K)",
                               "value": fmt(J5, 15), "E5_jump_exact": ws.frac(-120 * Q(84672, 78125))})
    need(abs(J5 + mp.mpf("6.58")) < 0.01, "P5 jump ~ -6.58")
    # (e) near-match kink
    a3 = mp.mpf(336) / 625; beta = mp.mpf(5) / 2; m = 3
    Dfun = lambda t: N.F.delay(4 * t) - N.G.delay(4 * t) - 1
    Dp = lambda t: 4 * (N.F.ddelay(4 * t) - N.G.ddelay(4 * t))
    nm = []
    def t_eps_near1(eps):
        """activation point of the branch through t=1: for eps>0 first t>1 with D=eps;
        for eps<0 last t<1 with D<=eps (scan step 1e-3, then bisection)."""
        step = mp.mpf("0.001") * (1 if eps > 0 else -1)
        t = mp.mpf(1)
        while True:
            tn = t + step
            if tn <= mp.mpf("0.05") or tn > 20:
                return None
            if (eps > 0 and Dfun(tn) >= eps) or (eps < 0 and Dfun(tn) <= eps):
                lo, hi = (t, tn) if eps > 0 else (tn, t)
                break
            t = tn
        for _ in range(200):
            mid = (lo + hi) / 2
            if Dfun(mid) < eps: lo = mid
            else: hi = mid
        return (lo + hi) / 2
    for es in ["1e-9", "-1e-9", "1e-6", "-1e-6", "1e-4", "-1e-4", "0.01", "-0.01", "0.1", "-0.1"]:
        eps = mp.mpf(es)
        te = t_eps_near1(eps)
        Fp, Gp = N.F.ddelay(4 * te), N.G.ddelay(4 * te)
        JE = te * Dp(te) ** 2 / (Fp + Gp)
        asym = beta / 2 * m ** 2 * a3 ** (mp.mpf(2) / m) * abs(eps) ** (2 - mp.mpf(2) / m)
        h = mp.mpf("1e-12")
        Ee = lambda t: E(t, eps)
        right = (-3 * Ee(te) + 4 * Ee(te + h) - Ee(te + h * 2)) / (h * 2)
        left = (3 * Ee(te) - 4 * Ee(te - h) + Ee(te - h * 2)) / (h * 2)
        fd = right - left
        need(abs(fd + JE) < mp.mpf(10) ** -9 * max(1, JE), f"FD jump matches formula eps={es}")
        if abs(eps) <= mp.mpf("1e-9"):
            need(abs(JE / asym - 1) < 0.02, f"asymptotic limit eps={es}")
        nF = sum(1 for v in N.F.loads_at(N.F.delay(4 * te)) if v > 0)
        nG = sum(1 for v in N.G.loads_at(N.G.delay(4 * te)) if v > 0)
        nm.append({"eps": es, "eps_percent_of_connector_freeflow": fmt(100 * eps, 6), "t_eps": fmt(te, 15),
                   "J_E_formula": fmt(JE, 15), "J_E_over_E'(1)=32.8": fmt(JE / mp.mpf("32.8"), 12),
                   "asymptotic_(beta/2)m^2a^(2/m)|eps|^(2-2/m)": fmt(asym, 15),
                   "formula_over_asymptotic": fmt(JE / asym, 12),
                   "finite_difference_slope_jump_E'(t+)-E'(t-)": fmt(fd, 15), "fd_h": "1e-12 (second-order one-sided)",
                   "active_bank_edges_at_t_eps": {"direct": nF, "shared": nG}})
    # minimum of D on the inactive window below 1 (eps below this: no near-1 kink)
    tg3 = (3 * mp.mpf(2) ** mp.mpf(-0.25)) / 4
    Dmin_t = min([mp.mpf(k) / 2000 for k in range(1170, 2000)] + [tg3], key=Dfun)
    rD = [mp.findroot(Dfun, (mp.mpf("0.2"), mp.mpf("0.22")), solver="anderson"),
          mp.findroot(Dfun, (mp.mpf("0.55"), mp.mpf("0.62")), solver="anderson")]
    out["global_connector_activity_scan"] = num({
        "finding": "equilibrium connector is ALSO active for t in (r1, r2) (low-demand window): free-flow 1+3/2 < F-bank delay once F0 congests; x=0 only on [0,r1] and [r2,1]",
        "r1": fmt(rD[0], 20), "r2": fmt(rD[1], 20),
        "D_min_on_[r2,1]": fmt(Dfun(Dmin_t), 15), "argmin_t": fmt(Dmin_t, 12),
        "G3_activation_t=3*2^(-1/4)/4": fmt(tg3, 15),
        "sample_t0.3_x_over_4t": fmt(N.solve(mp.mpf("1.2")) / mp.mpf("1.2"), 10),
        "consequence_for_eps<0": "for eps < D_min the near-1 activation point does not exist; the eps=-0.1 row is the low-demand activation at t~0.199, not a near-match kink"})
    out["e_near_match_kink"] = num({"definition": "connector cost 1+eps+x^4; J_E = t D'(t)^2/(F'(4t)+G'(4t)) at D(t_eps)=eps; E' drops by J_E",
                                    "rows": nm})
    # (f) BPR
    phys = [("F0", 2, 1, 1), ("F2", mp.mpf(8) / 3, mp.mpf(1) / 243, 3), ("G1", mp.mpf(3) / 2, mp.mpf(1) / 162, 3),
            ("G3", mp.mpf(7) / 4, mp.mpf(1) / 4, 1), ("h", 1, 1, 0)]
    bpr = []
    for nme, aa, bb, ld in phys:
        C = (mp.mpf("0.15") * aa / bb) ** mp.mpf(0.25)
        bpr.append({"edge": nme, "capacity_C": fmt(C, 10), "volume_over_capacity_t1": fmt(ld / C, 10),
                    "congested_over_freeflow_t1": fmt((aa + bb * mp.mpf(ld) ** 4) / aa, 10)})
    out["f_BPR"] = num({"convention": "c = a(1+0.15 (z/C)^4), C=(0.15 a/b)^(1/4)", "rows": bpr})
    # (g) per-commodity costs and (h) link loads
    A, B = tb.commodity_costs(xs); A0, B0 = tb.commodity_costs(z)
    dA, dB = ws.add(A, ws.scale(A0, -1, 8), 8), ws.add(B, ws.scale(B0, -1, 8), 8)
    kA, cA = ws.first_nonzero(dA); kB, cB = ws.first_nonzero(dB)
    d = mp.mpf("1e-3"); t = 1 + d; x = N.solve(4 * t)
    numA = (N.F.delay(4 * t - x) - N.F.delay(4 * t)) / d ** 3
    numB = (N.G.delay(4 * t + x) - N.G.delay(4 * t)) / d ** 3
    need(kA == 3 and kB == 3 and abs(numA - mp.mpf(cA.numerator) / cA.denominator) < 0.01, "commodity order 3")
    out["g_per_commodity_costs"] = {"label": "exact+numerical", "mp_dps": DPS,
        "A_cost_minus_inactive_continuation_leading": {"order": kA, "coefficient_exact": ws.frac(cA), "numerical_ratio_delta=1e-3": fmt(numA, 10)},
        "B_cost_minus_inactive_continuation_leading": {"order": kB, "coefficient_exact": ws.frac(cB), "numerical_ratio_delta=1e-3": fmt(numB, 10)},
        "A_minus_B_cost_series_delta^0..3": [ws.frac(a_ - b_) for a_, b_ in zip(A[:4], B[:4])],
        "smoothness": "per-commodity costs are C^2 but not C^3 at t=1 (jump in third derivative 6*coef)"}
    loads_rec = {"connector": {"order": 3, "coefficient_exact": ws.frac(xs[3]), "class": "C^2 not C^3"}}
    el, el0 = tb.edge_loads(xs), tb.edge_loads(z)
    names = {"direct_w1_u1": "F0", "direct_w3_u3": "F2", "shared_w3_u2": "G1", "shared_w1_u4": "G3"}
    xF, zFn, zGn = N.flows(4 * t, x=x); _, zF0, zG0 = N.flows(4 * t, x=mp.mpf(0))
    numl = dict(zip(["F0", "F2", "G1", "G3"], [(p_ - q_) / d ** 3 for p_, q_ in zip(zFn + zGn, zF0 + zG0)]))
    for (nmx, s1), (_, s0) in zip(el, el0):
        k_, c_ = ws.first_nonzero(ws.add(s1, ws.scale(s0, -1, 8), 8))
        nn = names[nmx]
        loads_rec[nn] = {"order": k_, "coefficient_exact": ws.frac(c_), "numerical_ratio_delta=1e-3": fmt(numl[nn], 10), "class": "C^2 not C^3"}
        need(k_ == 3, "edge order 3")
    out["h_link_flows_vs_aggregate"] = {"label": "exact+numerical", "mp_dps": DPS, "loads": loads_rec,
        "aggregate": "E-E0 = -K delta^5: E is C^4 not C^5 while every link load and both commodity costs are only C^2"}
    # E'(1) numerical
    hh = mp.mpf("1e-15")
    out["E_derivative_check"] = num({"E(1)": fmt(E(1), 20), "E'(1)_central_diff": fmt((E(1 + hh) - E(1 - hh)) / (hh * 2), 15), "expected": "164/5"})
    need(abs((E(1 + hh) - E(1 - hh)) / (hh * 2) - mp.mpf(164) / 5) < 1e-10, "E'(1)")
    out["status"] = "PASS" if not FAILS else "FAIL"; out["failures"] = FAILS
    args.output.write_text(json.dumps(out, sort_keys=True, indent=2) + "\n")
    print(("PASS" if not FAILS else "FAIL: " + "; ".join(FAILS)) + f" verify_practice -> {args.output}")
    sys.exit(0 if not FAILS else 1)


if __name__ == "__main__":
    main()
