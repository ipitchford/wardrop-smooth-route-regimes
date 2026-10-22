#!/usr/bin/env python3
"""Optimum transition at t = 1/gamma for the five-edge example (gamma = 5^(1/4)).

(i)   exact: Beckmann value V along the equilibrium ray, V'(T) = E(T)/T as a series
      identity, V - V0 = -(K/6) delta^6 + ...; same for the binomial family m=5 (p=4).
(ii)  exact rational interval proof that for t in [1/gamma - 1/100, 1) the
      equilibrium connector is strictly inactive and all four bank edges are loaded.
(iii) numerical values at 1/gamma and a cross-check of S - S0 near 1/gamma.
"""
import argparse, json, sys
from fractions import Fraction as Q
from math import isqrt, factorial
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import wardrop_series as ws
import mpmath as mp
from wardrop_numeric import Network, FIVE_DIRECT, FIVE_SHARED, fmt, DPS

FAILS = []


def need(c, m):
    if not c:
        FAILS.append(m)
    return bool(c)


def series_part(m, n):
    F, G = ws.binomial_family_banks(m)
    D = 2 ** (m - 1)
    tb = ws.TwoBank(F, 3, G, 2, D, n)
    x = tb.connector(); z = ws.zeros(n)
    E, E0 = tb.cost(x), tb.cost(z)
    V, V0 = tb.potential(x), tb.potential(z)
    # envelope identity  (1+delta) V'(delta) = E(delta)  (T = 1+delta is the demand multiplier)
    for VV, EE, nm in ((V, E, "active"), (V0, E0, "inactive")):
        lhs = ws.mul([Q(1), Q(1)], ws.deriv(VV), n)
        need(all(lhs[k] == EE[k] for k in range(n)), f"envelope identity m={m} {nm}")
    dE = [a - b for a, b in zip(E, E0)]; dV = [a - b for a, b in zip(V, V0)]
    kE, cE = ws.first_nonzero(dE); kV, cV = ws.first_nonzero(dV)
    K = -cE
    need(kE == 2 * m - 1 and kV == 2 * m, f"orders m={m}")
    need(cV == -K / (2 * m), f"V-V0 leading = -K/(2m) m={m}")
    return {"label": "exact", "m": m, "p": 4, "demand_each": D, "series_order": n,
            "envelope_identity_(1+delta)V'=E": "verified coefficientwise through delta^%d (active and inactive branches)" % (n - 1),
            "V(1)": ws.frac(V[0]), "E(1)": ws.frac(E[0]),
            "E_minus_E0_leading": {"order": kE, "coefficient": ws.frac(cE)}, "K": ws.frac(K),
            "V_minus_V0_leading": {"order": kV, "coefficient": ws.frac(cV), "equals_-K/(2m)": ws.frac(-K / (2 * m))},
            "V_minus_V0_next": ws.frac(dV[kV + 1]) if kV + 1 <= n else None,
            "x_leading": ws.frac(x[m])}, K, V[0]


# ---------- exact rational fourth-root bounds ----------
PREC = 60


def lo4(q):
    q = Q(q); ws.check(q >= 0, "neg")
    N = (q.numerator << (4 * PREC)) // q.denominator
    r = isqrt(isqrt(N))
    while (r + 1) ** 4 <= N: r += 1
    while r ** 4 > N: r -= 1
    out = Q(r, 1 << PREC)
    ws.check(out ** 4 <= q, "lo4")
    return out


def hi4(q):
    out = lo4(q) + Q(1, 1 << PREC)
    ws.check(out ** 4 >= q, "hi4")
    return out


def qF_hi(s):
    return hi4(1 + s) + 3 * hi4(1 + 3 * s)


def proof_part():
    r = lo4(Q(1, 5))                     # r <= 5^(-1/4) = 1/gamma
    t_lo = r - Q(1, 100)
    T_lo = 4 * t_lo
    steps = []
    # G3 loaded: s_G > -1/4  <=> T > q_G(-1/4) = 3*2^(-1/4) <=> T^4 > 81/2
    need(T_lo ** 4 > Q(81, 2), "G3 positive"); steps.append("T_lo^4 > 81/2 => q_G(-1/4)=3*2^(-1/4) < T_lo => s_G(T) > -1/4: G3 (and G1, threshold -1/2) loaded")
    need(T_lo ** 4 > Q(2, 3), "F2 positive"); steps.append("T_lo^4 > 2/3 => q_F(-1/3)=(2/3)^(1/4) < T_lo => s_F(T) > -1/3: F2 (and F0, threshold -1) loaded")
    # sigma: rational lower bound of s_F(T_lo): largest grid point with upper bound q_F(sigma) <= T_lo
    lo, hi = Q(-1, 3), Q(0)
    for _ in range(40):
        mid = (lo + hi) / 2
        if qF_hi(mid) <= T_lo: lo = mid
        else: hi = mid
    sigma = lo
    need(qF_hi(sigma) <= T_lo and sigma > Q(-1, 3), "sigma bound")
    steps.append("q_F(sigma) <= T_lo (upper root bounds) => s_F(T) >= sigma for T >= T_lo")
    # interval check of h(s) = q_F(s) - q_G(s) on [sigma, -1/4] where G3 term is 0
    Nsub = 400
    ok = True; minlb = None
    for i in range(Nsub):
        l = sigma + (Q(-1, 4) - sigma) * i / Nsub
        rr = sigma + (Q(-1, 4) - sigma) * (i + 1) / Nsub
        lb = lo4(1 + l) + 3 * lo4(1 + 3 * l) - 3 * hi4(1 + 2 * rr)   # G3 load max(0,1+4s)^(1/4) = 0 on this interval
        ok &= lb > 0
        minlb = lb if minlb is None or lb < minlb else minlb
    need(ok, "interval positivity of h on [sigma,-1/4]")
    steps.append(f"h(s)=q_F(s)-q_G(s) > 0 on [sigma,-1/4]: {Nsub} subintervals, monotone bounds h >= q_F_lo(l) - q_G_hi(r) (G3 inactive there)")
    steps.append("h(s) > 0 on (-1/4,0): h = -Delta^3 g(0) with g(y)=(1+(y+1)s)^(1/4); Delta^3 g(0)=g'''(xi)=(21/64) s^3 (1+(xi+1)s)^(-11/4) < 0 for s<0")
    steps.append("s_F(T) in [sigma,0) for T in [T_lo,4) and h(s_F)>0 means q_G(s_F) < T = q_G(s_G), so s_F < s_G, i.e. F(T) < G(T)+1+0^4: connector strictly unused; unique Wardrop costs => E(t)=E0(t)=4t[5+s_F(4t)+s_G(4t)] on [t_lo,1), analytic (inverse function theorem, q' > 0 with all four edges active)")
    return {"label": "exact", "method": "rational interval arithmetic (Fractions, fourth-root bounds with 2^-60 resolution) + finite-difference sign lemma",
            "t_lo_rational_lower_bound_of_1/gamma-1/100": f"{float(t_lo):.15f}", "T_lo=4 t_lo": f"{float(T_lo):.15f}",
            "q_G(-1/4)=3*2^(-1/4)": f"{3 * 2 ** -0.25:.15f}", "sigma_lower_bound_s_F(T_lo)": f"{float(sigma):.15f}",
            "min_interval_lower_bound_of_h": f"{float(minlb):.12f}", "subintervals": Nsub, "steps": steps,
            "conclusion": "for all t in [1/gamma-1/100, 1): x_eq(t)=0 strictly and F0,F2,G1,G3 all carry positive load; E analytic there (incl. a neighbourhood of 1/gamma)",
            "caveat": "the connector IS active in the separate low-demand window t in (0.2105, 0.5839) (see verify_practice.json global scan); this lies below 1/gamma-0.01=0.6587"}


def numeric_part(K, V1):
    mp.mp.dps = DPS
    N = Network(FIVE_DIRECT, FIVE_SHARED)
    M = Network(FIVE_DIRECT, FIVE_SHARED, k=5)
    Kf = mp.mpf(K.numerator) / K.denominator
    g = mp.mpf(5) ** mp.mpf(0.25); tg = 1 / g
    E = lambda t: N.eq_cost_formula(4 * t)

    def S_direct(t, conn=True):
        if conn:
            y = M.solve(4 * t)
        else:
            y = mp.mpf(0)
        yy, zF, zG = M.flows(4 * t, x=y)
        return M.physical_cost(yy, zF, zG), y

    def S_scale(t, conn=True):
        T = g * t
        x = N.solve(4 * T) if conn else mp.mpf(0)
        x, zF, zG = N.flows(4 * T, x=x)
        return N.potential(x, zF, zG) / g

    Eg = E(tg); Sg, yg = S_direct(tg); Sgs = S_scale(tg)
    V1f = mp.mpf(V1.numerator) / V1.denominator
    need(abs(Sg - V1f / g) < mp.mpf(10) ** -45 and abs(Sgs - Sg) < mp.mpf(10) ** -45, "S(1/gamma)=V(1)/gamma")
    need(yg == 0, "optimum connector off at 1/gamma")
    Pg = Eg / Sg
    J6 = Eg * 120 * Kf * g ** 5 / Sg ** 2
    rows = []
    for ds in ["1e-2", "1e-3"]:
        d = mp.mpf(ds); t = tg + d
        Sa, ya = S_direct(t); S0, _ = S_direct(t, conn=False)
        Sa2, S02 = S_scale(t), S_scale(t, conn=False)
        need(abs(Sa - Sa2) < mp.mpf(10) ** -45 and abs(S0 - S02) < mp.mpf(10) ** -45, "direct vs scaling near 1/gamma")
        pred = -Kf * g ** 5 / 6 * d ** 6
        Pdiff = E(t) / Sa - E(t) / S0
        Ppred = E(t) * Kf * g ** 5 / (6 * S0 ** 2) * d ** 6
        rows.append({"delta=t-1/gamma": ds, "S_minus_S0_direct_KKT": fmt(Sa - S0, 15), "prediction_-(K gamma^5/6) delta^6": fmt(pred, 15),
                     "ratio": fmt((Sa - S0) / pred, 12), "optimal_connector_load": fmt(ya, 12),
                     "P_minus_P0": fmt(Pdiff, 15), "P_ratio_to_E K gamma^5 delta^6/(6 S0^2)": fmt(Pdiff / Ppred, 12)})
        need(abs((Sa - S0) / pred - 1) < (0.2 if ds == "1e-2" else 0.03), f"S expansion ratio {ds}")
        # left side: optimum connector off
        Sl, yl = S_direct(tg - d)
        need(yl == 0, "optimum connector off left of 1/gamma")
    return {"label": "numerical", "mp_dps": DPS, "gamma": fmt(g, 25), "1/gamma": fmt(tg, 25),
            "E(1/gamma)": fmt(Eg, 25), "S(1/gamma)": fmt(Sg, 25), "S(1/gamma)_exact_form": f"V(1)/gamma = {ws.frac(V1)}/5^(1/4)",
            "P(1/gamma)": fmt(Pg, 20), "P6_jump_+E*5!*K*gamma^5/S^2": fmt(J6, 20),
            "S6_jump_-(K gamma^5/6)*720": fmt(-120 * Kf * g ** 5, 20),
            "S_expansion_crosscheck": rows,
            "note": "S computed two ways: direct KKT optimisation of sum z c(z) and V(gamma t)/gamma; S0 = connector forbidden"}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    r3, K3, V1 = series_part(3, 9)
    need(K3 == Q(84672, 78125) and V1 == 17, "m=3 K and V(1)")
    r5, K5, _ = series_part(5, 12)
    out = {"i_series_m3": r3, "i_series_m5": r5,
           "claims": {"S(t)=V(gamma t)/gamma": "checked numerically (direct KKT vs scaling) in iii",
                      "S-S0 = -(K gamma^5/6)(t-1/gamma)^6": "follows from V-V0=-(K/6)(T-1)^6 (exact, i) with T=gamma t; checked numerically in iii",
                      "S in C^5 not C^6 at 1/gamma; P in C^5 not C^6 with P6 jump +E 5! K gamma^5/S^2": "requires E analytic near 1/gamma: proved in ii"},
           "ii_inactive_interval_proof": proof_part(), "iii_numerics": numeric_part(K3, V1)}
    out["status"] = "PASS" if not FAILS else "FAIL"; out["failures"] = FAILS
    a.output.write_text(json.dumps(out, sort_keys=True, indent=2) + "\n")
    print(("PASS" if not FAILS else "FAIL: " + "; ".join(FAILS)) + f" verify_optimum_transition -> {a.output}")
    sys.exit(0 if not FAILS else 1)


if __name__ == "__main__":
    main()
