"""Semantic negative controls for the five-link quartic example (Section 5.2).

Standard library only (exact rational arithmetic). Run from the package root:
    python3 -B -O negative_controls.py

The check recomputes, from the physical link costs a + b z^4, the bank-delay jets at bank
load T = 4 (demand scale t = 1): F(4) = 3, G(4) = 2, F' = G' = 2/5, F'' = G'' = 42/125,
F''' - G''' = 63/1250, hence the route gap Delta(t) = F(4t) - G(4t) - 1 = (336/625)(t - 1)^3 + ...
It then confirms that plausible corruptions are REJECTED: a perturbed connector free-flow time,
perturbed bank coefficients or free-flow times (the published loads then
fail to be a Wardrop split or the contact breaks) and a wrong claimed gap
coefficient. Assertions are explicit (not `assert`), so the controls also run under -O.
Exit status is nonzero on any failure.
"""
from fractions import Fraction as Q
import sys


def binom_quarter(k):
    """Generalised binomial coefficient C(1/4, k)."""
    c = Q(1)
    for j in range(k):
        c = c * (Q(1, 4) - j) / (j + 1)
    return c


def bank_jets(links, loads, order=3):
    """Jets F(T0), F'(T0), ..., for a parallel bank of links (a, b) with cost a + b z^4,
    all used, at the Wardrop split with the given rational link loads (equal costs required)."""
    levels = {a + b * z ** 4 for (a, b), z in zip(links, loads)}
    if len(levels) != 1:
        return None  # loads are not a Wardrop split: costs differ
    c0 = levels.pop()
    # T(c0 + s) - T0 = sum_k tk s^k  with z_i(s) = z_i (1 + s/(b_i z_i^4))^(1/4)
    t = [Q(0)] * (order + 1)
    for (a, b), z in zip(links, loads):
        u = 1 / (b * z ** 4)
        for k in range(1, order + 1):
            t[k] += z * binom_quarter(k) * u ** k
    a1, a2, a3 = t[1], t[2], t[3]
    s1, s2, s3 = 1 / a1, -a2 / a1 ** 3, (2 * a2 ** 2 - a1 * a3) / a1 ** 5  # series reversion
    return [c0, s1, 2 * s2, 6 * s3]  # F, F', F'', F'''


def gap_coefficient(F, G, connector_a):
    """Return (gap at t=1, leading cubic coefficient of the gap in delta = t - 1) if the
    first two derivatives agree; otherwise (gap, None)."""
    gap0 = F[0] - G[0] - connector_a
    if F[1] != G[1] or F[2] != G[2]:
        return gap0, None
    return gap0, (F[3] - G[3]) * 4 ** 3 / 6


def evaluate(params):
    F = bank_jets(params["F"], [Q(1), Q(3)])
    G = bank_jets(params["G"], [Q(3), Q(1)])
    if F is None or G is None:
        return None
    return F, G, gap_coefficient(F, G, params["h"])


def accepts(params, claimed=Q(336, 625)):
    r = evaluate(params)
    if r is None:
        return False
    F, G, (gap0, cubic) = r
    return gap0 == 0 and cubic == claimed and F[0] == 3 and G[0] == 2


BASE = {"F": [(Q(2), Q(1)), (Q(8, 3), Q(1, 243))],
        "G": [(Q(3, 2), Q(1, 162)), (Q(7, 4), Q(1, 4))],
        "h": Q(1)}

failures = []


def need(ok, msg):
    print(("PASS  " if ok else "FAIL  ") + msg)
    if not ok:
        failures.append(msg)


F, G, (gap0, cubic) = evaluate(BASE)
need(F == [3, Q(2, 5), Q(42, 125), F[3]] and G[:3] == [2, Q(2, 5), Q(42, 125)], "bank jets F(4)=3, G(4)=2, F'=G'=2/5, F''=G''=42/125")
need(F[3] - G[3] == Q(63, 1250), "F''' - G''' = 63/1250")
need(gap0 == 0 and cubic == Q(336, 625), "route gap = (336/625)(t-1)^3 + O((t-1)^4)")
need(accepts(BASE), "published parameters accepted")

controls = {
    "connector free-flow time 1 -> 101/100": dict(BASE, h=Q(101, 100)),
    "F2 coefficient 1/243 -> 1/242": dict(BASE, F=[BASE["F"][0], (Q(8, 3), Q(1, 242))]),
    "G3 free-flow time 7/4 -> 9/5": dict(BASE, G=[BASE["G"][0], (Q(9, 5), Q(1, 4))]),
    "G1 coefficient 1/162 -> 1/160": dict(BASE, G=[(Q(3, 2), Q(1, 160)), BASE["G"][1]]),
    "F0 coefficient 1 -> 1/2": dict(BASE, F=[(Q(2), Q(1, 2)), BASE["F"][1]]),
}
for name, p in controls.items():
    need(not accepts(p), f"corruption rejected: {name}")
need(not accepts(BASE, claimed=Q(337, 625)), "corruption rejected: claimed gap coefficient 337/625")
if failures:
    print(f"\n{len(failures)} negative-control failure(s).")
    sys.exit(1)
print("\nAll negative controls behaved as expected.")
