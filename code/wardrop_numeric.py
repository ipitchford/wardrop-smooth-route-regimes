"""mpmath equilibrium / optimum solver for two-bank networks with a connector.

Commodity A (demand dA) from origin to destination: direct bank, or connector
then shared bank.  Commodity B (demand dB): shared bank only.  Edge costs
a + b z^4 (physical).  Bank splits equalise delays (Wardrop) or marginal
costs a + 5 b z^4 (system optimum).  All root finding is safeguarded Newton
on scalar monotone equations with bisection fallback.
"""
import mpmath as mp

DPS = 60
mp.mp.dps = DPS


def fmt(v, sig=25):
    return mp.nstr(mp.mpf(v), sig, min_fixed=-5, max_fixed=8)


def newton_decreasing(f, df, lo, hi, x0=None, tol_exp=None):
    """Root of a strictly decreasing f on [lo,hi] with f(lo)>=0>=f(hi)."""
    tol = mp.mpf(10) ** (-(tol_exp or (mp.mp.dps - 6)))
    x = x0 if (x0 is not None and lo < x0 < hi) else (lo + hi) / 2
    for _ in range(400):
        fx = f(x)
        if fx > 0:
            lo = x
        else:
            hi = x
        d = df(x)
        xn = x - fx / d if d != 0 else (lo + hi) / 2
        if not (lo < xn < hi):
            xn = (lo + hi) / 2
        if abs(xn - x) <= tol * max(abs(x), mp.mpf(10) ** (-mp.mp.dps + 15)):
            fx = f(xn); d = df(xn)
            return xn - fx / d if d != 0 else xn
        if hi - lo <= tol * max(abs(lo), abs(hi), mp.mpf(10) ** (-mp.mp.dps + 15)):
            return (lo + hi) / 2
        x = xn
    raise ArithmeticError("newton_decreasing did not converge")


class Bank:
    def __init__(self, edges, k=1):
        """edges: list of (a,b); k=1 equilibrium delays a+b z^4, k=5 marginal a+5b z^4."""
        self.e = [(mp.mpf(a), mp.mpf(b) * k) for a, b in edges]

    def loads_at(self, lam):
        return [max(mp.mpf(0), (lam - a) / b) ** mp.mpf(0.25) if lam > a else mp.mpf(0) for a, b in self.e]

    def total(self, lam):
        return sum(self.loads_at(lam))

    def dtotal(self, lam):
        return sum(((lam - a) / b) ** mp.mpf(-0.75) / (4 * b) for a, b in self.e if lam > a)

    def delay(self, T):
        """bank delay lam(T) with total(lam)=T, T>0."""
        T = mp.mpf(T)
        lo = min(a for a, b in self.e)
        hi = lo + 1
        while self.total(hi) < T:
            hi = lo + 2 * (hi - lo)
        return newton_decreasing(lambda l: T - self.total(l), lambda l: -self.dtotal(l), lo, hi)

    def ddelay(self, T):
        return 1 / self.dtotal(self.delay(T))


class Network:
    """direct edges, shared edges, connector (a_h, b_h) cost a_h + b_h x^4 (+eps)."""

    def __init__(self, direct, shared, conn=(1, 1), k=1):
        self.k = k
        self.F = Bank(direct, k)
        self.G = Bank(shared, k)
        self.h = (mp.mpf(conn[0]), mp.mpf(conn[1]) * k)
        self.direct, self.shared, self.conn = direct, shared, conn

    def hcost(self, x, eps=0):
        return self.h[0] + eps + self.h[1] * x ** 4

    def solve(self, dA, dB=None, eps=0):
        """Returns connector flow x (for k=5: system-optimal connector flow)."""
        dA = mp.mpf(dA); dB = dA if dB is None else mp.mpf(dB)
        phi = lambda x: self.F.delay(dA - x) - self.G.delay(dB + x) - self.hcost(x, eps)
        if phi(mp.mpf(0)) <= 0:
            return mp.mpf(0)
        dphi = lambda x: -self.F.ddelay(dA - x) - self.G.ddelay(dB + x) - 4 * self.h[1] * x ** 3
        x0 = phi(mp.mpf(0)) / (-dphi(mp.mpf(0)))
        return newton_decreasing(phi, dphi, mp.mpf(0), dA * (1 - mp.mpf(10) ** -30), x0)

    def flows(self, dA, dB=None, x=None, eps=0):
        dA = mp.mpf(dA); dB = dA if dB is None else mp.mpf(dB)
        if x is None:
            x = self.solve(dA, dB, eps)
        zF = self.F.loads_at(self.F.delay(dA - x))
        zG = self.G.loads_at(self.G.delay(dB + x))
        return x, zF, zG

    def physical_cost(self, x, zF, zG, eps=0):
        c = lambda a, b, z: z * (mp.mpf(a) + mp.mpf(b) * z ** 4)
        tot = sum(c(a, b, z) for (a, b), z in zip(self.direct, zF))
        tot += sum(c(a, b, z) for (a, b), z in zip(self.shared, zG))
        tot += x * (mp.mpf(self.conn[0]) + eps + mp.mpf(self.conn[1]) * x ** 4)
        return tot

    def potential(self, x, zF, zG):
        P = lambda a, b, z: mp.mpf(a) * z + mp.mpf(b) * z ** 5 / 5
        tot = sum(P(a, b, z) for (a, b), z in zip(self.direct, zF))
        tot += sum(P(a, b, z) for (a, b), z in zip(self.shared, zG))
        return tot + P(self.conn[0], self.conn[1], x)

    def eq_cost_formula(self, dA, dB=None, x=None, eps=0):
        """E = dA*cost_A + dB*cost_B with cost_A=F(dA-x), cost_B=G(dB+x)."""
        dA = mp.mpf(dA); dB = dA if dB is None else mp.mpf(dB)
        if x is None:
            x = self.solve(dA, dB, eps)
        return dA * self.F.delay(dA - x) + dB * self.G.delay(dB + x)


FIVE_DIRECT = [(2, 1), (mp.mpf(8) / 3, mp.mpf(1) / 243)]
FIVE_SHARED = [(mp.mpf(3) / 2, mp.mpf(1) / 162), (mp.mpf(7) / 4, mp.mpf(1) / 4)]
