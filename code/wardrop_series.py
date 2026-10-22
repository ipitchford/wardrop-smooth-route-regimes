"""Exact truncated power-series toolkit (fractions.Fraction only).

Shared by verify_practice.py, verify_optimum_transition.py and
verify_small_examples.py.  A series is a list of Fractions [c0, c1, ..., cn].
Bank model: an edge with cost a + b z^p, at a bank delay level+s, carries
z = w (1+u s)^(1/p) when a = level - 1/u and b = 1/(u w^p).
"""
from fractions import Fraction as Q
from math import comb


class CheckFailed(ArithmeticError):
    pass


def check(cond, msg):
    if not cond:
        raise CheckFailed(msg)


def frac(v):
    v = Q(v)
    return str(v.numerator) if v.denominator == 1 else f"{v.numerator}/{v.denominator}"


def zeros(n):
    return [Q(0)] * (n + 1)


def add(a, b, n):
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n + 1)]


def scale(a, c, n):
    return [Q(c) * (a[i] if i < len(a) else 0) for i in range(n + 1)]


def mul(a, b, n):
    out = zeros(n)
    for i, ai in enumerate(a[:n + 1]):
        if ai:
            for j, bj in enumerate(b[:n + 1 - i]):
                if bj:
                    out[i + j] += ai * bj
    return out


def power(a, k, n):
    out = [Q(1)] + [Q(0)] * n
    for _ in range(k):
        out = mul(out, a, n)
    return out


def compose(a, b, n):
    """a(b(d)) with b[0] == 0."""
    check(b[0] == 0, "composition needs zero inner constant")
    out = zeros(n)
    bp = [Q(1)] + [Q(0)] * n
    for ak in a[:n + 1]:
        if ak:
            for i in range(n + 1):
                out[i] += ak * bp[i]
        bp = mul(bp, b, n)
    return out


def deriv(a):
    return [i * a[i] for i in range(1, len(a))] + [Q(0)]


def binom(r, k):
    out = Q(1)
    for j in range(k):
        out *= (Q(r) - j) / (j + 1)
    return out


def one_plus_us_pow(u, r, n):
    """(1+u s)^r as a series in s."""
    return [binom(r, k) * Q(u) ** k for k in range(n + 1)]


def bank_load_series(edges, n, p=4):
    """q(s) = sum w (1+u s)^(1/p); edges = [(w,u),...]."""
    out = zeros(n)
    for w, u in edges:
        out = add(out, scale(one_plus_us_pow(u, Q(1, p), n), w, n), n)
    return out


def bank_potential_series(edges, level, n, p=4):
    """sum_e (a_e z_e + b_e z_e^(p+1)/(p+1)) as a series in s."""
    out = zeros(n)
    for w, u in edges:
        a = Q(level) - Q(1, u)
        out = add(out, scale(one_plus_us_pow(u, Q(1, p), n), a * w, n), n)
        out = add(out, scale(one_plus_us_pow(u, Q(p + 1, p), n), Q(w, u * (p + 1)), n), n)
    return out


def edge_load_series(w, u, n, p=4):
    return scale(one_plus_us_pow(u, Q(1, p), n), w, n)


def inverse_series(q, n):
    """Given q(s)=q0+q1 s+..., return s(y) with q(s(y)) = q0 + y, s(0)=0."""
    check(q[1] > 0, "inverse needs positive slope")
    sh = list(q[:n + 1]); sh[0] = Q(0)
    s = zeros(n)
    for k in range(1, n + 1):
        r = compose(sh, s, k)[k]
        s[k] = ((1 if k == 1 else 0) - r) / q[1]
    check(compose(sh, s, n) == [Q(0), Q(1)] + [Q(0)] * (n - 1), "inverse failed")
    return s


def physical_check(edges, level, p=4):
    """Return physical (a,b,ref_load) and check a+b z^p = level at ref load."""
    out = []
    for w, u in edges:
        a = Q(level) - Q(1, u); b = Q(1, u * w ** p)
        check(a + b * Q(w) ** p == level, "reference latency mismatch")
        out.append((a, b, Q(w)))
    return out


class TwoBank:
    """Direct bank (level Ld) for commodity A: origin->dest, shared bank (level Ls)
    used by commodity B and by A via the connector (cost 1 + x^p).
    Demand D(1+delta) for each commodity.  Requires Ld - Ls - 1 = 0."""

    def __init__(self, dedges, Ld, sedges, Ls, D, n, p=4):
        self.n, self.D, self.p = n, Q(D), p
        self.dedges, self.sedges, self.Ld, self.Ls = dedges, sedges, Q(Ld), Q(Ls)
        check(self.Ld - self.Ls - 1 == 0, "levels not matched")
        qd, qs = bank_load_series(dedges, n, p), bank_load_series(sedges, n, p)
        check(qd[0] == D and qs[0] == D, "reference demand mismatch")
        self.qd, self.qs = qd, qs
        self.sd, self.ss = inverse_series(qd, n), inverse_series(qs, n)   # s as function of load deviation y
        self.beta_d, self.beta_s = qd[1], qs[1]

    def bank_s(self, x):
        n, D = self.n, self.D
        yd = add([Q(0), D], scale(x, -1, n), n)
        ys = add([Q(0), D], x, n)
        return compose(self.sd, yd, n), compose(self.ss, ys, n)

    def connector(self):
        """x(delta) solving s_d(D delta - x) - s_s(D delta + x) - x^p = 0 (formal branch)."""
        n, p = self.n, self.p
        x = zeros(n)
        der = -(self.sd[1] + self.ss[1])

        def res(x, k):
            a, b = TwoBank.bank_s(self._k(k), x[:k + 1])
            return [a[i] - b[i] - power(x[:k + 1], p, k)[i] for i in range(k + 1)]
        for k in range(1, n + 1):
            x[k] = -res(x, k)[k] / der
        check(all(c == 0 for c in res(x, n)), "connector residual nonzero")
        return x

    def _k(self, k):
        o = object.__new__(TwoBank)
        o.n, o.D, o.sd, o.ss = k, self.D, self.sd[:k + 1], self.ss[:k + 1]
        return o

    def gap(self):
        """D(delta) = F(D(1+delta)) - G(D(1+delta)) - 1 with x = 0."""
        a, b = self.bank_s(zeros(self.n))
        return add(a, scale(b, -1, self.n), self.n)

    def cost(self, x):
        n = self.n
        a, b = self.bank_s(x)
        tot = add(a, b, n); tot[0] += self.Ld + self.Ls
        return mul([self.D, self.D], tot, n)

    def potential(self, x):
        n, p = self.n, self.p
        a, b = self.bank_s(x)
        Vd = compose_const(bank_potential_series(self.dedges, self.Ld, n, p), a, n)
        Vs = compose_const(bank_potential_series(self.sedges, self.Ls, n, p), b, n)
        conn = add(x, scale(power(x, p + 1, n), Q(1, p + 1), n), n)
        return add(add(Vd, Vs, n), conn, n)

    def commodity_costs(self, x):
        a, b = self.bank_s(x)
        A = list(a); A[0] += self.Ld
        B = list(b); B[0] += self.Ls
        return A, B

    def edge_loads(self, x):
        a, b = self.bank_s(x)
        out = []
        for name, edges, s in (("direct", self.dedges, a), ("shared", self.sedges, b)):
            for w, u in edges:
                out.append((f"{name}_w{w}_u{u}", compose_const(edge_load_series(w, u, self.n, self.p), s, self.n)))
        return out


def compose_const(a, b, n):
    """a(b) where b has zero constant term (alias)."""
    return compose(a, b, n)


def binomial_family_banks(m):
    """Five-edge-type family: weights C(m,j), u=j+1; even j -> level-3 bank, odd j -> level-2 bank."""
    F = [(comb(m, j), j + 1) for j in range(0, m + 1, 2)]
    G = [(comb(m, j), j + 1) for j in range(1, m + 1, 2)]
    return F, G


def first_nonzero(a, start=0):
    for i in range(start, len(a)):
        if a[i] != 0:
            return i, a[i]
    return None, Q(0)
