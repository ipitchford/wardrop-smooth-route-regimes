#!/usr/bin/env python3
"""Exact finite certificates for quartic routing examples.

This verifier uses only Python's standard library. It independently derives
the jets from the physical link coefficients by formal series inversion and
implicit-series substitution. The arbitrary-order theorem additionally
requires the finite-difference and implicit-function arguments in the paper.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from math import comb, factorial, isqrt
from pathlib import Path
import json


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def frac(value):
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def mul(a, b, n):
    out = [Q(0)] * (n + 1)
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
    check(b[0] == 0, 'Series composition requires zero inner constant')
    out = [Q(0)] * (n + 1)
    bp = [Q(1)] + [Q(0)] * n
    for ak in a[:n + 1]:
        if ak:
            for i in range(n + 1):
                out[i] += ak * bp[i]
        bp = mul(bp, b, n)
    return out


def binomial(r, k):
    out = Q(1)
    for j in range(k):
        out *= (r - j) / (j + 1)
    return out


def root4(q):
    q = Q(q)
    check(q >= 0, 'Negative rational fourth-root input')
    a = isqrt(isqrt(q.numerator))
    b = isqrt(isqrt(q.denominator))
    check(a ** 4 == q.numerator and b ** 4 == q.denominator,
          'Rational fourth root is not exact')
    return Q(a, b)


def physical_links(m):
    check(m >= 1 and m % 2 == 1, 'm must be positive and odd')
    edges = []
    for j in range(m + 1):
        w, a = comb(m, j), j + 1
        group = 'F' if j % 2 == 0 else 'G'
        level = 3 if group == 'F' else 2
        edges.append({'j': j, 'group': group, 'a': Q(level) - Q(1, a),
                      'b': Q(1, a * w ** 4), 'baseline_load': Q(w)})
    return edges


def throughput_series(edges, group, n):
    """Derive sum ((lambda-a)/b)^(1/4) from physical costs."""
    level = Q(3 if group == 'F' else 2)
    out = [Q(0)] * (n + 1)
    for edge in edges:
        if edge['group'] != group:
            continue
        a, b = edge['a'], edge['b']
        check(a > 0 and b > 0 and level > a, 'Admissibility failure')
        w = root4((level - a) / b)
        check(a + b * w ** 4 == level, 'Physical latency mismatch')
        for k in range(n + 1):
            out[k] += w * binomial(Q(1, 4), k) / (level - a) ** k
    return out


def inverse_series(q, demand, n):
    """Solve q(s(delta)) = demand*(1+delta) with s(0)=0."""
    check(q[0] == demand and q[1] > 0, 'Inverse base conditions failed')
    shifted = list(q); shifted[0] = Q(0)
    result = [Q(0)] * (n + 1)
    for k in range(1, n + 1):
        residual = compose(shifted, result, k)[k]
        desired = demand if k == 1 else Q(0)
        result[k] = (desired - residual) / q[1]
    composed = compose(shifted, result, n)
    check(composed == [Q(0), Q(demand)] + [Q(0)] * (n - 1),
          'Inverse does not reproduce its physical throughput')
    return result


def connector_series(f, g, demand, n):
    """Solve F(q-x)-G(q+x)-1-x^4=0 by implicit substitution."""
    x = [Q(0)] * (n + 1)
    derivative = -(f[1] + g[1]) / demand
    check(derivative < 0, 'Connector Jacobian has wrong sign')

    def residual(x, k):
        minus, plus = [Q(0)] * (k + 1), [Q(0)] * (k + 1)
        for i in range(1, k + 1):
            minus[i] = -x[i] / demand
            plus[i] = x[i] / demand
        minus[1] += 1; plus[1] += 1
        fp, gp, hp = compose(f, minus, k), compose(g, plus, k), power(x, 4, k)
        return [fp[i] - gp[i] - hp[i] for i in range(k + 1)]

    for k in range(1, n + 1):
        x[k] = -residual(x, k)[k] / derivative
    check(all(c == 0 for c in residual(x, n)), 'Connector residual is nonzero')
    return x


def cost_series(f, g, x, demand, n):
    minus = [-c / demand for c in x]
    plus = [c / demand for c in x]
    minus[1] += 1; plus[1] += 1
    fp, gp = compose(f, minus, n), compose(g, plus, n)
    total = [a + b for a, b in zip(fp, gp)]
    total[0] = Q(5)
    return mul([Q(demand), Q(demand)], total, n)


def inspect_m(m):
    n = 2 * m - 1
    demand = Q(2 ** (m - 1))
    edges = physical_links(m)
    qf, qg = throughput_series(edges, 'F', n), throughput_series(edges, 'G', n)
    f, g = inverse_series(qf, demand, n), inverse_series(qg, demand, n)
    x = connector_series(f, g, demand, n)
    E = cost_series(f, g, x, demand, n)
    E0 = cost_series(f, g, [Q(0)] * (n + 1), demand, n)
    moments = [sum((-1) ** (m - j) * comb(m, j) * (j + 1) ** k
                   for j in range(m + 1)) for k in range(m + 1)]
    check(moments == [0] * m + [factorial(m)], 'Binomial moment mismatch')
    # These identities validate the finite-difference coefficient independently
    # from the throughput inversion above; the sign for all z uses the integral
    # of the mth derivative, whose factors are recorded here.
    factors = [1 - 4 * j for j in range(m)]
    P = 1
    for a in factors: P *= a
    check(P > 0, 'Odd-order fourth-root derivative sign failed')

    if m >= 3:
        d = Q(16 * P, (m + 2) ** (m + 1))
        alpha = Q(2 ** (m - 1) * P, (m + 2) ** m)
        ec = -Q(m * 2 ** (m + 3) * P ** 2, (m + 2) ** (2 * m + 1))
        check(all(f[k] == g[k] for k in range(m)), 'Earlier mismatch in macro latencies')
        check(f[m] - g[m] == d, 'Wrong gap leading coefficient')
        check(all(x[k] == 0 for k in range(m)) and x[m] == alpha,
              'Wrong connector onset order or coefficient')
        check(all(E[k] == E0[k] for k in range(2 * m - 1)),
              'Social cost differs before the claimed order')
        check(E[-1] - E0[-1] == ec, 'Wrong social-cost contact coefficient')
    else:
        d, alpha, ec = f[1] - g[1], x[1], E[1] - E0[1]
        check((d, alpha, ec) == (Q(2), Q(1, 3), -Q(2, 3)),
              'Transverse quartic control failed')

    out = {'m': m, 'physical_edge_count_including_connector': m + 2,
           'demand_each': frac(demand),
           'moment_sign_convention':'inverse-throughput G minus F; signs (-1)^(m-j)',
           'moments_through_m': moments,
           'fourth_root_derivative_factors': factors,
           'gap_coefficient': frac(d), 'connector_coefficient': frac(alpha),
           'equilibrium_cost_difference_coefficient': frac(ec),
           'E_at_1': frac(E[0]), 'E_left_derivative_at_1': frac(E0[1]),
           'E_right_derivative_at_1': frac(E[1]),
           'physical_edges': [{k: frac(v) if isinstance(v, Q) else v
                              for k, v in edge.items()} for edge in edges]}
    if m == 3:
        check([(e['a'], e['b']) for e in edges] ==
              [(Q(2),Q(1)), (Q(3,2),Q(1,162)),
               (Q(8,3),Q(1,243)), (Q(7,4),Q(1,4))],
              'Five-edge physical example differs from stated costs')
        jets = {}
        for name, qs in [('F',qf),('G',qg)]:
            first, second, third = qs[1], 2*qs[2], 6*qs[3]
            inv3 = 3*second**2 / first**5 - third/first**4
            jets[name] = {'inverse_first':frac(first), 'inverse_second':frac(second),
                          'inverse_third':frac(third), 'latency_first':frac(1/first),
                          'latency_second':frac(-second/first**3),
                          'latency_third':frac(inv3)}
        check(qf[1] == qg[1] == Q(5,2), 'First inverse jets wrong')
        check(2*qf[2] == 2*qg[2] == -Q(21,4), 'Second inverse jets wrong')
        check(6*qf[3] == Q(861,32) and 6*qg[3] == Q(231,8),
              'Third inverse jets wrong')
        check((d, alpha, ec, E[0], E[1]) ==
              (Q(336,625),Q(84,125),-Q(84672,78125),Q(20),Q(164,5)),
              'Named exact five-edge quantities wrong')
        out['jets'] = jets
        out['equilibrium_cost_fifth_derivative_jump'] = frac(factorial(5)*ec)
    return out


def perturbation_control():
    """Keep the base latency/load but perturb one quartic edge's slope."""
    edges = physical_links(3)
    # f1(z)=19/10+(11/10)z^4, keeping f1(1)=3.
    edges[0]['a'],edges[0]['b'] = Q(19,10),Q(11,10)
    qf,qg = throughput_series(edges,'F',1),throughput_series(edges,'G',1)
    f,g = inverse_series(qf,Q(4),1),inverse_series(qg,Q(4),1)
    x = connector_series(f,g,Q(4),1)
    E = cost_series(f,g,x,Q(4),1)
    E0 = cost_series(f,g,[Q(0),Q(0)],Q(4),1)
    d = f[1]-g[1]; jump = E[1]-E0[1]
    check(d == Q(8,545) and x[1] == Q(4,219), 'Perturbation lost transversality')
    check(jump == -Q(32,119355), 'Perturbed first derivative jump wrong')
    check(jump == -d*d/((f[1]+g[1])/4), 'Quadratic jump law wrong')
    return {'f1_a':frac(edges[0]['a']),'f1_b':frac(edges[0]['b']),
            'gap_first_derivative':frac(d),'connector_first_derivative':frac(x[1]),
            'equilibrium_cost_first_derivative_jump':frac(jump)}


def main():
    parser = argparse.ArgumentParser()
    default = Path(__file__).resolve().parents[1] / 'evidence' / 'exact_verification.json'
    parser.add_argument('--output',type=Path,default=default)
    args = parser.parse_args()
    records = [inspect_m(m) for m in (1,3,5,7)]
    report = {'status':'passed','arithmetic':'fractions.Fraction; exact integers',
              'coverage':'Finite certificates m=1,3,5,7; not a substitute for the all-odd-m proof',
              'cases':records,'transverse_perturbation':perturbation_control(),
              'checks':['Physical link positivity and exact baseline latencies',
                        'Formal inversion reproduces physical throughput series',
                        'Connector implicit residual vanishes to retained order',
                        'First nonzero cost difference occurs at the stated order',
                        'Independent integer binomial moments and derivative-factor signs',
                        'Transverse m=1 and perturbed m=3 first-derivative controls']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'cases':[c['m'] for c in records],
                      'output':str(args.output)},indent=2))


if __name__=='__main__': main()
