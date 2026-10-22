#!/usr/bin/env python3
"""Exact finite fingerprints for the classical weighted-basis identity.

This is a reproducibility check, not a proof of the all-matrix lemma.
All checks use Fraction arithmetic and remain enabled under python -O -S.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rank(matrix):
    a = [list(map(F, row)) for row in matrix]
    if not a:
        return 0
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][j]
        a[r] = [x / scale for x in a[r]]
        for i in range(r + 1, len(a)):
            scale = a[i][j]
            a[i] = [x - scale * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def solve(matrix, rhs):
    n = len(matrix)
    a = [list(map(F, row)) + [F(rhs[i])] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        require(pivot is not None, "Singular square system")
        a[j], a[pivot] = a[pivot], a[j]
        scale = a[j][j]
        a[j] = [x / scale for x in a[j]]
        for i in range(n):
            if i != j:
                scale = a[i][j]
                a[i] = [x - scale * y for x, y in zip(a[i], a[j])]
    return [row[-1] for row in a]


def determinant(matrix):
    n = len(matrix)
    a = [list(map(F, row)) for row in matrix]
    result = F(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        result *= a[j][j]
        for i in range(j + 1, n):
            scale = a[i][j] / a[j][j]
            a[i] = [x - scale * y for x, y in zip(a[i], a[j])]
    return result


def matvec(matrix, vector):
    return [sum((F(x) * y for x, y in zip(row, vector)), F(0)) for row in matrix]


def check_instance(matrix, diagonal, seed):
    """Compare direct KKT elimination with independently enumerated bases."""
    n = len(diagonal)
    rhs = matvec(matrix, seed)
    selected = []
    for i in range(len(matrix)):
        trial = selected + [i]
        if rank([matrix[j] for j in trial]) > len(selected):
            selected.append(i)
    reduced = [matrix[i] for i in selected]
    reduced_rhs = [rhs[i] for i in selected]
    k = len(selected)
    if k == 0:
        require(not any(rhs), "A zero-rank matrix has a nonzero residual")
        return 0
    inverse_diagonal = [1 / F(d) for d in diagonal]
    gram = [[sum((F(reduced[i][j]) * inverse_diagonal[j] * reduced[h][j]
                  for j in range(n)), F(0)) for h in range(k)] for i in range(k)]
    multiplier = solve(gram, reduced_rhs)
    direct = [inverse_diagonal[j] * sum((F(reduced[i][j]) * multiplier[i]
              for i in range(k)), F(0)) for j in range(n)]
    require(matvec(matrix, direct) == rhs, "Direct KKT solution is infeasible")

    total_weight = F(0)
    numerator = [F(0)] * n
    largest_inverse_frobenius_squared = F(0)
    basis_count = 0
    for indices in combinations(range(n), k):
        square = [[row[j] for j in indices] for row in reduced]
        det = determinant(square)
        if not det:
            continue
        basic = solve(square, reduced_rhs)
        weight = det * det
        for j in indices:
            weight *= inverse_diagonal[j]
        total_weight += weight
        basis_count += 1
        for j, value in zip(indices, basic):
            numerator[j] += weight * value
        frobenius_squared = F(0)
        for i in range(k):
            inverse_column = solve(square, [F(h == i) for h in range(k)])
            frobenius_squared += sum((x * x for x in inverse_column), F(0))
        largest_inverse_frobenius_squared = max(
            largest_inverse_frobenius_squared, frobenius_squared)
    average = [x / total_weight for x in numerator]
    require(direct == average, "Convex basis identity differs from direct KKT solve")
    require(total_weight == determinant(gram), "Cauchy-Binet denominator mismatch")
    lhs = sum((x * x for x in direct), F(0))
    rhs_bound = largest_inverse_frobenius_squared * sum((x * x for x in rhs), F(0))
    require(lhs <= rhs_bound, "Uniform Frobenius upper bound fails")
    return basis_count


def main():
    matrices = [
        [[1, 1]],
        [[1, 0, 1], [0, 1, 1]],
        [[1, 0, 1, 2], [0, 1, 1, 1]],
        [[1, 1, 0], [2, 2, 0]],
    ]
    weight_patterns = [
        [F(1), F(1), F(1), F(1)],
        [F(1, 10**12), F(10**12), F(7, 11), F(23, 5)],
        [F(10**12), F(1, 10**12), F(10**6), F(1, 10**6)],
    ]
    seed = [F(-3), F(4), F(5), F(-6)]
    instances = 0
    basis_comparisons = 0
    cases = []
    for matrix in matrices:
        n = len(matrix[0])
        for size in range(n + 1):
            for indices in combinations(range(n), size):
                restricted = [[row[j] for j in indices] for row in matrix]
                for pattern in weight_patterns:
                    basis_comparisons += check_instance(
                        restricted, [pattern[j] for j in indices], [seed[j] for j in indices])
                    instances += 1
        cases.append({"matrix": matrix, "column_subsets": 2**n})
    report = {
        "status": "PASS",
        "arithmetic": "exact fractions; checks remain active under -O -S",
        "scope": "Finite fingerprints only; the general lemma is proved in fixed_regime.tex",
        "checks": ["KKT feasibility", "convex basis identity", "Cauchy-Binet denominator",
                   "weight-independent norm bound", "dependent rows", "zero columns",
                   "all column restrictions, including empty restrictions"],
        "instances": instances,
        "basis_comparisons": basis_comparisons,
        "weight_patterns": [[str(d) for d in pattern] for pattern in weight_patterns],
        "cases": cases,
    }
    output = Path(__file__).with_name("weighted_basis_verification.json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
