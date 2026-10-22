#!/usr/bin/env python3
"""Adversarial exact aggregate-series reconstruction of the multicontact examples.

Uses the supplied physical coefficients and rational physical boundary loads,
not the claimed route-gap or aggregate coefficients, to construct inverse
throughput series and solve the scalar Wardrop equation. The final comparison
is against the separately generated multicontact certificate.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location("main_exact_arithmetic", ROOT/"code"/"verify_exact.py")
exact = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exact)
check, mul, compose, power = exact.check, exact.mul, exact.compose, exact.power


def reconstruct(case, cj):
    p, m = case["p"], case["m"]
    n = 2*m-1
    contact = case["contacts"][cj]
    t0 = Q(contact["t"])
    qseries = {}
    for bank in ["F", "G"]:
        level = Q(case["bank_reference_"+bank+"_delay"]) + cj
        qs = [Q(0)]*(n+1)
        for e in case["physical_edges"]:
            if e["bank"] != bank:
                continue
            ae, be = Q(e["free_flow"]), Q(e["congestion_coefficient"])
            z0 = Q(e["reference_load"] if cj==0 else e["second_boundary_load"])
            check(ae+be*z0**p==level, "Boundary load violates physical cost")
            for k in range(n+1):
                qs[k] += z0*exact.binomial(Q(1,p),k)/(level-ae)**k
        check(qs[0] == t0, "Boundary total demand mismatch")
        qseries[bank] = qs
    # delta here is relative demand displacement: t=t0*(1+delta).
    f = exact.inverse_series(qseries["F"], t0, n)
    g = exact.inverse_series(qseries["G"], t0, n)
    x = [Q(0)]*(n+1)
    jac = -(f[1]+g[1])/t0

    def bank_costs(xx,k):
        minus, plus = [-z/t0 for z in xx[:k+1]], [z/t0 for z in xx[:k+1]]
        minus[1]+=1
        plus[1]+=1
        return compose(f,minus,k), compose(g,plus,k)

    def residual(xx,k):
        ff,gg = bank_costs(xx,k)
        hp = power(xx,p,k)
        return [ff[j]-gg[j]-hp[j] for j in range(k+1)]

    for k in range(1,n+1):
        x[k] = -residual(x,k)[k]/jac
    check(all(z==0 for z in residual(x,n)), "New implicit residual check")

    def cost(xx):
        ff,gg = bank_costs(xx,n)
        return mul([t0,t0],[ff[j]+gg[j] for j in range(n+1)],n)

    active, inactive = cost(x), cost([Q(0)]*(n+1))
    diff = [active[j]-inactive[j] for j in range(n+1)]
    a_actual = (f[m]-g[m])/t0**m
    x_actual = x[m]/t0**m
    coefficient = diff[n]/t0**n
    check(all(z==0 for z in x[:m]), "Earlier connector term")
    check(all(z==0 for z in diff[:n]), "Earlier aggregate term")
    check(a_actual == Q(contact["route_gap_leading_coefficient"]), "Gap coefficient disagrees")
    check(x_actual == Q(contact["signed_connector_leading_coefficient"]), "Connector coefficient disagrees")
    check(coefficient == -Q(contact["aggregate_coefficient_magnitude"]), "Aggregate coefficient disagrees")
    orientation = 1 if a_actual>0 else -1
    jump = orientation*factorial(n)*coefficient
    check(jump == Q(contact["right_minus_left_aggregate_derivative_jump"]), "Derivative jump disagrees")
    return {"p":p,"m":m,"contact":cj+1,"t0":str(t0),
            "flow_order":m,"aggregate_order":n,
            "signed_gap_coefficient":str(a_actual),
            "signed_flow_coefficient":str(x_actual),
            "aggregate_active_coefficient":str(coefficient),
            "right_minus_left_derivative_jump":str(jump),
            "status":"PASS"}


def main():
    source = ROOT/"checks"/"moments"/"two_contacts_exact.json"
    cases = json.loads(source.read_text())["cases"]
    results = [reconstruct(case,j) for case in cases for j in range(2)]
    (HERE/"multicontact_reconstruction.json").write_text(json.dumps({"method":"independent exact physical-series reconstruction","results":results},indent=2)+"\n")
    print(json.dumps({"status":"PASS","reconstructed_contacts":len(results),"examples":[[x["p"],x["m"]] for x in cases]},indent=2))


if __name__=="__main__":
    main()
