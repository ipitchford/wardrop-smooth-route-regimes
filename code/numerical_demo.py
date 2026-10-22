#!/usr/bin/env python3
"""Numerical sanity checks, explicitly separate from exact proof acceptance.

Requires NumPy and SciPy. The Wardrop side uses scalar inverse-throughput
bisection; the social optimum independently solves the physical marginal
cost equations in three free flow variables. Their comparison checks the
common-quartic scaling without building that identity into the optimiser.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import root


# Order: direct f1, direct f2, shared g1, shared g2, connector h.
A = np.array([2.,8/3,3/2,7/4,1.])
B = np.array([1.,1/243,1/162,1/4,1.])


def checked(condition, message):
    if not condition:
        raise ArithmeticError(message)


def bisect(fun, lo, hi):
    fl,fh=fun(lo),fun(hi)
    checked(fl*fh<=0, 'Bisection not bracketed')
    if fl==0: return lo
    if fh==0: return hi
    for _ in range(150):
        mid=(lo+hi)/2
        if mid==lo or mid==hi: break
        fm=fun(mid)
        if fl*fm<=0:
            hi=mid
        else:
            lo=mid; fl=fm
    return (lo+hi)/2


def throughput(level, group, a=A, b=B):
    idx=(0,1) if group=='F' else (2,3)
    return sum(max((level-a[i])/b[i],0.)**.25 for i in idx)


def macro(total,group,a=A,b=B):
    idx=(0,1) if group=='F' else (2,3)
    low=min(a[i] for i in idx)
    high=max(a[i]+b[i]*total**4 for i in idx)+1.
    return bisect(lambda level:throughput(level,group,a,b)-total,low,high)


def macro_slope(total,group,a=A,b=B):
    level=macro(total,group,a,b)
    idx=(0,1) if group=='F' else (2,3)
    active=[i for i in idx if level>a[i]]
    inverse_slope=sum(1/(4*b[i]*((level-a[i])/b[i])**.75) for i in active)
    return 1/inverse_slope


def wardrop(t,a=A,b=B):
    total=4*t
    def gap(x):
        return macro(total-x,'F',a,b)-a[4]-b[4]*x**4-macro(total+x,'G',a,b)
    if t==1 and np.array_equal(a,A) and np.array_equal(b,B):
        return np.array([1.,3.,3.,1.,0.]),3.,2.
    g0=gap(0.)
    x=0. if g0<=0 else bisect(gap,0.,total)
    lf,lg=macro(total-x,'F',a,b),macro(total+x,'G',a,b)
    z=np.array([max((lf-a[0])/b[0],0.)**.25,
                max((lf-a[1])/b[1],0.)**.25,
                max((lg-a[2])/b[2],0.)**.25,
                max((lg-a[3])/b[3],0.)**.25,x])
    checked(np.max(np.abs([z[0]+z[1]+z[4]-total,
                           z[2]+z[3]-z[4]-total]))<2e-10,
            'Wardrop conservation failed')
    return z,lf,lg


def direct_social_optimum(t):
    total=4*t
    def loads(v):
        u,vv,x=v
        return np.array([u,total-x-u,vv,total+x-vv,x])
    def kkt(v):
        z=loads(v); marginal=A+5*B*z**4
        return np.array([marginal[0]-marginal[1],
                         marginal[2]-marginal[3],
                         marginal[4]+marginal[3]-marginal[1]])
    initial=np.array([.85*t,2.9*t,.007*t])
    solution=root(kkt,initial,tol=1e-11)
    z=loads(solution.x)
    residual=float(np.max(np.abs(kkt(solution.x))))
    checked(residual<1e-9 and np.min(z)>0,
            'Direct social KKT solution is not an interior optimum')
    # Each edge objective a*z+b*z^5 has positive second derivative for z>0;
    # convexity plus feasible interior stationarity makes this the global
    # optimum. The residual check remains numerical, not a certificate.
    value=float(np.sum(z*(A+B*z**4)))
    marginal=A+5*B*z**4
    derivative=4*float(marginal[0]+marginal[2])
    return z,value,derivative,residual


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,
                        default=Path(__file__).resolve().parents[1]/'evidence')
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    kappa=5**.25
    zscaled,lf,lg=wardrop(kappa)
    zopt,C,Cprime,residual=direct_social_optimum(1.)
    discrepancy=float(np.max(np.abs(zopt-zscaled/kappa)))
    checked(discrepancy<1e-8,'Common-power scaling numerical check failed')
    E,Eprime=20.,164/5
    checked(0<C<E,'Social value must strictly improve the Wardrop value')
    report={'status':'passed','evidence_type':'floating-point sanity checks only',
            'exact_proof_dependency':False,'kappa':kappa,
            'wardrop_at_kappa_loads':zscaled.tolist(),
            'direct_social_at_1_loads':zopt.tolist(),
            'max_load_difference_direct_vs_scaling':discrepancy,
            'direct_social_KKT_max_residual':residual,
            'C_at_1':C,'Cprime_at_1':Cprime,
            'E_at_1':E,'Eprime_at_1':Eprime,'PoA_at_1':E/C,
            'PoA_common_first_derivative_at_1':(Eprime*C-E*Cprime)/C**2,
            'PoA_fifth_derivative_jump':-120*84672/(78125*C)}
    rows=[]
    for t in [.75,.9,.95,.99,1.,1.01,1.05,1.1,1.2,1.4,1.5]:
        ze,lf,lg=wardrop(t)
        zo,c,cp,r=direct_social_optimum(t)
        e=float(np.sum(ze*(A+B*ze**4)))
        e_routes=4*t*(lf+lg)
        checked(abs(e-e_routes)<1e-9,'Edge and path social accounting disagree')
        rows.append({'t':t,'equilibrium_connector':float(ze[4]),
                     'optimal_connector':float(zo[4]),'equilibrium_social_cost':e,
                     'optimal_social_cost':c,'price_of_anarchy':e/c})
    # A constant-plus-quartic perturbation preserving the transition point
    # changes one physical slope; exact first-derivative values are verified
    # separately in verify_exact.py.
    ap,bp=A.copy(),B.copy();ap[0]=1.9;bp[0]=1.1
    pert=[]
    for t in [.9999,1.0001]:
        zz,ll,gg=wardrop(t,ap,bp)
        pert.append({'t':t,'connector':float(zz[4]),'cost':4*t*(ll+gg)})
    checked(pert[0]['connector']==0 and pert[1]['connector']>0,
            'Transverse perturbed support control failed')
    report['fixed_breakpoint_direct_edge_perturbation']={
        'experiment':'Preserve f1(1)=3 while changing its physical slope',
        'changed_edge':'f1',
        'free_flow_constant':float(ap[0]),
        'quartic_coefficient':float(bp[0]),
        'fixed_breakpoint':1.,
        'samples':pert}
    # Main-paper sensitivity experiment: change only h's free-flow constant
    # to 1+epsilon. The breakpoint solves D(t)=epsilon. Its exact jump formula
    # is evaluated independently from the asymptotic power law.
    sensitivity=[]
    onset_coefficient=336/625
    jump_constant=(5/4)*9*onset_coefficient**(2/3)
    for eps in [-1e-4,-1e-6,-1e-8,1e-8,1e-6,1e-4]:
        fun=lambda t:macro(4*t,'F')-macro(4*t,'G')-1-eps
        te=bisect(fun,.8,1.3)
        fs,gs=macro_slope(4*te,'F'),macro_slope(4*te,'G')
        dt=4*(fs-gs)
        jump=te*dt*dt/(fs+gs)
        shift_pred=np.copysign((abs(eps)/onset_coefficient)**(1/3),eps)
        jump_pred=jump_constant*abs(eps)**(4/3)
        sensitivity.append({'connector_constant_change':eps,'breakpoint':te,
                            'breakpoint_shift':te-1,'leading_shift':float(shift_pred),
                            'shift_ratio':(te-1)/shift_pred,
                            'positive_cost_derivative_drop':jump,
                            'leading_derivative_drop':jump_pred,
                            'drop_ratio':jump/jump_pred})
    report['connector_free_flow_sensitivity']=sensitivity
    checked(all(row['positive_cost_derivative_drop']>0 for row in sensitivity),
            'Connector perturbation failed to restore a strict derivative drop')
    (args.output_dir/'numerical_sanity.json').write_text(json.dumps(report,indent=2)+'\n')
    with (args.output_dir/'numerical_curve.csv').open('w',newline='') as fp:
        writer=csv.DictWriter(fp,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
