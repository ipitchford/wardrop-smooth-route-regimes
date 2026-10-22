#!/usr/bin/env python3
"""Exact two-contact common-power network checks; standard library only.

Run from any directory. All mathematical checks use Fraction and explicit
exceptions, including under python -O -S. The JSON includes exact physical
edge coefficients, all derivative contacts, and local asymptotic constants.
These finite checks corroborate the accompanying general analytic proof.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def falling(a, k):
    ans = Q(1)
    for j in range(k):
        ans *= a-j
    return ans

def times_linear(v, root):
    out = [Q(0)]*(len(v)+1)
    for j, x in enumerate(v):
        out[j] -= root*x
        out[j+1] += x
    return out

def peval(v, z):
    ans = Q(0)
    for x in reversed(v):
        ans = ans*z+x
    return ans

def rank(matrix):
    a = [list(map(Q, row)) for row in matrix]
    nr, nc = len(a), len(a[0])
    rr = 0
    for cc in range(nc):
        pivot = next((j for j in range(rr,nr) if a[j][cc]),None)
        if pivot is None:
            continue
        a[rr],a[pivot] = a[pivot],a[rr]
        d = a[rr][cc]
        a[rr] = [x/d for x in a[rr]]
        for j in range(nr):
            if j!=rr and a[j][cc]:
                fac = a[j][cc]
                a[j] = [x-fac*y for x,y in zip(a[j],a[rr])]
        rr += 1
        if rr == nr:
            break
    return rr

def encode(obj):
    if isinstance(obj,Q):
        return str(obj)
    if isinstance(obj,dict):
        return {k:encode(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):
        return [encode(v) for v in obj]
    return obj

def one_case(p,m,q,L):
    alpha = Q(1,p)
    exponents = [p*k for k in range(m)] + [1-p*k for k in range(m)]
    require(len(set(exponents))==2*m,'distinct polynomial roots')
    v = [Q(1)]
    for k in exponents:
        v = times_linear(v, q**k)
    require(all(x*((-1)**i)>0 for i,x in enumerate(v)),'strict alternating signs')
    require(all(peval(v,q**k)==0 for k in exponents),'polynomial vanishing')
    M = sum(v[::2])
    require(M==-sum(v[1::2]) and M>0,'equal positive bank masses')
    require(2*M==peval(v,Q(-1)),'mass formula')
    n = [q**(i+1) for i in range(2*m+1)]
    u = [ni**p-1 for ni in n]
    w = [abs(vi)/M for vi in v]
    require(all(ui>0 for ui in u),'positive nodes')
    require(sum(w[::2])==sum(w[1::2])==1,'unit reference demands')
    rows = []
    contacts = []
    for s in [0,1]:
        Hder = []
        Fder = []
        Gder = []
        for k in range(m+1):
            vals = [falling(alpha,k)*ui**k*(ni**(1-p*k) if s else 1) for ni,ui in zip(n,u)]
            if k<m:
                rows.append(vals)
            Fder.append(sum(w[i]*vals[i] for i in range(0,2*m+1,2)))
            Gder.append(sum(w[i]*vals[i] for i in range(1,2*m+1,2)))
            Hder.append(Fder[-1]-Gder[-1])
        require(all(x==0 for x in Hder[:m]),f'vanishing through order {m-1}, s={s}')
        require(Hder[m]!=0,f'exact order {m}, s={s}')
        require((Hder[m]>0)==(s==0),'contraction then expansion orientation')
        require(Fder[0]==Gder[0] and Fder[1]==Gder[1],'contact load and derivative match')
        t = Fder[0]
        B = Fder[1]
        hm = Hder[m]/factorial(m)
        a = -hm/B**(m+1)
        connector = a*B/2
        K = t*m*a*a*B/2
        jump = factorial(2*m-1)*K*(1 if s==0 else -1)
        contacts.append({'s':s,'t':t,'H_derivatives':Hder,'B':B,'inverse_difference_leading_coefficient':hm,'route_gap_leading_coefficient':a,'signed_connector_leading_coefficient':connector,'aggregate_coefficient_magnitude':K,'right_minus_left_aggregate_derivative_jump':jump,'regularity_E_P':f'C^{2*m-2} not C^{2*m-1}','flow_regularity':f'C^{m-1} not C^{m}'})
    require(rank(rows)==2*m,'full confluent row rank')
    t2 = contacts[1]['t']
    require(contacts[0]['t']==1 and t2==q**(m+1),'exact second demand threshold')
    require(t2**p<p+1,'scaled optimum above last contact')
    require(q*peval(v,-q)/peval(v,Q(-1))==t2,'reciprocal polynomial threshold formula')
    require(all(v[2*m-i]==q**(i-m)*v[i] for i in range(2*m+1)),'reciprocal coefficient identity')
    edges = []
    for i,(ui,wi) in enumerate(zip(u,w)):
        bank = 'F' if i%2==0 else 'G'
        lev0 = Q(L+1 if bank=='F' else L)
        ae = lev0-1/ui
        be = 1/(ui*wi**p)
        require(ae>0 and be>0,'positive physical cost coefficients')
        require(ae+be*wi**p==lev0,'reference physical delay')
        z2 = wi*n[i]
        require(ae+be*z2**p==lev0+1,'second-boundary physical delay')
        edges.append({'index':i,'bank':bank,'reference_load':wi,'second_boundary_load':z2,'inverse_node':ui,'free_flow':ae,'congestion_coefficient':be,'power':p})
    edges.append({'index':2*m+1,'bank':'connector','reference_load':Q(0),'second_boundary_load':Q(0),'free_flow':Q(1),'congestion_coefficient':Q(1),'power':p})
    return {'p':p,'m':m,'q':q,'bank_reference_G_delay':L,'bank_reference_F_delay':L+1,'polynomial_root_exponents':exponents,'polynomial_coefficients':v,'normalising_mass':M,'physical_edges':edges,'contacts':contacts,'maximal_zero_budget':2*m,'used_zero_multiplicity':2*m,'row_rank':2*m,'all_checks_passed':True}

def main():
    cases = [one_case(4,3,Q(11,10),3),one_case(2,3,Q(11,10),5),one_case(4,5,Q(101,100),30)]
    out = {'status':'PASS','arithmetic':'fractions.Fraction only for certificate checks','cases':cases,'scope':'Finite reconstruction of explicit examples; the general zero budget and synthesis rely on the analytic proofs.'}
    (HERE/'two_contacts_exact.json').write_text(json.dumps(encode(out),indent=2)+'\n')
    for c in cases:
        print(f"p={c['p']}, m={c['m']}, q={c['q']}: PASS; edges={len(c['physical_edges'])}, thresholds=1,{c['contacts'][1]['t']}, rank={c['row_rank']}")
        for j,x in enumerate(c['contacts']):
            print(f"  contact {j+1}: B={float(x['B']):.12g}, a={float(x['route_gap_leading_coefficient']):.12g}, derivative jump={float(x['right_minus_left_aggregate_derivative_jump']):.12g}")
    print('All exact checks passed.')

if __name__=='__main__':
    main()
