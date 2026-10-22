# Agent navigation

## Identity

- Release slug: `wardrop-smooth-route-regimes`; version `1.0.0-candidate`; DOI 10.5281/zenodo.23552630 (concept DOI 10.5281/zenodo.23552629)
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md))
- Manuscript: [PDF](paper/wardrop_smooth_transitions.pdf), [LaTeX source](paper/wardrop_smooth_transitions.tex) with [sections](paper/sections/) and the generated [bibliography](paper/bibliography.tex)
- Machine-readable claims: [CLAIMS.json](CLAIMS.json) and the per-claim ledger [evidence/claim_ledger.json](evidence/claim_ledger.json)
- Assurance limits: [ASSURANCE.md](ASSURANCE.md); internal review summary: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)
- Provenance: [PROVENANCE.md](PROVENANCE.md); prior art and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json)
- Citation: [CITATION.cff](CITATION.cff), [Zenodo metadata](.zenodo.json); licences: [LICENSES.md](LICENSES.md)

## Scientific map

Object: static nonatomic multi-commodity Wardrop equilibrium with separable nondecreasing link costs and demands t·r. E(t) is total equilibrium cost, S(t) optimum cost, P = E/S the price of anarchy. A regime is the set of cheapest routes; a breakpoint is where it changes. Source: Cominetti, Dose and Scarsini (CDS), Transportation Research Part B 182, 102922 (2024), arXiv:2305.03459v2 numbering.

- **CDS Question 1 (after Remark 3.5), answered yes at interior points** (Theorem 3.1, proof): C¹ costs with c'(z) > 0 for z > 0 give C¹ link flows on every open interval with a fixed regime, any network. Endpoint t0 = 0: general extension false (Proposition 3.3); BPR case open.
- **CDS Question 2 (before Example 4.9), answered no**: five links, costs F0 = 2+z⁴, F2 = 8/3+z⁴/243 (s→v), G1 = 3/2+z⁴/162, G3 = 7/4+z⁴/4 (u→v), connector 1+z⁴ (s→u), demands 4t: expansion at t = 1 with E, P ∈ C⁴ \ C⁵ (Section 5.2, proof; exact). m+2 links: C^(2m−2) \ C^(2m−1) (Theorem 5.3). Three links suffice with convex or mixed-exponent costs (Section 5.1).
- **Minimality** (Theorem 6.1, Appendix B, proof): five links minimal for a + b z^p, a, b > 0, common p > 1, all simple paths, proportional positive demands, nested breakpoints.
- **Weak comparison** (Theorem 4.1, Corollary 4.2): CDS Theorem 4.4(i) holds for BPR costs, expansions and contractions, zero free-flow allowed. **Genericity** (Theorem 4.5): strict jumps for almost every free-flow vector at nondegenerate nested breakpoints.
- **Practice** (Section 8, numerical): hidden activation below assignment tolerance; link flows are only C² at the switch. **Sioux Falls** (Section 8.5, exploratory): transverse activations with visible kinks.

Do not read: the minimality outside its class; "exact" as a proof of the general theorems; numerical values as enclosures; the Sioux Falls scan as evidence about prevalence; Theorem 3.1 as covering zero demand.

## Replay

Python 3.11+; exact programs need the standard library (SymPy for two); numerical programs need [requirements-numerics.txt](requirements-numerics.txt) plus mpmath:

```sh
shasum -a 256 -c SHA256SUMS
python3 -B -O negative_controls.py     # exact jets of the five-link example; corruptions rejected
python3 -B verify_all.py --exact       # exact programs, every record byte for byte
python3 tools/fetch_siouxfalls.py      # third-party data for the Sioux Falls step (SHA-256 checked)
python3 -B verify_all.py               # everything, about 5 minutes
```

[verify_all.py](verify_all.py) reruns every program in a temporary copy and compares outputs byte for byte. Programs: [code/](code/), [checks/](checks/); records: [evidence/](evidence/); figures: [figures/](figures/); data fetcher: [tools/fetch_siouxfalls.py](tools/fetch_siouxfalls.py).

## Useful next research question

Decide whether CDS Theorem 3.2 extends to BPR costs at zero demand, or whether four links suffice for an equal-slope expansion when BPR exponents differ between links.
