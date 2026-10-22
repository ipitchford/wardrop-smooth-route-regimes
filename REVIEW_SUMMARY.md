# Review summary (rating-free)

This candidate was checked internally before release. The checks were AI-assisted: a first round of internal reviews in five roles (editor, methodology, domain, practitioner perspective and devil's advocate), an AI-assisted adversarial review by a different model (OpenAI gpt-5.6-sol; a first review and a targeted audit of the revision) and a verification re-review. **None of this is peer review, independent reproduction or formal verification.** The original reports are kept privately; this file records their substantive objections and how they were handled. It contains no scores or ratings.

## Main objections and responses

| Objection | Response in this version |
|---|---|
| Definitions of feasible sets, regimes and route sets were not stated precisely enough for the theorems to be checked | Section 2 defines the feasible set and feasible-load set (same route sets for equilibrium and optimum) and gives a decomposition-independent definition of the regime, valid at zero demand |
| The price-of-anarchy part of Theorem 3.1 depended on a lemma of the source paper whose hypotheses were not checked | A self-contained envelope proposition (Danskin form) was added. The adversarial audit found that a stated specialisation (S is C¹ whenever z·c(z) is convex and c is strictly increasing) did not follow, with an explicit counterexample to uniqueness of the optimum; the hypothesis and its uses were corrected. BPR costs are unaffected |
| Scope of the five-link minimality was stated too broadly | The class (a + b z^p, a, b > 0, common p > 1, all simple paths, proportional positive demands, nested breakpoints) is stated wherever minimality is claimed; the derivative count 2 + Σ m_j is called specific to the two-bank architecture |
| The devil's advocate asked whether fewer links suffice outside the common-exponent class | Three-link examples were added: strictly convex polynomial costs (expansion) and mixed BPR exponents (contraction), verified exactly |
| The four-edge exclusion was a sketch | Rewritten as Appendix B with a case table and a separate exact enumeration (`checks/small_networks/verify_minimality.py`); the re-review audited it line by line and found no defect |
| Positioning against recent and classical work | Comparisons added with classical sensitivity analysis (Hall; Dafermos and Nagurney; Tobin and Friesz; Qiu and Magnanti; Patriksson; Lu; Josefsson and Patriksson) and with Wang, Shao and Wu (2026), whose C¹ result for total cost and the price of anarchy overlaps but does not give C¹ link flows; D'Ambrosio (2026) is distinguished. Every reference was checked against Crossref or arXiv |
| Questions must be quoted faithfully and answered as posed | Both questions are quoted verbatim from arXiv v2 with their locations; Question 1 is answered "yes at interior points", with the endpoint stated separately |
| Readability for non-specialists | New plain abstract, a "main idea" section with a figure, a results-at-a-glance table with status labels, a notation dictionary and a conclusion with open problems |
| What the results mean in practice | Section 8 added: hidden activation below assignment tolerance, link flows reveal the change earlier, near-match kinks computed from the exact formula, sign of the cost change (no Braess effect at the expansion), the optimum-side transition, and an exploratory Sioux Falls scan |
| Extensions requested | Genericity of strict jumps proved at nondegenerate breakpoints (degenerate case open); weak slope comparison for BPR costs proved for expansions and contractions; optimum-side transition proved; contractions also need five links with a common exponent |

## Corrections found during revision

- The five-link example has further breakpoints at low demand (connector used for t in (0.2105, 0.5839); one bank link unused below t ≈ 0.631). All statements about the example are local to t = 1 and t = 1/γ.
- A near-match figure based on the asymptotic law was replaced by values from the exact formula, because the asymptotic law is inaccurate at those perturbation sizes.

## Open after review

The BPR endpoint case at zero demand; genericity at degenerate breakpoints; four versus five links for mixed-exponent expansions; non-nested breakpoints; prevalence of near-tangent activations in calibrated networks.
