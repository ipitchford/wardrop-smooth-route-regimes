# Smoothness within and across Wardrop route regimes

Reproducibility package for the paper *Smoothness within and across Wardrop route regimes* (Anonymous, 22 October 2026). Version 1.0.0-candidate; DOI [10.5281/zenodo.23552630](https://doi.org/10.5281/zenodo.23552630). Release records: [STATUS.md](STATUS.md), [ASSURANCE.md](ASSURANCE.md), [PROVENANCE.md](PROVENANCE.md), [AI_INDEX.md](AI_INDEX.md), [LICENSES.md](LICENSES.md).

## In plain English

In a traffic assignment every driver takes a cheapest route. As demand grows, new routes become cheapest from time to time. Cominetti, Dose and Scarsini (Transportation Research Part B, 2024; arXiv:2305.03459v2) asked two questions about this.

1. **While the set of cheapest routes stays the same, do equilibrium link flows change smoothly with demand, even when link costs are flat at zero load, as standard BPR costs are?** Yes, at every interior demand level, on any network. At zero demand a general version fails; whether it fails for BPR costs is open.
2. **When a new route becomes cheapest, must the slope of total travel time, or of the price of anarchy, jump?** No. A five-link network with ordinary quartic BPR costs has a genuine new route at which total travel time and the price of anarchy have four continuous derivatives. Five links are the fewest possible with a common BPR exponent.

The paper also shows that the slope never jumps the wrong way, even with zero free-flow times, and that for almost every choice of free-flow times it does jump. Smooth route changes are therefore possible but exceptional, at least at breakpoints where every link on the old cheapest routes carries flow. When they happen, the new route carries so little flow that an assignment run at ordinary convergence tolerances cannot see it: after 10% more demand the new route carries about 0.01% of the traffic.

## What kind of evidence supports each result

| Result | Where | Status |
|---|---|---|
| Fixed-regime $C^1$ link flows for costs with positive slope at positive load (Question 1) | Thm 3.1 | proof |
| General endpoint extension fails; BPR endpoint case open | Prop 3.3 | proof; open |
| Slope comparison of CDS Thm 4.4(i) for BPR costs, expansions and contractions, zero free-flow allowed | Thm 4.1, Cor 4.2 | proof |
| Strict slope jumps for almost every free-flow vector (nondegenerate breakpoints) | Thm 4.5 | proof |
| Three-link examples: strictly convex costs; mixed BPR exponents (contraction) | §5.1 | proof; exact check |
| Five-link quartic BPR expansion and contraction with $E,P\in C^4\setminus C^5$ | §5.2 | proof; exact check |
| $m+2$ links give $C^{2m-2}\setminus C^{2m-1}$ for every $p>1$ | Thm 5.3 | proof |
| Five links minimal for a common exponent (expansions and contractions) | Thm 6.1, App. B | proof; finite enumeration as a check |
| Mixed exponents: contractions need 3 links; expansions 4 or 5 | Table 3 | proof; open |
| Optimum-side transition at $t=1/\gamma$: $S\in C^{2m-1}\setminus C^{2m}$; same for $P$ when $E$ is $C^{2m}$ there (checked exactly for the quartic example); jump ≈ +89.7 | Thm 7.4 | proof; exact interval check; numerical value |
| Perturbation, demand-ratio and transverse-direction results | §7 | proof |
| Hidden activation below tolerance, near-match kinks, optimum values | §8 | numerical (mpmath, 60 digits) |
| Sioux Falls scan: activations are transverse with visible kinks | §8.5 | exploratory |
| Several contacts: sharp budget and synthesis in the two-bank architecture | App. A | proof; exact check |

"Proof" means a complete analytic proof in the paper. "Exact" means rational arithmetic that confirms constants; it never replaces a proof. "Numerical" means floating-point computation, which is not a certificate. "Exploratory" means an illustration with no claim of generality.

Scope reminders:
- The five-link minimality holds among costs $a_e+b_ez^p$ with $a_e,b_e>0$ and one common $p>1$, with every simple path available, proportional positive demands, and nested breakpoints.
- The order-dependent count $2+\sum m_j$ is specific to the two-bank architecture.
- Theorem 4.5 covers only nondegenerate breakpoints.

## Replay everything

```bash
python3 -B -O negative_controls.py  # standard library: exact jets of the five-link example; corruptions rejected
python3 verify_all.py --exact        # exact programs only
python3 tools/fetch_siouxfalls.py    # third-party Sioux Falls data (not redistributed; SHA-256 checked)
python3 verify_all.py                # about 5 minutes; all checks
python3 verify_all.py --quick        # skips the two slow numerical programs
```

`verify_all.py` copies the package to a temporary directory, reruns every program and compares every generated record byte for byte with the copy supplied here. It changes nothing in the package. Requirements: Python 3.11 or later. The exact programs use the standard library; `verify_small_examples.py` and `verify_optimum_transition.py` also need `sympy`. The numerical programs need `mpmath`, `numpy`, `scipy` and `matplotlib` (tested versions in `requirements-numerics.txt`; the replay was run with Python 3.13).

| Program | Output | Type |
|---|---|---|
| `code/verify_exact.py` | `evidence/exact_verification.json` | exact |
| `code/verify_small_examples.py` | `evidence/small_examples.json` | exact |
| `code/verify_optimum_transition.py` | `evidence/optimum_transition.json` | exact; numerical |
| `code/verify_practice.py` | `evidence/practice.json` | numerical |
| `checks/moments/verify_two_contacts.py` | `checks/moments/two_contacts_exact.json` | exact |
| `checks/transversality/verify_multicontact_reconstruction.py` | `checks/transversality/multicontact_reconstruction.json` | exact |
| `checks/transversality/verify_transversality.py` | `checks/transversality/exact_checks.json` | exact |
| `checks/small_networks/verify_small_network_identities.py` | `checks/small_networks/small_network_verification.json` | exact |
| `checks/small_networks/verify_weighted_basis.py` | `checks/small_networks/weighted_basis_verification.json` | exact |
| `checks/small_networks/verify_minimality.py` | `evidence/minimality_checks.json` | exact |
| `code/numerical_demo.py`, `code/plot_results.py` | `evidence/numerical_*`, `evidence/smooth_transition_*`, `figures/smooth_transition.*` | numerical |
| `code/sioux_falls.py` | `evidence/sioux_falls_results.json`, `figures/sioux_falls.pdf` | exploratory |

The Sioux Falls files are **not redistributed**: the Transportation Networks for Research repository (https://github.com/bstabler/TransportationNetworks; LeBlanc, Morlok and Pierskalla 1975) states that its data sets are for academic research purposes only and must be cited. `tools/fetch_siouxfalls.py` downloads them from commit d1639b4ef218c17928ba573e806ddf8ba5e7ae6d into `data/siouxfalls/` and refuses any file whose SHA-256 differs: `SiouxFalls_net.tntp` ace99b24cec69c273ff0cf3d6d074110177f0cc0ae24b0c7a9f4f4cb5e27635c, `SiouxFalls_trips.tntp` 56f9566857f3f66730fd5c4232258d7ee3ac2931a476526331afd062f4958de7.

## Build the paper

```bash
cd paper
pdflatex -interaction=nonstopmode wardrop_smooth_transitions.tex
pdflatex -interaction=nonstopmode wardrop_smooth_transitions.tex
```

The bibliography is in `paper/bibliography.tex`, generated from `paper/references.bib` by `python3 tools/make_bibliography.py`, which includes only cited entries. No BibTeX run is needed. Every reference was checked against Crossref or arXiv on 22 October 2026. PDF bytes can differ between builds because of build metadata.

## Files

- `paper/` — LaTeX source (`wardrop_smooth_transitions.tex`, `sections/`), bibliography and the built PDF.
- `code/`, `checks/` — verification programs (see table above).
- `evidence/` — generated records and `claim_ledger.json`, which lists each claim with its hypotheses, status and evidence.
- `figures/` — figures; `data/siouxfalls/` — where `tools/fetch_siouxfalls.py` places the third-party test-network data.
- `negative_controls.py` — semantic negative controls; `tools/` — bibliography generator and data fetcher.
- Release records: `STATUS.md`, `ASSURANCE.md`, `PROVENANCE.md`, `AI_INDEX.md`, `CLAIMS.json`, `RESEARCH_GATES.json`, `REVIEW_SUMMARY.md`, `CITATION.cff`, `.zenodo.json`, `LICENSES.md` (CC0-1.0 prose and original data, MIT code).
- `verify_all.py`, `SHA256SUMS` — replay script and checksums of every file in the package.

## Status

This is an unrefereed research paper. AI tools were used to help develop and check proofs, write code, check the literature and edit the text; the author takes responsibility for the content. No empirical data were collected.
