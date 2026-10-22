# Assurance

| Dimension | State | Note |
|---|---|---|
| Written proofs | present | Complete proofs of every theorem in the manuscript; not machine-checked |
| Exact checks | passed | Rational and symbolic programs listed in [README.md](README.md) (standard library; SymPy for two programs); they confirm constants and instances, not the general theorems |
| Producer replay | passed | `verify_all.py` reruns every program in a temporary copy and compares every generated record byte for byte. A clean full replay of the final research package ran in a separate copy; the released programs and records were bound to it by SHA-256, and checks affected by publication changes were rerun |
| Negative controls | passed | `negative_controls.py` (standard library, exact): recomputes the bank-delay jets of the five-link example from the physical costs and confirms that a perturbed connector free-flow time, perturbed bank coefficients or free-flow times, and a wrong claimed gap coefficient are rejected |
| Independent rerun | not assessed | No unaffiliated rerun |
| Independent reimplementation | not assessed | Exact and numerical checks are producer-side; the negative-control jet computation is a second implementation inside the same workflow, not external reproduction |
| Formal verification | not assessed | No proof assistant |
| Specialist review | not assessed | AI-assisted internal and adversarial review is not external specialist review |
| Editorial peer review | not assessed | None |
| Prior art | bounded | Dated bounded searches summarised in [RESEARCH_GATES.json](RESEARCH_GATES.json); no priority claim |

Trust boundary: the general statements rest on the written proofs. Exact programs check coefficients and instances. Values labelled numerical are deterministic high-precision floating-point computations, not enclosures. The Sioux Falls scan is exploratory and uses third-party data fetched by `tools/fetch_siouxfalls.py` and checked by SHA-256.
