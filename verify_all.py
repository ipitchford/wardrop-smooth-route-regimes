"""Replay every verification program and compare its records byte for byte.

Usage (from the package root):
    python3 verify_all.py            # all checks (about 5 minutes)
    python3 verify_all.py --quick    # skip the two slow numerical programs
    python3 verify_all.py --exact    # only the exact programs (used by CI)

The Sioux Falls step needs the third-party test-network files, which are not
redistributed: fetch them first with `python3 tools/fetch_siouxfalls.py`.

The package is copied to a temporary directory, every program is run there,
and each generated record is compared with the supplied copy. Nothing in the
package is modified. Exit status is nonzero if any program fails or any
record differs. Exact programs need only the Python standard library (plus
sympy for verify_small_examples.py and verify_optimum_transition.py);
numerical programs need mpmath, numpy, scipy and matplotlib
(see requirements-numerics.txt).
"""
import argparse
import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# (label, command relative to package root, records produced, slow?)
STEPS = [
    ("exact family series", ["code/verify_exact.py", "--output", "evidence/exact_verification.json"],
     ["evidence/exact_verification.json"], False),
    ("small examples", ["code/verify_small_examples.py", "--output", "evidence/small_examples.json"],
     ["evidence/small_examples.json"], False),
    ("optimum-side transition", ["code/verify_optimum_transition.py", "--output", "evidence/optimum_transition.json"],
     ["evidence/optimum_transition.json"], False),
    ("two contacts", ["checks/moments/verify_two_contacts.py"],
     ["checks/moments/two_contacts_exact.json"], False),
    ("multicontact reconstruction", ["checks/transversality/verify_multicontact_reconstruction.py"],
     ["checks/transversality/multicontact_reconstruction.json"], False),
    ("transversality and codimension", ["checks/transversality/verify_transversality.py"],
     ["checks/transversality/exact_checks.json"], False),
    ("small-network identities", ["checks/small_networks/verify_small_network_identities.py"],
     ["checks/small_networks/small_network_verification.json"], False),
    ("weighted least squares", ["checks/small_networks/verify_weighted_basis.py"],
     ["checks/small_networks/weighted_basis_verification.json"], False),
    ("minimality route spaces", ["checks/small_networks/verify_minimality.py", "--output", "evidence/minimality_checks.json"],
     ["evidence/minimality_checks.json"], False),
    ("five-link numerics", ["code/numerical_demo.py", "--output-dir", "evidence"],
     ["evidence/numerical_sanity.json", "evidence/numerical_curve.csv"], False),
    ("five-link figure", ["code/plot_results.py", "--output-dir", "figures", "--evidence-dir", "evidence"],
     ["evidence/smooth_transition_flow.csv", "evidence/smooth_transition_cost.csv",
      "evidence/smooth_transition_metadata.json", "figures/smooth_transition.pdf",
      "figures/smooth_transition.png"], False),
    ("practice numbers", ["code/verify_practice.py", "--output", "evidence/practice.json"],
     ["evidence/practice.json"], True),
    ("Sioux Falls (exploratory)", ["code/sioux_falls.py", "data/siouxfalls/SiouxFalls_net.tntp",
                                   "data/siouxfalls/SiouxFalls_trips.tntp",
                                   "evidence/sioux_falls_results.json", "figures/sioux_falls.pdf"],
     ["evidence/sioux_falls_results.json", "figures/sioux_falls.pdf"], True),
]

EXACT = {"exact family series", "small examples", "optimum-side transition", "two contacts",
         "multicontact reconstruction", "transversality and codimension", "small-network identities",
         "weighted least squares", "minimality route spaces"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quick", action="store_true", help="skip the slow numerical programs")
    parser.add_argument("--exact", action="store_true", help="run only the exact programs")
    args = parser.parse_args()
    failures = []
    with tempfile.TemporaryDirectory(prefix="wardrop_replay_") as tmp:
        work = Path(tmp) / "package"
        shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns("paper", "__pycache__", "*.pyc"))
        for label, cmd, records, slow in STEPS:
            if slow and (args.quick or args.exact):
                print(f"SKIP  {label}")
                continue
            if args.exact and label not in EXACT:
                print(f"SKIP  {label} (not exact)")
                continue
            missing = [c for c in cmd if c.startswith("data/") and not (work / c).exists()]
            if missing:
                failures.append(f"{label}: missing {missing}; run python3 tools/fetch_siouxfalls.py")
                print(f"FAIL  {label}: missing third-party data; run python3 tools/fetch_siouxfalls.py first")
                continue
            for rec in records:
                target = work / rec
                if target.exists():
                    target.unlink()
            proc = subprocess.run([sys.executable, *cmd], cwd=work, capture_output=True, text=True)
            if proc.returncode != 0:
                failures.append(f"{label}: exit {proc.returncode}")
                print(f"FAIL  {label} (exit {proc.returncode})\n{proc.stderr[-2000:]}")
                continue
            bad = []
            for rec in records:
                new, old = work / rec, ROOT / rec
                if not new.exists():
                    bad.append(f"{rec} not produced")
                elif not old.exists():
                    bad.append(f"{rec} missing from package")
                elif new.read_bytes() != old.read_bytes():
                    bad.append(f"{rec} differs (new {sha(new)[:12]}, supplied {sha(old)[:12]})")
            if bad:
                failures.extend(f"{label}: {b}" for b in bad)
                print(f"DIFF  {label}: " + "; ".join(bad))
            else:
                print(f"PASS  {label} ({len(records)} record(s) byte-identical)")
    if failures:
        print(f"\n{len(failures)} problem(s).")
        sys.exit(1)
    print("\nAll replayed records are byte-identical to the supplied copies.")


if __name__ == "__main__":
    main()
