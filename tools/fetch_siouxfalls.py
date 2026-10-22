"""Fetch the Sioux Falls test-network files used by code/sioux_falls.py and check their SHA-256.

The files are not redistributed in this package: the Transportation Networks for Research
repository states that its data sets are for academic research purposes only and must be
cited (https://github.com/bstabler/TransportationNetworks, README, "License"). This script
downloads them from a pinned commit into data/siouxfalls/ and refuses any file whose
SHA-256 differs from the copy used for the published results.

Usage (from the package root):  python3 tools/fetch_siouxfalls.py
Citation: Transportation Networks for Research Core Team, Transportation Networks for Research,
https://github.com/bstabler/TransportationNetworks (accessed 22 October 2026); network from
LeBlanc, Morlok and Pierskalla, Transportation Research 9 (1975) 309-318.
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

COMMIT = "d1639b4ef218c17928ba573e806ddf8ba5e7ae6d"
BASE = f"https://raw.githubusercontent.com/bstabler/TransportationNetworks/{COMMIT}/SiouxFalls/"
FILES = {
    "SiouxFalls_net.tntp": "ace99b24cec69c273ff0cf3d6d074110177f0cc0ae24b0c7a9f4f4cb5e27635c",
    "SiouxFalls_trips.tntp": "56f9566857f3f66730fd5c4232258d7ee3ac2931a476526331afd062f4958de7",
}


def main():
    out = Path(__file__).resolve().parent.parent / "data" / "siouxfalls"
    out.mkdir(parents=True, exist_ok=True)
    for name, digest in FILES.items():
        target = out / name
        if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest() == digest:
            print(f"ok    {name} (already present)")
            continue
        req = urllib.request.Request(BASE + name, headers={"User-Agent": "wardrop-smooth-route-regimes replay"})
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read()
        got = hashlib.sha256(data).hexdigest()
        if got != digest:
            sys.exit(f"SHA-256 mismatch for {name}: expected {digest}, got {got}")
        target.write_bytes(data)
        print(f"ok    {name} (downloaded, SHA-256 verified)")


if __name__ == "__main__":
    main()
