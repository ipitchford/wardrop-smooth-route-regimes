"""Convert paper/references.bib into paper/bibliography.tex (a thebibliography block).

Standard library only; deterministic output. Run from the package root.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def parse(text):
    entries = []
    for m in re.finditer(r"@(\w+)\{([^,]+),(.*?)\n\}", text, re.S):
        kind, key, body = m.group(1).lower(), m.group(2).strip(), m.group(3)
        fields = {}
        for f in re.finditer(r"(\w+)\s*=\s*\{(.*?)\}\s*,?\s*\n", body + "\n", re.S):
            fields[f.group(1).lower()] = " ".join(f.group(2).split())
        entries.append((kind, key, fields))
    return entries

def authors(a):
    names = [n.strip() for n in a.split(" and ")]
    out = []
    for n in names:
        if n.startswith("{"):
            out.append(n.strip("{}"))
            continue
        if "," in n:
            last, first = [x.strip() for x in n.split(",", 1)]
        else:
            parts = n.split(); last, first = parts[-1], " ".join(parts[:-1])
        initials = " ".join(p[0] + "." for p in re.split(r"[\s]+", first) if p)
        out.append(f"{last}, {initials}".strip(", "))
    if len(out) == 1:
        return out[0]
    return ", ".join(out[:-1]) + ", \\& " + out[-1]

def clean(t):
    return t.replace("{", "").replace("}", "") if t else t

def fmt(kind, key, f):
    s = authors(f.get("author", "")) + f" ({f.get('year') or 'n.d.'}). "
    title = clean(f.get("title", ""))
    if kind == "article":
        s += f"{title}. \\emph{{{f.get('journal','')}}}"
        if f.get("volume"): s += f", {f['volume']}"
        if f.get("number"): s += f"({f['number']})"
        if f.get("pages"): s += f", {f['pages']}"
        s += "."
    elif kind == "inproceedings":
        s += f"{title}. In \\emph{{{f.get('booktitle','')}}}"
        if f.get("series"): s += f", {f['series']}"
        if f.get("volume"): s += f" {f['volume']}"
        if f.get("pages"): s += f", pp. {f['pages']}"
        s += f". {f.get('publisher','')}."
    elif kind == "book":
        s += f"\\emph{{{title}}}"
        if f.get("volume"): s += f", Vol. {f['volume']}"
        s += f". {f.get('publisher','')}."
    else:
        s += f"\\emph{{{title}}}."
        if f.get("howpublished"): s += " " + f["howpublished"] + "."
    if f.get("note"): s += " " + f["note"] + "."
    if f.get("doi"): s += f" \\url{{https://doi.org/{f['doi']}}}."
    elif f.get("eprint") and kind == "misc" and "arXiv" not in f.get("note",""):
        s += f" arXiv:{f['eprint']}."
    return f"\\bibitem{{{key}}}\n{s}\n"

def main():
    entries = parse((ROOT / "paper" / "references.bib").read_text(encoding="utf-8"))
    tex = "".join(p.read_text(encoding="utf-8") for p in sorted((ROOT / "paper").rglob("*.tex")) if p.name != "bibliography.tex")
    cited = set()
    for m in re.finditer(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", tex):
        cited.update(k.strip() for k in m.group(1).split(","))
    missing = cited - {e[1] for e in entries}
    if missing:
        raise SystemExit(f"cited but not in references.bib: {sorted(missing)}")
    entries = [e for e in entries if e[1] in cited]
    entries.sort(key=lambda e: (e[2].get("author", "").lower(), e[2].get("year", "")))
    body = "\\begin{thebibliography}{99}\n\n" + "\n".join(fmt(*e) for e in entries) + "\n\\end{thebibliography}\n"
    (ROOT / "paper" / "bibliography.tex").write_text(body, encoding="utf-8")
    print(f"wrote {len(entries)} entries")

if __name__ == "__main__":
    main()
