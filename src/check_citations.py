"""Citation-integrity check for the DNA data-storage paper.

Asserts every entry in paper/references.bib carries the fields a real citation
needs (author, title, year, and a venue + a DOI so it can be verified), and, once
a manuscript exists, that every \\cite key in paper/*.tex is defined in the bib.

Run: python src/check_citations.py
"""
import glob
import os
import re

ROOT = os.path.join(os.path.dirname(__file__), os.pardir)
BIB = os.path.join(ROOT, "paper", "references.bib")
REQUIRED = ("author", "title", "year", "doi")


def parse_bib(text):
    """Return {key: {field: value}} for each @article-style entry."""
    entries = {}
    for m in re.finditer(r"@\w+\s*\{\s*([^,]+),(.*?)\n\}", text, re.DOTALL):
        key, body = m.group(1).strip(), m.group(2)
        fields = dict(
            (f.lower(), v.strip())
            for f, v in re.findall(r"(\w+)\s*=\s*\{(.*?)\}", body, re.DOTALL)
        )
        entries[key] = fields
    return entries


def cite_keys_in_tex():
    keys = set()
    for tex in glob.glob(os.path.join(ROOT, "paper", "*.tex")):
        with open(tex, encoding="utf-8") as f:
            for m in re.finditer(r"\\cite[a-z]*\{([^}]+)\}", f.read()):
                keys.update(k.strip() for k in m.group(1).split(","))
    return keys


def check():
    with open(BIB, encoding="utf-8") as f:
        entries = parse_bib(f.read())
    assert entries, "no entries parsed from references.bib"

    for key, fields in entries.items():
        missing = [r for r in REQUIRED if r not in fields or not fields[r]]
        assert not missing, f"{key}: missing required field(s): {missing}"
        venue = fields.get("journal") or fields.get("booktitle") or fields.get("publisher")
        assert venue, f"{key}: no venue (journal/booktitle/publisher)"

    undefined = cite_keys_in_tex() - set(entries)
    assert not undefined, f"\\cite keys not in references.bib: {sorted(undefined)}"

    print(f"ok  {len(entries)} bib entries, all with {', '.join(REQUIRED)} + venue")
    print(f"ok  {len(cite_keys_in_tex())} \\cite keys, all defined")
    return entries


if __name__ == "__main__":
    check()
