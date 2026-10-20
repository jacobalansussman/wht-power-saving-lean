#!/usr/bin/env python3
"""extract_lists.py : read the two data lists of Braverman and He's causal code out of a plain-text copy of their
preprint, and write them as JSON for the other programs in this folder.

The preprint is M. Braverman and Z. He, "The Multiple Unicast Conjecture is False", arXiv:2610.10108, version 1
(7 October 2026). Its appendix lists (1) the 175 deleted edge indices R and (2) the vertex-processing order pi of
the 182 vertices. Neither list is copied into this repository; make your own text copy of the paper, for example

    pdftotext -layout 2610.10108v1.pdf bh.txt        (pdftotext is part of poppler)
    python3 extract_lists.py bh.txt bh_lists.json

The program looks for the sentence that introduces each list, keeps only lines that consist of comma-separated
whole numbers (so page numbers and headings are skipped), and checks what the preprint says about the lists:
175 distinct indices below 910, and a permutation of 0..181. If a check fails, it says so and writes nothing.
"""
import json, re, sys

if len(sys.argv) != 3:
    sys.exit("usage: python3 extract_lists.py <text copy of arXiv:2610.10108v1> <output.json>")
text = open(sys.argv[1], encoding="utf-8", errors="replace").read().split("\n")

def find(pattern, start=0):
    for k in range(start, len(text)):
        if re.search(pattern, text[k], flags=re.I):
            return k
    sys.exit("could not find a line matching %r; is this the version 1 text?" % pattern)

def numbers(a, b):
    # a list line: whole numbers separated by commas, possibly with a leading "R =" or "pi =", a final full stop and
    # the bracket glyphs a PDF-to-text tool leaves (characters of the private-use area); page numbers have no comma
    out = []
    for line in text[a:b]:
        s = re.sub(r"[\ue000-\uf8ff]", " ", line)
        s = re.sub(r"^\s*(R|\u03c0|pi)\s*=", " ", s).replace(".", " ")
        if "," in s and re.fullmatch(r"[\d,\s]+", s):
            out += [int(x) for x in re.findall(r"\d+", s)]
    return out

r0 = find(r"set of 175 deleted edge indices")
p0 = find(r"vertex-processing order", r0 + 1)
p1 = find(r"elementary proof of the dimension bound|^\s*B\s*$", p0 + 1)
R, PI = numbers(r0 + 1, p0), numbers(p0 + 1, p1)
ok = True
if not (len(R) == 175 and len(set(R)) == 175 and max(R) < 910):
    print("R: expected 175 distinct indices below 910, found", len(R), "numbers,", len(set(R)), "distinct"); ok = False
if sorted(PI) != list(range(182)):
    print("pi: expected a permutation of 0..181, found", len(PI), "numbers"); ok = False
if not ok:
    sys.exit("not written")
json.dump({"source": "Braverman and He, arXiv:2610.10108v1, appendix: deleted edges R and vertex order pi",
           "deleted_edges": R, "vertex_order": PI}, open(sys.argv[2], "w"))
print("wrote", sys.argv[2], ": 175 deleted edges, vertex order of 182 vertices")
