#!/usr/bin/env python3
"""summ.py (pk2-orbits): tally orb.py trial lines (final verdicts of each design) from the out_*.txt files given.
Final status of a design (its last evaluation), per the exact verdicts:
  FOUND            : every output exact YES (re-verified)       -- none expected
  one output exact : some output exact YES
  dead (all c)     : linearised psi or chart relaxation refuses some output: no choice of orbit scalars works
  exact no         : exact end solve NO for some output at the drawn scalars (full solve or a layered screen)
  < W top states   : fewer than W states reach level m: both outputs impossible for every scalar choice
  cap              : neither; every failing output exceeded the closure cap at every layer tried
"""
import sys, ast, re, collections

def classify(h):
    ps, ch, ex = h[0], h[1], h[2]
    top = h[5] if len(h) > 5 else ""
    W = len(ps)
    if all(x == "Y" for x in ex): return "FOUND"
    if any(x == "Y" for x in ex): return "one output exact"
    if "0" in ps or "n" in ch: return "dead for every scalar choice (relaxations)"
    if isinstance(top, str) and top.startswith("T"): return "fewer than W states reach level m (dead for every scalar)"
    if "n" in ex: return "exact no at drawn scalars (full solve)"
    if len(h) > 6 and "0" in h[6]: return "exact no at drawn scalars (psi at those scalars)"
    if isinstance(top, str) and "n" in top: return "exact no at drawn scalars (layered screen)"
    return "cap (undecided)"

def main(files):
    grand = collections.Counter(); gn = 0
    for fn in files:
        T = collections.Counter(); n = 0; kinds = collections.Counter(); found = 0; H = collections.Counter()
        layer = collections.Counter(); norb = 0; best = [None]
        for line in open(fn):
            if line.startswith("FRAME-GAME PROGRAM FOUND"): found += 1
            if not line.startswith("trial "): continue
            m = re.match(r"trial (\d+): a=(.*?) dup=(\d+) orbits=(\[.*?\]) hist=(\[.*\])$", line.strip())
            if not m: continue
            orbits = ast.literal_eval(m.group(4)); hist = ast.literal_eval(m.group(5))
            n += 1; h = hist[-1]
            T[classify(h)] += 1
            H[h[3]] += 1
            if "0" in h[0]: H["psilin refuses"] += 1
            # layered depth: per output, the layer at which it stopped (n = fails exactly, c = capped); smaller = deeper
            if len(h) > 5 and isinstance(h[5], str) and h[5].startswith("L") or (len(h) > 5 and "," in str(h[5])):
                parts = str(h[5]).split(",")
                if (all(x.startswith("L") and x.endswith("c") for x in parts) and "n" not in h[2]
                        and not (len(h) > 6 and "0" in h[6])
                        and "0" not in h[0] and "n" not in h[1]):
                    worst = max(int(x[1]) for x in parts)      # every output passed the layers above 'worst' exactly
                    key = (worst, str(h[5]), h[2], line[:60].strip())
                    if best[0] is None or worst < best[0][0]: best[0] = key
            if "0" not in h[0]:
                H["psi-pass"] += 1
                if h[3] == "H": H["psi-pass and height system"] += 1
            if len(h) > 5 and h[5] not in ("-", "-,-", "-,-,-"):
                for x in str(h[5]).split(","):
                    if x != "-": layer[x] += 1
            for o in orbits: kinds[o[0]] += 1
            norb += len(orbits)
        print(f"== {fn}: {n} designs, {norb} orbit injections; FOUND lines: {found}")
        for k, v in sorted(T.items()): print(f"   {k}: {v}")
        print(f"   heightfilter at the end: height system exists (cannot save) {H.get('H', 0)}, none {H.get('-', 0)}; "
              f"of the {H.get('psi-pass', 0)} that pass linearised psi for every output, {H.get('psi-pass and height system', 0)} have a height system")
        if layer: print(f"   layered-screen verdicts on capped outputs: {dict(sorted(layer.items()))}")
        print(f"   linearised psi refuses some output: {H.get('psilin refuses', 0)} of {n}")
        if best[0]: print(f"   deepest undecided design (no screen refuses; every output passes the layers above d, then caps at d): (d, ...): {best[0]}")
        print(f"   orbit kinds: {dict(kinds)}")
        grand.update(T); gn += n
    print(f"== ALL: {gn} designs")
    for k, v in sorted(grand.items()): print(f"   {k}: {v}")

if __name__ == "__main__":
    main(sys.argv[1:])
