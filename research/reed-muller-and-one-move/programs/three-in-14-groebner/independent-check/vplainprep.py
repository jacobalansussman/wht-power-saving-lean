"""vplainprep.py: must-pass control data. For each of my 44 sets: the plain 3 in 15 (deleted move restored),
U+ = every copy of the certificate's directions and of the deleted direction, and every form support inside U+."""
import json
from vcore import Inst, dir_vec, tgt_vec, EXP
import vscan_frames as VF
sets = json.load(open("vsets.json")); certs = {c["n"]: c for c in json.load(open("vcerts_local.json"))}
out = {}
for n, s in enumerate(sets):
    fr = [VF.frames[k][f] for k, f in zip(s["ks"], s["frames"])]
    moves = [x for F in fr for x in F]
    I = Inst([tgt_vec(EXP[(2 * k) % 31]) for k in s["ks"]], [dir_vec(x) for x in moves])
    Ud = set(certs[n]["U_dirs"]) | {s["delete"]}
    U = [i for i, x in enumerate(moves) if x in Ud]; Um = sum(1 << i for i in U); sups = {}
    for q in range(1, 1 << 15):
        dm = int(I.Dmask[q])
        if dm and dm & ~Um == 0: sups[dm] = [int(I.Pmask[q]) >> l & 1 for l in range(3)]
    out[n] = {"moves": moves, "U": U, "sups": [[[i for i in U if dm >> i & 1], p] for dm, p in sorted(sups.items())]}
json.dump(out, open("vplain_sups.json", "w"))
from collections import Counter
print("U+ sizes", Counter(len(v["U"]) for v in out.values()), "supports", Counter(len(v["sups"]) for v in out.values()))
