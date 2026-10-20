"""vbpm.py (r3-classes): INDEPENDENT check of the new exclusions.  Zero cells of E are re-derived by the skeptic's own
code (vrec_z.py = v-rm2-exact vrec.py unchanged except that it keeps its zero set; own field, self-dual basis
[3,5,12,17,26], flats, recursion, pooling); the block-plus-matching test is then run by bpm.c (fresh C code):
  F even: rank E[D_F,D_F] = |D_F|;  G any form, S = D_F n D_G: rank E[S,S] <= s_G;  E[D_F,D_F] = (E[S,S] padded) + rest,
  so |D_F| <= s_G + termrank(D_F x D_F minus S x S, avoiding zeros).  A violation excludes the target set.
usage: vbpm.py [--maxd 16] [--rounds 3] < lines "k1,...,kW"."""
import sys, time, subprocess, os
import numpy as np
import vrec_z as vr
HERE = os.path.dirname(os.path.abspath(__file__))
a = sys.argv[1:]
maxd = int(a[a.index("--maxd") + 1]) if "--maxd" in a else 16
rounds = int(a[a.index("--rounds") + 1]) if "--rounds" in a else 3
for line in sys.stdin:
    line = line.strip()
    if not line: continue
    ks = [int(x) for x in line.split(",")]
    t0 = time.time()
    v = vr.run(ks, maxd=maxd, rounds=rounds); Z = vr.LASTZ
    w = vr.CLSPAR[:, ks].sum(axis=1); r = w[vr.VID]; s = vr.WT - r
    Zm = [0] * 32
    for (x, y) in Z: Zm[x] |= 1 << (y - 1)
    txt = "%d\n" % (vr.NF - 1) + "".join("%d %d %d\n" % (int(vr.MASK[q]), int(vr.WT[q]), int(s[q])) for q in range(1, vr.NF))
    txt += " ".join("%x" % Zm[x] for x in range(1, 32)) + "\n"
    out = subprocess.run([os.path.join(HERE, "bpm"), "20"], input=txt, capture_output=True, text=True).stdout.strip()
    print(line, "| vrec", v, "| zeros", len(Z), "| BPM", out, "| %.0fs" % (time.time() - t0), flush=True)
