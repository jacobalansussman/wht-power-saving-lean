#!/usr/bin/env python3
"""fp-measure: distinct POSITIONS (vertices of the frame graph of the window) among the raw labels on each paid move,
for one F1 call of each local word (fm_ledger.run with pos=True).  usage: python3 -I pos_run.py NAME=word.gcert.json.gz [NAME=...]"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fm_ledger as F
words = [tuple(a.split("=", 1)) for a in sys.argv[1:]]   # NAME=path/to/word.gcert.json.gz ...
out = {}
for n, p in words:
    r = F.run(F.load(p), pos=True)
    out[n] = dict(M=r["M"], P=r["P"], P_pos=r["P_pos"], pos_hist=r["pos_hist"], maxlab=r["maxlab"])
    print("%s: M=%d  label packing P=%d (%.2f per move)  position packing P_pos=%d (%.2f per move)  positions per move %s"
          % (n, r["M"], r["P"], r["P"] / r["M"], r["P_pos"], r["P_pos"] / r["M"], r["pos_hist"]))
json.dump(out, open("out/pos_run.json", "w"), indent=1)
