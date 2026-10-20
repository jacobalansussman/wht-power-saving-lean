#!/usr/bin/env python3
"""ctr-search / census.py: label-spanned subspaces W of F_2^9 and the ports they hold.
For each dimension d: number of label-spanned W, max |S_W|, histogram of |S_W|, and for the largest ones the
F2 rank of B[S_W, S_W] (the pairing cost of a term with R = A = S_W) and how many 8-cliques S_W meets / contains."""
import sys, os, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *

def rref_key(vs):
    rows = []
    for v in vs:
        for r in rows:
            if (v >> (r.bit_length() - 1)) & 1: v ^= r
        if v:
            for i, r in enumerate(rows):
                if (r >> (v.bit_length() - 1)) & 1: rows[i] = r ^ v
            rows.append(v)
    return tuple(sorted(rows))

def elements(key):
    els = [0]
    for r in key: els += [x ^ r for x in els]
    return els

lab2port = {u: i for i, u in enumerate(LABS)}

def ports_in(key):
    return sorted(lab2port[x] for x in elements(key) if x in lab2port)

if __name__ == '__main__':
    om, cl = max_cliques()
    clmask = [0] * N
    for ci, c in enumerate(cl):
        for p in c: clmask[p] |= 1 << ci
    layer = {rref_key([u]) for u in LABS}
    out = {}
    for d in range(1, 10):
        hist = collections.Counter()
        best = []
        for key in layer:
            S = ports_in(key)
            hist[len(S)] += 1
            best.append((len(S), key))
        best.sort(reverse=True)
        top = []
        seen_sizes = set()
        for sz, key in best[:2000]:
            S = ports_in(key)
            m = 0
            for p in S: m |= clmask[p]
            pr = pair_cost(S, S)
            top.append((sz, pr, bin(m).count('1')))
        tops = collections.Counter(top)
        out[d] = dict(count=len(layer), hist=sorted(hist.items()), top=sorted(tops.items(), reverse=True)[:8])
        print(f"d={d}: {len(layer)} label-spanned W; |S_W| histogram (top) {sorted(hist.items())[-6:]}; "
              f"(|S|, pair rank B[S,S], cliques met) among largest: {sorted(tops.items(), reverse=True)[:5]}", flush=True)
        if d == 9: break
        nxt = set()
        if len(layer) > 400000:
            print('  stopping growth (too many)'); break
        for key in layer:
            els = set(elements(key))
            for u in LABS:
                if u not in els: nxt.add(rref_key(list(key) + [u]))
        layer = nxt
    json.dump({str(k): v for k, v in out.items()}, open(os.path.join(HERE, 'census.json'), 'w'), indent=1)
