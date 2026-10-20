#!/usr/bin/env python3
"""d-bit-exits / kernel_relax.py: how far is the response-kernel entrance from feasible if frames are ignored?
For every helper p: is K_p in the F2 span of R_d(k) over ALL z registers d != p and ALL times k <= first write of p
(frames ignored)? And with donors restricted to registers whose frame at k is nested with E (inside E or containing E)?
usage: python3 -I kernel_relax.py WORD.json.gz"""
import gzip, json, os, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames_lib import Frames
w = json.load(gzip.open(sys.argv[1], 'rt'))
V, R = w['v'], w['R']; N = 2*V+R; F = Frames(w)
A, B = w['A'], w['B']; helpers = range(2*V, N)
isz = lambda r: not (V <= r < 2*V)
first_write = {}; first_touch = {}; hist = defaultdict(list); cur = list(w['start'])
for k, (f, t, ss) in enumerate(A):
    for r in [t, *ss]:
        if cur[r] != f: hist[r].append((k, f)); cur[r] = f
        first_touch.setdefault(r, k)
    first_write.setdefault(t, k)
def frame_at(r, k):
    fr = w['start'][r]
    for kk, f in hist[r]:
        if kk < k: fr = f
        else: break
    return fr
resp = [0]*N
for t in range(V): resp[V+t] = 1 << t
for f, t, ss in reversed(B):
    for s in ss: resp[s] ^= resp[t]
for t, ks in enumerate(w['scatter']):
    for k in ks: resp[w['ret'][k]] ^= 1 << t
respA = [None]*(len(A)+1); respA[len(A)] = list(resp)
for k in range(len(A)-1, -1, -1):
    f, t, ss = A[k]
    for s in ss: resp[s] ^= resp[t]
    respA[k] = list(resp)
K = {h: respA[0][h] for h in helpers}
def span_has(vs, x):
    b = {}
    for v in vs:
        while v:
            p = v.bit_length()-1
            if p in b: v ^= b[p]
            else: b[p] = v; break
    while x:
        p = x.bit_length()-1
        if p not in b: return False
        x ^= b[p]
    return True
st = Counter()
for p in helpers:
    E = A[first_touch[p]][0]; fw = first_write.get(p, len(A))
    allv = set(); nested = set()
    for k in range(0, fw+1):
        for d in range(N):
            if d == p or not isz(d): continue
            allv.add(respA[k][d])
            fd = frame_at(d, k)
            if F.le(fd, E) or F.le(E, fd): nested.add(respA[k][d])
    st['any donor, any frame: %s' % span_has(allv, K[p])] += 1
    st['donor frame nested with E: %s' % span_has(nested, K[p])] += 1
print(json.dumps(dict(st), indent=1))
