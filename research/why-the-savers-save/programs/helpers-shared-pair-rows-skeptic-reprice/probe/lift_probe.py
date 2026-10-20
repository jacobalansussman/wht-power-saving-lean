#!/usr/bin/env python3
"""v-pair-market-reprice lift probe (frame-only).  The pair-market design serves each target's own-block need FIRST:
first read of a pair row at own line L1, then the raised read of the shared second-half row, so the target sits at
its own plane P_own(t) (rank 3) before any face read.  U8 instead reads the second half last (final reads at U_t).
This replays U8's B and tail phases with join semantics (gate frame := F[f] + cur[a] + cur[b], both registers climb),
once unchanged (sanity: must equal U8) and once with each target lifted to P_own(t) at its rank-2 stop, and reports
targets whose frames leave u_t^perp, low-end helpers that leave their final frame, and helper frame growth.
Usage: python3 -I lift_probe.py WORD DEMAND BLOCKS"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
d = json.load(open(sys.argv[2])); U8 = json.load(open(sys.argv[3]))['blocks']
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']; h = c['h']
blk = {s: i for i, x in enumerate(U8) for s in x}
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    # reduced canonical form
    bb = sorted(bb, reverse=True)
    out = []
    for y in bb:
        for i in range(len(out)):
            if out[i] ^ y < out[i]: pass
        out.append(y)
    return frozenset(span_canon(out))
def span_canon(gens):
    piv = {}
    for a in gens:
        while a:
            bit = a.bit_length() - 1
            if bit in piv: a ^= piv[bit]
            else: piv[bit] = a; break
    for k in sorted(piv):
        for j in piv:
            if j > k and piv[j] >> k & 1: piv[j] ^= piv[k]
    return tuple(sorted(piv.values(), reverse=True))
def B_(gens): return span_canon(list(gens))
def cont(a, b): return B_(list(a) + list(b)) == B_(b)
perp = lambda p: B_([x for x in range(1, 1 << h) if bin(x & p).count('1') % 2 == 0])
dem = {int(t): {int(s): x for s, x in dd.items()} for t, dd in d['demand'].items()}
# own plane per target: the 4-subset (within one U8 block) that has more than one sign pattern among its readers
sub_readers = defaultdict(list)
for t in range(v):
    per = defaultdict(list)
    for s in dem[t]: per[blk[s]].append(s)
    for b, ss in per.items(): sub_readers[tuple(sorted(ss))].append(t)
own = {}
for ss, ts in sub_readers.items():
    if len(ts) == 4:
        for t in ts: own[t] = B_([ports[s] for s in ss])
print('targets with an own plane:', len(own), ' own-plane dims:', Counter(len(x) for x in own.values()))
print('own plane inside u_t^perp:', sum(cont(own[t], perp(ports[t])) for t in own), 'of', len(own))
endf = [B_(F[f]) for f in c['final']]
def replay(aug):
    cur = [B_(F[f]) for f in c['start']]
    for ph in ['A', 'precondition', 'cut', 'compensation']:
        for e in c.get(ph) or []:
            f = e[-1]
            if f is None: continue
            for r in (e[1:3] if e[0] == 'add' else e[1:2]):
                if cont(cur[r], F[f]): cur[r] = B_(F[f])
    for r, f in zip(c['centers'], c['center_frames']): cur[r] = B_(F[f])
    rank2 = set(); viol = Counter(); grown = Counter(); orig_mismatch = 0
    for ph in ['B', 'tail']:
        for e in c.get(ph) or []:
            f = e[-1]
            if f is None: continue
            regs = e[1:3] if e[0] == 'add' else e[1:2]
            G = B_(list(F[f]) + [x for r in regs for x in cur[r]])
            if G != B_(F[f]): orig_mismatch += 1
            for r in regs:
                if aug and v <= r < 2*v and len(G) >= 2 and r not in rank2:
                    rank2.add(r); G2 = B_(list(G) + list(own[r - v])); 
                    cur[r] = G2
                else: cur[r] = G
            if aug:
                for r in regs:
                    if v <= r < 2*v and r in rank2: cur[r] = B_(list(cur[r]) + list(own[r - v]))
    bad_t = [t for t in range(v) if not cont(cur[v + t], endf[v + t])]
    full = B_([1 << i for i in range(h)])
    low = [r for r in range(2*v, n) if endf[r] != full]
    bad_h = [r for r in low if not cont(cur[r], endf[r])]
    return cur, orig_mismatch, bad_t, bad_h
cur0, mm0, bt0, bh0 = replay(False)
print('unaugmented: gates whose frame had to grow:', mm0, ' targets outside u^perp:', len(bt0), ' low-end helpers outside final:', len(bh0))
cur1, mm1, bt1, bh1 = replay(True)
print('augmented (targets at P_own from their rank-2 stop): gates whose frame grows:', mm1,
      ' targets outside u^perp:', len(bt1), ' low-end helpers outside final:', len(bh1))
hg = Counter(len(cur1[r]) - len(cur0[r]) for r in range(2*v, n))
print('helper end-of-B frame growth (dims added -> count):', dict(sorted(hg.items())))
json.dump(dict(mm0=mm0, mm1=mm1, bad_targets=bt1, bad_lowend_helpers=bh1, helper_growth=dict(hg)),
          open(sys.argv[4] if len(sys.argv) > 4 else '/dev/null', 'w'))
