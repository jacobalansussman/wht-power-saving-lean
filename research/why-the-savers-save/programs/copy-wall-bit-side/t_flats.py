"""t_flats.py (d-bitspan): enumerate all flats of rank <= 7 exactly (rank by rank, lean memory), and for each flat
the number of bit cliques it meets on P120 and on the two P111 punctures.  Saves per rank: the met distribution, the
best flats, and the masks of rank 6/7 flats missing few cliques (families for the SAT).  usage: venv -I t_flats.py"""
import sys, time, collections, pickle; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import numpy as np
from e8geo import *
from e8geo import _normrows
UNITS = [('P111-s2', PUNCT_S2), ('P111-word', PUNCT_WORD), ('P120', [])]
CL = {}
for name, P in UNITS:
    live = [i for i in range(NL) if i not in P]
    C = cliques(live, 8)
    inc = np.zeros((NL, len(C)), dtype=np.float32)
    for j, c in enumerate(C):
        for T in c: inc[T, j] = 1
    CL[name] = (C, inc, live)
    print(name, 'cliques', len(C), flush=True)
LM = [1 << i for i in range(NL)]
level = {0: np.eye(8, dtype=np.int8).tobytes()}
SUMMARY = {}; KEEP = {}
for r in range(1, 8):
    t0 = time.time(); new = {}
    for mask, Nb in level.items():
        N = np.frombuffer(Nb, dtype=np.int8).reshape(-1, 8).astype(np.int64)
        C = _normrows(L @ N.T)
        idx = np.nonzero(np.any(C != 0, axis=1))[0]
        keys, inv = np.unique(C[idx], axis=0, return_inverse=True)
        inv = np.asarray(inv).reshape(-1)
        groups = [[] for _ in range(len(keys))]
        for a, gi in zip(idx.tolist(), inv.tolist()): groups[gi].append(a)
        for mem in groups:
            m2 = mask
            for i in mem: m2 |= LM[i]
            if m2 in new: continue
            l = L[mem[0]]; t = N @ l
            nzt = [j for j in range(len(t)) if t[j] != 0]
            piv = min(nzt, key=lambda j: abs(t[j]))
            rows = [t[piv] * N[j] - t[j] * N[piv] for j in range(len(N)) if j != piv]
            N2 = _normrows(np.array(rows, dtype=np.int64).reshape(len(rows), 8)) if rows else np.zeros((0, 8), np.int64)
            assert np.abs(N2).max(initial=0) < 127
            new[m2] = N2.astype(np.int8).tobytes()
    level = new
    masks = list(level.keys())
    print('rank %d: %d flats (%.1f s)' % (r, len(masks), time.time() - t0), flush=True)
    for name, (C, inc, live) in CL.items():
        met = np.zeros(len(masks), dtype=np.int32)
        nlive = np.zeros(len(masks), dtype=np.int32)
        for a in range(0, len(masks), 4000):
            chunk = masks[a:a + 4000]
            Fm = np.zeros((len(chunk), NL), dtype=np.float32)
            for k, m in enumerate(chunk):
                Fm[k, bits(m)] = 1
            met[a:a + 4000] = ((Fm @ inc) > 0.5).sum(1)
            nlive[a:a + 4000] = Fm[:, live].sum(1)
        miss = len(C) - met
        dist = collections.Counter(miss.tolist())
        SUMMARY[(name, r)] = dict(n=len(masks), maxmet=int(met.max()), ncl=len(C), maxlive=int(nlive.max()),
                                  missdist=sorted(dist.items()))
        print('   %-10s most met %5d of %5d (%.1f per dim), most live lines %3d, flats meeting all %d, fewest misses %s'
              % (name, met.max(), len(C), met.max() / r, nlive.max(), dist.get(0, 0), sorted(dist.items())[:5]), flush=True)
        if r >= 5:
            order = np.argsort(miss)[:20000]
            KEEP[(name, r)] = [(masks[k], int(miss[k])) for k in order]
    pickle.dump(dict(SUMMARY=SUMMARY, KEEP=KEEP), open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/flats_met.pkl'), 'wb'))
    if r == 7:
        pickle.dump(level, open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/flats7.pkl'), 'wb'))
print('done')
