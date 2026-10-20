"""For each line system: v, Q-rank n, ip distribution, and a large pairwise NON-orthogonal set (lower bound for the
F_2 minrank of the orthogonality pattern = lower bound for the complex-side window h when the decoder is the Gram)."""
import sys, os, time, json, random, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import linesys

def clique(L, seeds=(1, 2, 3), restarts=20, cap=None):
    v = len(L); best = []
    for sd in seeds:
        rng = random.Random(sd)
        for _ in range(restarts):
            s = rng.randrange(v)
            C = [s]; cand = np.nonzero(L @ L[s])[0]; cand = cand[cand != s]
            while len(cand):
                if len(cand) > 3000: pick = cand[rng.randrange(len(cand))]
                else:
                    sub = L[cand]; M = (sub @ sub.T) != 0
                    deg = M.sum(1); top = np.nonzero(deg == deg.max())[0]
                    pick = cand[top[rng.randrange(len(top))]]
                C.append(int(pick)); ip = L[cand] @ L[pick]
                cand = cand[(ip != 0) & (cand != pick)]
            if len(C) > len(best): best = C
    Lb = L[best]; G = Lb @ Lb.T
    assert np.all(G != 0), 'not a clique'
    return sorted(best)

if __name__ == '__main__':
    out = {}
    for name in sys.argv[1:]:
        t = time.time(); L = getattr(linesys, name)(); v = len(L)
        n = int(np.linalg.matrix_rank(L.astype(float)))
        norms = set((L * L).sum(1).tolist())
        rec = dict(v=v, n=n, norms=sorted(norms))
        if v <= 5000:
            G = L @ L.T; nrm = G[0, 0]
            vals, cnt = np.unique(np.abs(G[np.triu_indices(v, 1)]) * 4 // nrm, return_counts=True)
            rec['abs_ip_norm4_counts'] = {int(a): int(b) for a, b in zip(vals, cnt)}
        C = clique(L, restarts=(20 if v < 5000 else 4))
        rec['clique'] = len(C); rec['seconds'] = round(time.time() - t, 1)
        out[name] = rec; print(name, rec, flush=True)
        np.save('clique_%s.npy' % name, np.array(C))
    json.dump(out, open('survey_%s.json' % '_'.join(sys.argv[1:]), 'w'), indent=1)
