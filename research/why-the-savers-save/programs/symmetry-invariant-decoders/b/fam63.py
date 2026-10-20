"""fam63.py (t-km): rebuild the C7-invariant cost-63 family that passes the octad Hall test (km_p7_d7_c63.json),
verify Hall on ALL 1260 octads with d-R7's hall_deficit, report its structure, then masked ALS (d-R7 als.py) with
restarts: is there a real fitting decoder B = X R (read masks = ports of the frames)?  usage: python -I fam63.py RESTARTS ITERS"""
import sys, json
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
from als import als
import numpy as np
z = np.load('gens.npz'); perm = z['perm7'].astype(np.int64); inv = np.argsort(perm)
km = json.load(open('km_p7_d7_c63.json'))
fixed = [(sorted(f), dd) for f, dd in json.load(open('fixed_flats_g7.json'))]
fixed = [(f, dd) for f, dd in fixed if dd < 9]
fams = []
for r in km['results']:
    for rep in r['passing_reps']:
        fams.append((rep, r['fixed_idx']))
print('passing families:', len(fams))
out = []
for rep, fidx in fams:
    S = []
    m = np.zeros(NP, bool); m[rep] = True
    for k in range(7):
        S.append(m.copy()); m = m[inv]
    for fi in fidx:
        mm = np.zeros(NP, bool); mm[fixed[fi][0]] = True; S.append(mm)
    S = np.array(S)
    dims = [f2dim([LAB[i] for i in np.nonzero(s)[0]]) for s in S]
    hd = hall_deficit(S)
    print('dims', dims, 'cost', sum(dims), 'ports per frame', S.sum(1).tolist(), 'fixed idx', fidx,
          'Hall deficit sum over all 1260 octads', int(hd.clip(0).sum()), 'octads with slack 0:', int((hd == 0).sum()))
    # every port read by how many frames
    print('  frames per port:', sorted(__import__('collections').Counter(S.sum(0).tolist()).items()))
    rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    best = 1e9
    for rs in range(int(sys.argv[1])):
        res, X, R, h = als(S, ORTH, iters=int(sys.argv[2]), rng=rng)
        best = min(best, res)
        print('  restart', rs, 'residual %.3e' % res, 'hist', ['%.2e' % x for x in h[::max(1, len(h) // 6)]], flush=True)
    out.append(dict(rep=rep, fixed_idx=fidx, dims=dims, frames=[np.nonzero(s)[0].tolist() for s in S], best_residual=best))
json.dump(out, open('fam63.json', 'w'), indent=1)
