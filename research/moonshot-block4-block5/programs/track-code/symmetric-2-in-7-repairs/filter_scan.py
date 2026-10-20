"""filter_scan.py -- mw2-absorb4: scan the free letter of one summand of a small limit skeleton (full exact test each).
Tests the pair mechanism found in the seed-66 limit: a two-term filter b1 + a b2 seen through two paths needs
(1 + a u) invertible at the K pair and a^2 = 1 at the junk pairs.  Here b2 = b1 w for many free letters w (Paulis,
diagonal Z4 phases, address maps, random engine letters, random Wd monomials); for each, the full residual.
Usage: filter_scan.py lim.pkl move summand nletters starts [twist]"""
import sys, pickle, time, itertools
from twist import *

d = pickle.load(open(sys.argv[1], 'rb'))
t_f, k_f, nlet, starts = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
m, W, H = d['m'], d['W'], d['H']
rec = d['recipe']; tw = d.get('twist', {})
base = build_T(m, W, H, rec, d['skseed'], tw)
G, R = base.G, base.R
n = G.n
kind, summ, Dm = rec[t_f]
src0, b1 = summ[0]
rng = random.Random(17)
ws = []
for a in range(n):
    for b in range(n):
        ws.append(('pauli', G.mono([x ^ a for x in range(n)], [2 * (bin(x & b).count('1') & 1) for x in range(n)])))
for _ in range(nlet):
    ws.append(('engine', rand_free_b(G, rng, m)))
for _ in range(nlet):
    sig = list(range(n)); rng.shuffle(sig)
    ws.append(('wd', G.mono(sig, [rng.randrange(4) for _ in range(n)])))
for _ in range(nlet // 2):
    ws.append(('diag', G.mono(list(range(n)), [rng.randrange(4) for _ in range(n)])))
rng.shuffle(ws)
stats = {}
t0 = time.time(); best = (np.inf, None)
for idx, (typ, w) in enumerate(ws):
    if is_scalar(w):
        continue
    b2 = G.mul(b1, w)
    rec2 = [(k, list(sm), dd) for (k, sm, dd) in rec]
    sm = rec2[t_f][1]
    if k_f < len(sm):
        sm[k_f] = (sm[k_f][0], b2)
    else:
        sm.append((src0, b2))
    try:
        sk = build_T(m, W, H, rec2, d['skseed'], tw, R)
    except Exception:
        continue
    if [1 + len(mv['images']) for mv in sk.moves] != [1 + len(mv['images']) for mv in base.moves]:
        continue
    if any(mv.get('absorbed') is False for mv in sk.moves if mv['kind'] == 'T'):
        continue
    if not sk.setup_final(1500):
        stats.setdefault(typ, []).append(None); continue
    prep_fast(sk)
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
    bb = np.inf
    for st in range(starts):
        x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
        x, c = vp_lm(sk, N, free_idx, x0, iters=80)
        bb = min(bb, c)
    o = order_of(G, w)
    stats.setdefault(typ, []).append(round(bb, 4))
    if bb < best[0]:
        best = (bb, (typ, o, idx))
        pickle.dump({"m": m, "W": W, "H": H, "recipe": rec2, "twist": tw, "skseed": d['skseed'], "res": bb}, open(sys.argv[1].replace('.pkl', f'_fs{t_f}{k_f}_best.pkl'), 'wb'))
    if bb < 1e-12:
        print(f"*** FULL HIT: letter {idx} ({typ}, order {o}) ***", flush=True)
        pickle.dump({"m": m, "W": W, "H": H, "recipe": rec2, "twist": tw, "skseed": d['skseed']}, open(f"HIT_filter_{idx}.pkl", "wb"))
    if idx % 50 == 0:
        print(f"{idx}/{len(ws)}: best {best} ({time.time()-t0:.0f}s)", flush=True)
for typ, v in stats.items():
    vv = [x for x in v if x is not None]
    print(typ, "tested", len(v), "capped", sum(1 for x in v if x is None), "floors (value: count)",
          sorted({x: vv.count(x) for x in set(vv)}.items())[:12])
print("best", best)
