"""u-helpers-ours: frame-chain census of the U8 local word that #942 uses (unchanged U8 word, 333 calls, T37 owner packing).
Adapted from bound-check/g-x4-order/census.py (#799 X4). Prices with the m-anatomy pricer (exact re-implementation of the hub
engine) on #942's global decomposition: 333 x local {1:2548,2:338,3:446,4:162,5:72,8:8} + idle {48:11086, 1:6}, m 73, W 33152.
A per-call local delta dl moves -333*sum(r*dl[r]) units to/from the idle remainder at rank 48 (leftover as one child).
Usage: python3 -I census942.py u8-FULL-WORD.json.gz census942.json"""
import sys, json, time, os, math
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import engine as E
from pricer import price_b
from fractions import Fraction as F
LOC = {1: 2548, 2: 338, 3: 446, 4: 162, 5: 72, 8: 8}; IDLE = {48: 11086, 1: 6}; CALLS = 333; M = 73; W = 33152; D = 11840
B942 = F(1664727258740678855772, 10**24)
def global_profile(dl):
    G = {}
    for r in set(LOC) | set(dl): G[r] = CALLS * (LOC.get(r, 0) + dl.get(r, 0))
    dw = sum(r * n for r, n in dl.items()); moved = -CALLS * dw
    for r, n in IDLE.items(): G[r] = G.get(r, 0) + n
    q, rem = divmod(moved, 48); G[48] += q
    if rem: G[rem] = G.get(rem, 0) + 1
    assert all(n >= 0 for n in G.values()), G; return {r: n for r, n in G.items() if n}
def price_delta(dl, certify=False): return price_b(global_profile(dl), M, W, certify=certify)
def pct(b): return float(b / B942 - 1) * 100
if __name__ == '__main__':
    t0 = time.time()
    w = E.load(sys.argv[1]); caps = E.cone_caps(w); res = E.replay(w, caps=caps, record=True)
    assert res['ok'], res['why']
    stored = {int(k): v for k, v in w['report']['local_profile'].items()}
    print('profile', dict(sorted(res['paid'].items())), 'equal stored', dict(res['paid']) == stored, 'equal LOC', dict(res['paid']) == LOC)
    ch = res['chains']; R = w['report']['R']; n = 2 * E.V + R; ctr = set(w['centers'])
    print('V', E.V, 'R', R, 'centres', sorted(ctr))
    def cls(q):
        if q in ctr: return 'centres'
        if q < E.V: return 'data x'
        if q < 2 * E.V: return 'data y'
        h = q - 2 * E.V
        return 'factory' if h < 224 else ('face' if h < 392 else 'M')
    runs = {q: [x[3] for x in ch[q][1:]] for q in range(n)}
    tot = Counter(x for q in range(n) for x in runs[q]); assert dict(tot) == dict(res['paid'])
    bf = float(B942); phif = lambda r: r * ((M / r) ** bf - 1)
    out = {}
    for k in ('data x', 'centres', 'data y', 'factory', 'face', 'M'):
        regs = [q for q in range(n) if cls(q) == k]
        h = Counter(x for q in regs for x in runs[q]); nr = sum(h.values()); mass = sum(r * z for r, z in h.items())
        nn = sum(1 for q in regs if runs[q]); single = Counter(sum(runs[q]) for q in regs if runs[q])
        dl = Counter(single); dl.subtract(h); dl = {r: z for r, z in dl.items() if z}
        b2 = price_delta(dl) if k != 'centres' else B942
        widths = Counter(sum(runs[q]) for q in regs)
        out[k] = dict(n=len(regs), climbing=nn, runs=nr, mean_run=round(mass / nr, 3) if nr else 0,
                      stops_per_register=round((nr - nn) / nn, 3) if nn else 0,
                      F_share_pct=round(100 * CALLS * sum(z * phif(r) for r, z in h.items()) / D, 2),
                      all_single_run_ceiling_pct_b=round(pct(b2), 3), hist=dict(sorted(h.items())), widths=dict(sorted(widths.items())))
        print('%-8s n %3d runs %4d mean run %.2f stops/reg %.2f F-share %.1f%% single-run ceiling %+.2f%% widths %s' %
              (k, len(regs), nr, mass / max(nr, 1), (nr - nn) / max(nn, 1), out[k]['F_share_pct'], out[k]['all_single_run_ceiling_pct_b'], dict(sorted(widths.items()))))
    idle_share = 100 * sum(z * phif(r) for r, z in IDLE.items()) / D
    print('idle F-share %.2f%%' % idle_share)
    dump = {q: dict(cls=cls(q), runs=runs[q], frames=[len(Fr) for s, i, Fr, d in ch[q]], stages=[s for s, i, Fr, d in ch[q]]) for q in range(n)}
    json.dump(dict(classes=out, idle_F_share_pct=idle_share, profile=dict(sorted(res['paid'].items())), registers=dump), open(sys.argv[2], 'w'))
    print('seconds', round(time.time() - t0, 1))
