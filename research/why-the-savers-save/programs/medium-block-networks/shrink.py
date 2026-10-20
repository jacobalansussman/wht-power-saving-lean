#!/usr/bin/env python3
"""cx-medium / shrink.py (own code, stdlib; python3 -I).  What shrinking the hub's template costs, with the hub's U8
local word (our lowering E1, fsp-price/data/words.json) and a GENEROUS hosting model.
Template with P bank pairs: 2P+1 stages GF,GB,GF,... (pair 1 in stages 0,1,2; pair q >= 2 in stages 2q-1, 2q, joined by
the free twin gates of Gates.lean), windows of sigma = 8 fresh coordinates + the shared line, m = 1 + 8(2P+1).
Per replica: 2Pv data arrays (v = 112), each paying m - 1: 8 per active stage inside the call, the rest as idle
children (the hub's CH table generalised: pair q >= 2 has A_q = 8(2q-2) [stages 0..2q-3] and B_q = 8 [stage 2q-2];
pair p < P has F_p = 8(2P-2p) [stages 2p+1..2P]); 2P+1 calls of the U8 word (profile and 8 copies of rank 8 each);
504 helper leftover pieces per call ({4:80, 5:32, 7:54, 8:122, 9:216}, 3778 levels).
Generous hosting ("exact fill", an-ours knobs.py generalised): a piece of call k can go to any idle child whose window
contains stage k; every piece that fits by capacity is hosted (computed by a max-flow on aggregated capacities per
(child type, call)); the idle remainder is cut into pieces of the largest child rank of the layout, 8(2P-2) (48 at
P = 4); unhosted pieces go to banks with IDEAL packing (no pads; banks are D-neutral).  At P = 4 this must reproduce
an-ours' exact fill 1.6647327e-3.  These are UPPER bounds on what the U8 word can give in each shrunk layout."""
import sys, os, json
from collections import Counter
from fractions import Fraction
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../anatomy-our-lean-words')))
from anatomy_ours import price, rank, narrow

PROF = {1: 2548, 2: 338, 3: 446, 4: 162, 5: 72, 8: 8}
WID = {4: 80, 5: 32, 7: 54, 8: 122, 9: 216}
V, CST, SIG = 112, 64, 8
LEFT = sum(w * n for w, n in WID.items())

def layout(P):
    S = 2 * P + 1; m = 1 + SIG * S
    ch = []                                        # (name, rank, set of stages it covers), 2v instances each
    for q in range(2, P + 1):
        ch.append((f'A{q}', SIG * (2 * q - 2), set(range(0, 2 * q - 2))))
        ch.append((f'B{q}', SIG, {2 * q - 2}))
    for p in range(1, P):
        ch.append((f'F{p}', SIG * (2 * P - 2 * p), set(range(2 * p + 1, S))))
    return S, m, ch

def host(P):
    """greedy max hosting on aggregated capacities: per child type, capacity 2v * rank levels and at most 2v pieces per
    call (one piece per instance per call).  Pieces largest first; a piece of call k goes to the covering child type
    with the most remaining capacity.  Returns hosted levels per child type and the unhosted pieces."""
    S, m, ch = layout(P)
    cap = {nm: 2 * V * r for nm, r, _ in ch}
    per_call = {(nm, k): 2 * V for nm, _, st in ch for k in st}
    unhosted = Counter(); hosted = Counter()
    for k in range(S):
        for w in sorted(WID, reverse=True):
            for _ in range(WID[w]):
                best = None
                for nm, r, st in ch:
                    if k in st and cap[nm] >= w and per_call[(nm, k)] > 0 and r >= w:
                        if best is None or cap[nm] > cap[best]: best = nm
                if best is None: unhosted[w] += 1
                else:
                    cap[best] -= w; per_call[(best, k)] -= 1; hosted[best] += w
    return S, m, ch, hosted, unhosted

def evaluate(P, T=None):
    S, m, ch, hosted, unhosted = host(P)
    data = 2 * P * V
    idle_total = sum(2 * V * r for _, r, _ in ch)
    rem = idle_total - sum(hosted.values())
    unh = sum(w * n for w, n in unhosted.items())
    # replicas T so that banks are integral under ideal packing
    T = T or next(t for t in range(8, 4000, 8) if (t * unh) % m == 0)
    H = Counter({r: S * T * n for r, n in PROF.items()})
    rmax = max([r for _, r, _ in ch], default=0)
    if rem:
        q, rr = divmod(T * rem, rmax)
        H[rmax] += q
        if rr: H[rr] += 1
    banks = T * unh // m
    W = T * data + banks
    D = m * W - rank(dict(H))
    assert D == T * (data - S * CST), (D, T * (data - S * CST))
    a = float(price(dict(H), W, m=m))
    return dict(P=P, stages=S, m=m, T=T, W=W, banks=banks, D=D, D_per_replica=data - S * CST, idle_per_call=idle_total / S,
                hosted_per_call=sum(hosted.values()) / S, unhosted_per_call=unh / S, a=a, N=narrow(dict(H), m))

if __name__ == '__main__':
    out = []
    for P in (1, 2, 3, 4, 5):
        r = evaluate(P); out.append(r)
        print(f"P={P} stages={r['stages']} m={r['m']:3d}: D/replica {r['D_per_replica']:4d}; idle cover/call {r['idle_per_call']:7.1f},"
              f" hosted/call {r['hosted_per_call']:7.1f} of {LEFT}, to banks/call {r['unhosted_per_call']:7.1f};"
              f" T={r['T']} W={r['W']} banks={r['banks']}  price <= {r['a']:.7e} (x{r['a']/1.6797550e-3:.3f} of the hub)")
    p4 = out[3]['a']
    print('calibration: P=4 must equal an-ours exact fill 1.6647327e-3:', f"{p4:.7e}", abs(p4 - 1.6647327e-3) < 1e-10)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'shrink.json'), 'w'), indent=1)
