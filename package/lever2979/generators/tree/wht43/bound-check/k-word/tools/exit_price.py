"""x-exit-words: price hub complex words through Shea's complex588 route WITH tier-2 exit frames.

Route (Shea, shea256/fourier-transform-below-nlogn lean/complex588, on Kothe's hub #588 circuit):
  * 112 source-copy helpers (width 9; local profile +112 rank-1 copies, +112 rank-8 climbs), target
    preconditioning moved across scatter (no change of profile), 72 calls, dedicated #401/#377 banks of m = 73
    coordinates; idle charge per source {48:32, 32:32, 16:32, 8:48} (GX.bankPrice377); W = 64*112 + #banks;
    deficit D = 73 W - mass = 2560 must hold; Fourier figure = (largest Lean-plannable a at s = 40) - 1e-10.
TIER-2 CHARGING RULE (the exit-frame generalisation of GInvE / bank_frames):
  helper q with entrance frame E_q and exit frame U_q, E_q <= U_q (checked over GF(2)^9), contributes the
  residual r_q = dim U_q - dim E_q (NOT 9 - dim E_q) to its bank in every call; the bank's final pad is ONE block
  of rank 73 - load (load in [1, 73]); per call each bank serves at most one helper (use k injective).
  The hub's local profile already charges every helper's climb E_q -> U_q, so the deficit identity
  72 (local mass - total residual) = 73*7168 - idle - D is unchanged (checked: D = 2560).
Packings priced:
  cyclic : Shea's cyclic bins (route_price.cyclic_schedule), as in the Lean certificate today;
  ceil   : the pad-minimal table: B = ceil(72 S / 73) banks, all full but one of pad P = 73 B - 72 S
           (= S mod 73 since 72 = -1 mod 73); realisable whenever the multiset {72 x w_q} splits into
           B bins of capacity 73 with that one pad (Koenig edge-colouring gives the per-call table).
  paid   : today's engine (no tier 2): every exit paid up to the top frame (extra rank 9 - dim U climb per
           exit helper, residual 9 - dim E).
"""
import sys, os, json, gzip
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from route_price import unit, price, IDLE, widths_from

M, CALLS, V = 73, 72, 112

def read(p):
    raw = open(p, 'rb').read()
    return json.loads(gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw)

def rank(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)

def word_data(c):
    v, h, R = c['v'], c['h'], c['R']; F = c['frames']; st, fi = c['start'], c['final']
    assert (v, h) == (112, 9)
    E = [F[st[r]] for r in range(2*v, 2*v+R)]; U = [F[fi[r]] for r in range(2*v, 2*v+R)]
    for e, u in zip(E, U):
        assert rank(e) == len(e) and rank(u) == len(u)
        assert rank(list(u) + list(e)) == len(u), 'exit frame does not contain entrance frame'
    src = Counter((len(F[st[r]]), len(F[fi[r]])) for r in range(v))
    tgt = Counter((len(F[st[r]]), len(F[fi[r]])) for r in range(v, 2*v))
    assert src == Counter({(1, 9): 112}) and tgt == Counter({(0, 8): 112}), (src, tgt)
    assert sum(len(F[f]) for f in c['center_frames']) == 64
    prof = {int(k): n for k, n in c['local_profile'].items()}
    es = [len(e) for e in E]; us = [len(u) for u in U]
    return prof, es, us

def ceil_unit(profile, widths, copies=True, extra=None):
    prof = Counter({int(k): n for k, n in profile.items()})
    w = list(widths)
    if copies:
        w += [9] * V; prof[1] += V; prof[8] += V
    if extra:
        for r, n in extra.items(): prof[int(r)] += n
    S = sum(w)
    B = -(-CALLS * S // M)
    P = B * M - CALLS * S
    hist = Counter({r: CALLS * n for r, n in prof.items() if n})
    for r, n in IDLE.items(): hist[r] += n * V
    if P: hist[P] += 1
    W = 64 * V + B
    D = M * W - sum(r * n for r, n in hist.items())
    return dict(hist=dict(sorted(hist.items())), W=W, D=D, banks=B, pad=P, width=S)

def price_all(prof, es, us, kernel=True):
    out = {}
    w_exit = [u - e for e, u in zip(es, us)]
    # cyclic (Shea's schedule) with exit frames
    u = unit(prof, sorted(w_exit, reverse=True)); assert u['D'] == 2560
    p = price(u, kernel=kernel)
    out['exit_cyclic'] = dict(W=u['W'], banks=u['banks'], pads=u['pad_hist'], width=u['width'], **p)
    u = ceil_unit(prof, w_exit); assert u['D'] == 2560, u['D']
    p = price(u, kernel=kernel)
    out['exit_ceil'] = dict(W=u['W'], banks=u['banks'], pad=u['pad'], width=u['width'], **p)
    # today's engine: exits paid up to the top frame
    extra = Counter(9 - x for x in us if x < 9)
    w_top = [9 - e for e in es]
    if extra:
        u = unit(prof, sorted(w_top, reverse=True), extra=extra); assert u['D'] == 2560
        p = price(u, kernel=kernel)
        out['paid_cyclic'] = dict(W=u['W'], banks=u['banks'], width=u['width'], **p)
        u = ceil_unit(prof, w_top, extra=extra); assert u['D'] == 2560
        p = price(u, kernel=kernel)
        out['paid_ceil'] = dict(W=u['W'], banks=u['banks'], pad=u['pad'], width=u['width'], **p)
    else:
        out['paid_cyclic'] = out['exit_cyclic']; out['paid_ceil'] = out['exit_ceil']
    return out

if __name__ == '__main__':
    # usage: exit_price.py OUT.json [--model name profile_json widths_json] [word paths or INDEX dir]...
    outp = sys.argv[1]; args = sys.argv[2:]
    res = {}
    i = 0
    items = []
    while i < len(args):
        a = args[i]
        if a == '--model':
            name, prof, wid = args[i+1], json.loads(args[i+2]), json.loads(args[i+3])
            es, us = [], []
            for key, n in wid.items():   # key "e,u"
                e, u_ = map(int, key.split(','))
                es += [e]*n; us += [u_]*n
            items.append((name, {int(k): n for k, n in prof.items()}, es, us)); i += 4; continue
        if os.path.isdir(a):
            idx = json.load(open(os.path.join(a, 'INDEX.json')))
            for k, s in idx.items():
                if s['v'] != 112 or s['h'] != 9 or s['cst'] != 64: continue
                c = read(os.path.join(a, k + '.json.gz'))
                try: prof, es, us = word_data(c)
                except AssertionError as e: print('skip', k, e); continue
                items.append((k + ' ' + s['origins'][0].split('research/', 1)[-1][:100], prof, es, us))
        else:
            c = read(a); prof, es, us = word_data(c); items.append((a, prof, es, us))
        i += 1
    for name, prof, es, us in items:
        r = price_all(prof, es, us)
        r['profile'] = prof; r['pairs'] = sorted([list(k) + [n] for k, n in Counter(zip(es, us)).items()])
        res[name] = r
        print(f"{name[:70]:70s} exit_cyc={r['exit_cyclic']['dft']} exit_ceil={r['exit_ceil']['dft']} (B={r['exit_ceil']['banks']} P={r['exit_ceil']['pad']}) paid_cyc={r['paid_cyclic']['dft']} paid_ceil={r['paid_ceil']['dft']}", flush=True)
    json.dump(res, open(outp, 'w'), indent=1)
