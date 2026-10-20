# lever-reprice: price the recursion lever (rl-plan's separated variants) on a word given its paid profile, the counted
# literal head blocks (rank 48) and the other idle rank-48 blocks. Same mixed-target root as rl-plan/tools/rescue.py.
# Controls first: T74 (#1839) must reproduce rl-plan's 1.6948358e-3 (full cut), 1.6968621e-3 (half cut), break-even 9,443.
# Run: python3 -I -B tools/lever_price.py
import sys, json
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mixed import root, add, sep, D

def variants(name, prof, W, heads, non, net=None, c_bit=None, quiet=False):
    """prof: full paid profile (including the idle rank-48 blocks); heads/non: literal heads / other rank-48 idle.
    net: optional {rank: count} of HEAD-net(argued) remainders (cut fully in the optional rows)."""
    base = dict(prof); base[48] -= heads + non; assert base[48] >= 0
    if not base[48]: del base[48]
    I = heads + non
    out = {}
    def row(key, ch, tg, label):
        b = root(ch, tg); s = sep(ch, tg)
        out[key] = dict(b=b, separated=s, label=label)
        return b
    hub = row('uniform', add(base, 48, I), {73: W}, 'uniform (the hub word)')
    row('full_cut', add(base, 24, 2 * non), {73: W - heads, 25: heads}, 'heads cut to target 25, other 48s split 24+24')
    row('half_cut', add(add(base, 24, heads), 48, non), {73: W - heads, 49: heads}, 'heads cut by 24 axes (target 49), rest paid')
    row('unsep_cut', add(base, 48, non), {73: W - heads, 25: heads}, 'heads cut to 25, others kept at 48 (NOT separated)')
    row('all_cut', dict(base), {73: W - I, 25: I}, 'all idle 48s cut (structure false)')
    # best separated d (rescue.py's sweep)
    best = None
    for d in range(1, 49):
        tgt = 73 - d; ch = dict(base)
        if 48 - d: ch = add(ch, 48 - d, heads)
        cap = tgt - 1
        if 48 <= cap: ch = add(ch, 48, non)
        else:
            rest = 48; pieces = []
            while rest: p = min(cap, rest); pieces.append(p); rest -= p
            for p in pieces: ch = add(ch, p, non)
        tg = {73: W - heads, tgt: heads}
        if not sep(ch, tg): continue
        b = root(ch, tg)
        if best is None or b > best[1]: best = (d, b)
    out['best_d'] = dict(d=best[0], b=best[1])
    if net:
        # optional [argued]: also cut the net heads fully (each remainder r: target 73-r, block removed), on top of half cut
        ch = add(add(base, 24, heads), 48, non); tg = {73: W - heads - sum(net.values()), 49: heads}
        for r, n in net.items():
            ch[r] -= n
            if not ch[r]: del ch[r]
            tg[73 - r] = tg.get(73 - r, 0) + n
        row('half_cut_plus_net', ch, tg, 'half cut + net heads (argued) cut fully')
    # break-even heads for the full cut (others split 24+24) at fixed total idle I
    def full(t): return root(add(base, 24, 2 * (I - t)), {73: W - t, 25: t})
    lo, hi = 0, I
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if full(mid) > hub: hi = mid
        else: lo = mid
    out['breakeven_full'] = hi
    if not quiet:
        print(f'== {name}: W={W}, idle rank-48 = {I}, literal heads = {heads} ({heads / I * 100:.2f} %), others = {non}')
        for k in ('uniform', 'full_cut', 'half_cut', 'unsep_cut', 'all_cut', 'half_cut_plus_net'):
            if k not in out: continue
            r = out[k]; m = r['b'] - hub
            print(f"  {r['label']:52s} b={r['b']:.10e}  {float(m / hub) * 100:+.3f} % ({float(m):+.4e})  separated={r['separated']}")
        bd = out['best_d']
        print(f"  best separated d = {bd['d']}: b={bd['b']:.10e} ({float((bd['b'] - hub) / hub) * 100:+.3f} %)")
        print(f"  full-cut break-even: {out['breakeven_full']} heads ({out['breakeven_full'] / I * 100:.2f} % of {I}); have {heads}")
        if c_bit is not None:
            for k in ('full_cut', 'half_cut'):
                b = out[k]['b']
                print(f"  kappa with bit c={c_bit}: {k}: min(b, c) = {min(b, D(c_bit)):.10e} (uniform: {min(hub, D(c_bit)):.10e})")
    return out

if __name__ == '__main__':
    print('##### controls: T74 (#1839), rl-plan values 1.6948358e-3 / 1.6968621e-3 / break-even 9,443')
    T74 = {1: 1560432, 2: 336336, 3: 278724, 4: 84912, 5: 54612, 8: 5328, 48: 22748}
    o = variants('T74 #1839', T74, 66304, 11010, 11738)
    print('##### #3945 T360 staged544 terminal-gap fill (heads counted by tools/heads_fill.py)')
    P3945 = {1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920, 48: 137952}
    o = variants('#3945 T360', P3945, 322560, 57312, 80640, net={3: 23328, 5: 4320}, c_bit='0.001733427039589844234749')

# ---- other words: same parent template (T37 envelope, sha 21297b20; pinned in #3945 owner/SOURCES.json and in
# #3148's SOURCE-PINS.json), so every FULLY idle rank-48 bin of children 8/10 (pair 3) is a literal head and every idle
# bin of 12/13 (pair 0) is not [argued from the schedule; confirmed with 0 exceptions on #3945 by heads_fill.py].
# Only the idle profile is public in the bodies, so heads lie in [idle48 - 224 T, min(idle48, 224 T)]; on #3945 the
# measured count equals the lower end (pair 0 never hosted). Half-cut prices increase with heads (break-even 1 head).
def bounds_rows():
    WORDS = [
        ('#3148 T740 (merged, archive)', {1: 16376940, 2: 4821840, 3: 1567431, 4: 679320, 5: 147075, 8: 53280, 48: 282569}, 740),
        ('#3605 T740 (open; histogram #3659)', {1: 16417896, 2: 4835160, 3: 1540174, 4: 679321, 5: 143992, 6: 3, 8: 53280, 48: 283185}, 740),
        ('#3642 T530 (open)', {1: 11747920, 2: 3453480, 3: 1112042, 4: 486540, 5: 103690, 8: 38160, 48: 202828}, 530),
        ('#3653 T488 (open)', {1: 10817496, 2: 3179808, 3: 1024434, 4: 447984, 5: 94794, 8: 35136, 48: 186782}, 488),
        ('#3682 T360 (open)', {1: 7980120, 2: 2345760, 3: 755730, 4: 330480, 5: 69930, 8: 25920, 48: 137790}, 360),
        ('#3815 T360 (open; same idle-48 bins as #3945)', {1: 7989840, 2: 2349000, 3: 749088, 4: 330480, 5: 69120, 8: 25920, 48: 137952}, 360),
        ('#3945 T360 (open)', {1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920, 48: 137952}, 360),
    ]
    print('##### lever on each word, heads at the lower / upper end of the public bound')
    print(f"{'word':48s} {'b (hub word)':>15s} {'heads lo..hi':>17s} {'half cut lo':>15s} {'half cut hi':>15s} {'full cut lo':>15s} {'BE full':>8s}")
    for name, prof, T in WORDS:
        W = 896 * T; I = prof[48]; lo = max(0, I - 224 * T); hi = min(I, 224 * T)
        a = variants(name, prof, W, lo, I - lo, quiet=True)
        z = variants(name, prof, W, hi, I - hi, quiet=True)
        print(f"{name:48s} {a['uniform']['b']:.10e} {lo:7d}..{hi:7d} {a['half_cut']['b']:.10e} {z['half_cut']['b']:.10e} {a['full_cut']['b']:.10e} {a['breakeven_full']:8d}")
bounds_rows()
