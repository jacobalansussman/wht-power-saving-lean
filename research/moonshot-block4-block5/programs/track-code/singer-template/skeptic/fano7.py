"""fano7.py -- ms-fano: Fano-gadget PREFIXES of k moves (sampled) + an EXACT completion test that allows ONE extra
move per output (its content solved, not sampled): output l passes if K x_l is in
    Span_{Q[B']}(inputs, prefix generators)  [+ E . Span_{Q[B']}(same) for some candidate paid element E].
E candidates: (a) private last moves (bm.Game.E_candidates: terms one move from the kernel), (b) cleanup moves
(stag2.junk_E: E maps a generator term exactly onto a junk term riding with x_l's kernel-frame term).
Moves of the word = k + (number of outputs that needed an extra move).  Each pass is an exact solve over Q of the
formal identity (closure over frames).  Block 3, W = 2: k = 3 must give no 5-move word (14W/5); k = 4 controls.
Block 4, W = 2: k = 5 with both outputs needing one extra move = 7 moves (2 in 7, the target); k = 6 = 8 (control).
Prefix kinds:
  PL   private chains (x1: ceil(k/2) moves, x2: the rest) -- control.
  FS   Fano start: moves 1..g (g = 2 or 3) are FAT raw-input contents: each spreads both inputs over 2-4 of the seven
       frames N<v>, v in <e1,e2,e3>, with signs (frustration); later moves are route sums (<= 3 summands, exact
       alignment onto one element when possible = the 3-route halving shape).
  FN   (block 4) shared rank-1 move u1 = D b_x (x1 + s x2) with a second raw entry of opposite sign later; then
       gadget moves routing u1 into several of the seven N<x,v> at once, plus route sums.
usage: python3 fano7.py m kind k n seed [maxE] [out.json]"""
import sys, time, json, random
from tsearch import *
from stag2 import junk_E

def fat_raw(ft, rng, W, pts):
    G = ft.G; reps = ft.geo_rep_letters
    nt = rng.choice([2, 3, 4])
    vs = rng.sample(pts, min(nt, len(pts)))
    summ = []
    for j, v in enumerate(vs):
        l = j % W if j < W else rng.randrange(W)
        summ.append((rng.choice([1, -1]), G.lmul(rng.choice(reps[v]), G.inp(l))))
    return summ

def prefix(ft, kind, W, k, rng):
    G = ft.G; m = ft.m
    gens = [G.inp(l) for l in range(W)]
    if kind == "PL":
        per = [k // W + (1 if l < k % W else 0) for l in range(W)]
        for l in range(W):
            cur = l
            for d in range(per[l]):
                g = ft.lead_terms(gens[cur])[0][1]
                R = ft.routes(g)
                if not R:
                    return None
                b = rng.choice(list(R.values()))[0]
                gens.append(ft.move([(1, G.lmul(b, gens[cur]))])); cur = len(gens) - 1
        return gens
    pts = list(range(1, 8))
    if kind == "FS":
        g = rng.choice([2, 3]) if k >= 4 else 2
        for t in range(g):
            u = ft.move(fat_raw(ft, rng, W, pts))
            if u is None:
                return None
            gens.append(u)
    elif kind == "FN":
        bx = ft.x_rep; s = rng.choice([1, -1])
        gens.append(ft.move([(1, G.lmul(bx, G.inp(0))), (s, G.lmul(bx, G.inp(1)))]))
        u1 = gens[-1]
        lt = ft.lead_terms(u1)[0][1]
        R = list(ft.routes(lt).values()); rng.shuffle(R)
        # gadget move: u1 into 2-3 of the seven N<x,v>, plus a raw entry with the opposite relative sign
        summ = [(rng.choice([1, -1]), G.lmul(b, u1)) for (b, T) in R[:rng.choice([2, 3])]]
        v = rng.choice(pts); b = rng.choice(ft.geo_rep_letters[v])
        summ.append((1, G.lmul(b, G.inp(0)))); summ.append((-s, G.lmul(b, G.inp(1))))
        gens.append(ft.move(summ))
    elif kind == "FE":
        # mirrored end gadget: private chains to depth m-2 for both arrays, then shared route sums (3-route
        # sums onto one element when possible) into the frames next to the kernel frame
        for l in range(W):
            cur = l
            for d in range(m - 2):
                g = ft.lead_terms(gens[cur])[0][1]
                R = ft.routes(g)
                if not R:
                    return None
                b = rng.choice(list(R.values()))[0]
                gens.append(ft.move([(1, G.lmul(b, gens[cur]))])); cur = len(gens) - 1
    while len(gens) - W < k:
        summ = route_sum(ft, gens, rng)
        if summ is None:
            return None
        u = ft.move(summ)
        if u is None:
            return None
        gens.append(u)
    return gens

def complete(ft, gens, l, maxE, st):
    G = ft.G
    r, sol = G.span_test(gens, l, 1500)
    if r == "yes":
        return 0
    Es = []
    try:
        Es = G.E_candidates(gens, l)[:maxE]
    except Exception:
        Es = []
    Es += junk_E(ft, gens, l, maxE)
    for E in Es:
        st["E_tested"] += 1
        r, sol = G.span_test(gens, l, 1500, E)
        if r == "cap":
            st["cap"] += 1
        if r == "yes":
            return 1
    return None

def _deadline():
    import os
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'DEADLINE')
    try:
        return float(open(f).read().split()[0])
    except Exception:
        return None

def _budget():
    import os
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BUDGET')
    try:
        return float(open(f).read().split()[0])
    except Exception:
        return None

def main():
    m, kind, k, n, seed = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    maxE = int(sys.argv[6]) if len(sys.argv) > 6 else 10
    out = sys.argv[7] if len(sys.argv) > 7 else None
    W = int(__import__('os').environ.get('FANO_W', '2'))
    t0 = time.time(); rng = random.Random(seed)
    ft = setup(m, seed)
    st = {"prefixes": 0, "fail_build": 0, "E_tested": 0, "cap": 0, "x1": 0, "x2": 0, "words": {}}
    hits = []
    dl = _deadline()
    for it in range(n):
        bu = _budget()
        if (dl is not None and time.time() > dl) or (bu is not None and time.time() - t0 > bu):
            st['stopped_at_deadline_after'] = it; break
        if it and it % 10 == 0:
            print('progress', it, st, flush=True)
        if len(ft.G.MATS) > 40000:
            ft = setup(m, seed + 1000 + it)
        gens = prefix(ft, kind, W, k, rng)
        if gens is None:
            st["fail_build"] += 1; continue
        st["prefixes"] += 1
        extra = []
        for l in range(W):
            e = complete(ft, gens, l, maxE, st)
            extra.append(e)
        if extra[0] is not None:
            st["x1"] += 1
        if extra[1] is not None:
            st["x2"] += 1
        for l in range(2, W):
            if extra[l] is not None:
                st["x%d" % (l + 1)] = st.get("x%d" % (l + 1), 0) + 1
        if None not in extra:
            tot = k + sum(extra)
            st["words"][tot] = st["words"].get(tot, 0) + 1
            hits.append((it, tot))
    st["seconds"] = round(time.time() - t0, 1)
    print(f"fano7 m={m} {kind} k={k} n={n} seed={seed}: {st}", flush=True)
    if out:
        json.dump({"args": sys.argv[1:], "stats": st, "hits": hits}, open(out, "w"))

if __name__ == "__main__":
    main()
