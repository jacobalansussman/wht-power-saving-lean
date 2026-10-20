"""staggered.py -- ms-fano: EXACT search of STAGGERED-LINEAGE templates (the shared-base / halving shape of the Fano
blueprint, cut down to its smallest form) at block m, 2 arrays.
Template ST(L, c): one shared lineage of L = m+1 moves
     G1 = D b1 x1 ;  G2 = D (r2 G1 + s b2 x2) ;  G_t = D r_t G_{t-1}  (t = 3..L)
with every r_t a DOUBLE ADVANCE (both inputs' lead terms one step further on a geodesic, as long as possible), so
x1 reaches the kernel frame at move m and x2 at move m+1; then c CLEANUP moves, each D(s1 r G_i [+ s2 r' G_j])
with routes or exact alignments onto existing junk terms; final test = exact span test (closure over frames, exact
solve over Q, so the final combination may halve).  Moves = m + 1 + c: block 3 c = 1 is 2-in-5 (excluded by 14W/5,
must give 0), c = 2 is 2-in-6 (plain count); block 4 c = 2 is 2-in-7 (the target), c = 3 is 2-in-8 (plain count).
usage: python3 staggered.py m c n seed [out.json]"""
import sys, time, json, random
from tsearch import *

def double_advance(ft, gen, rng, samples=6000):
    """free letters b such that D b advances the deepest geodesic term of each input that is not yet at the kernel
    frame (and keeps it on a geodesic); returns a list (possibly empty)"""
    G = ft.G
    leads = {}
    for (l, g) in ft.lead_terms(gen):
        d0, dk = ft.depth(g)
        if dk == 0:
            continue
        if l not in leads or ft.depth(leads[l])[0] < d0:
            leads[l] = g
    if not leads:
        return []
    out = []
    for _ in range(samples):
        b = ft.geo.rand_free()
        ok = True
        for l, g in leads.items():
            d0, dk = ft.depth(g)
            h = G.intern(G.mul(G.D1, b, G.MATS[g]))
            if ft.depth(h) != (d0 + 1, dk - 1):
                ok = False; break
        if ok:
            out.append(b)
            if len(out) >= 4:
                break
    return out

def cleanup_content(ft, gens, rng, maxs=2):
    """1-2 summands; the first routes or aligns a random term of a random earlier generator onto a junk term of
    some generator (exact alignment when the needed letter is free), else a random route"""
    G = ft.G
    i = rng.randrange(2, len(gens))          # skip raw inputs as the first summand
    terms = list(gens[i].keys())
    (l, g) = rng.choice(terms)
    b = None
    if rng.random() < 0.6:
        # try to align D b g onto a term of another generator (same input): cancellation material
        cands = [(j, t) for j in range(2, len(gens)) if j != i for t in gens[j] if t[0] == l]
        rng.shuffle(cands)
        for (j, (l2, g2)) in cands[:12]:
            bb = ft.align(g2, g)
            if bb is not None:
                b = bb; break
    if b is None:
        R = ft.routes(g) if ft.on_geo(g) else {}
        if R:
            b = rng.choice(list(R.values()))[0]
        else:
            b = ft.geo.rand_free()
    summ = [(1, G.lmul(b, gens[i]))]
    if maxs > 1 and rng.random() < 0.6:
        j = rng.randrange(0, len(gens))
        if j != i:
            tj = list(gens[j].keys()); (l3, g3) = rng.choice(tj)
            T = G.intern(G.mul(b, G.MATS[g]))   # put the second summand's term on top of the first's (cancel or add)
            bb = None
            # second summand: b' with b' g3 == b g  (same element before D): b' = (b g) g3^-1, free iff monomial
            cand = G.mul(G.MATS[T], np.conj(G.MATS[g3]).T)
            if G.is_monomial(cand):
                bb = cand
            if bb is None:
                bb = ft.geo.rand_free()
            summ.append((rng.choice([1, -1]), G.lmul(bb, gens[j])))
    return summ

def build(ft, m, c, rng):
    G = ft.G
    reps = ft.geo_rep_letters
    gens = [G.inp(0), G.inp(1)]
    v1 = rng.randrange(1, 1 << m); b1 = rng.choice(reps[v1])
    gens.append(ft.move([(1, G.lmul(b1, gens[0]))]))
    prog = [("G1", v1)]
    for t in range(2, m + 2):
        prev = gens[-1]
        if t == 2:
            R = list(ft.routes(ft.lead_terms(prev)[0][1]).values())
            r2 = rng.choice(R)[0]
            v2 = rng.randrange(1, 1 << m); b2 = rng.choice(reps[v2]); s = rng.choice([1, -1])
            u = ft.move([(1, G.lmul(r2, prev)), (s, G.lmul(b2, gens[1]))])
            prog.append(("G2", v2, s))
        else:
            B = double_advance(ft, prev, rng)
            if not B:
                return None
            u = ft.move([(1, G.lmul(rng.choice(B), prev))])
            prog.append(("G%d" % t, "double"))
        if u is None:
            return None
        gens.append(u)
    for k in range(c):
        summ = cleanup_content(ft, gens, rng)
        u = ft.move(summ)
        if u is None:
            return None
        gens.append(u); prog.append(("clean", len(summ)))
    return gens, prog

def main():
    m, c, n, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    out = sys.argv[5] if len(sys.argv) > 5 else None
    t0 = time.time(); rng = random.Random(seed)
    ft = setup(m, seed)
    st = {"built": 0, "fail_build": 0, "hits": 0, "x1_only": 0, "x2_only": 0, "cap": 0}
    hits = []
    for it in range(n):
        if len(ft.G.MATS) > 40000:
            ft = setup(m, seed + 1000 + it)
        r = build(ft, m, c, rng)
        if r is None:
            st["fail_build"] += 1; continue
        gens, prog = r
        st["built"] += 1
        res = []
        for l in range(2):
            s_, sol = ft.G.span_test(gens, l, 1500)
            res.append(s_)
        if "cap" in res:
            st["cap"] += 1
        if res == ["yes", "no"] or (res[0] == "yes" and res[1] != "yes"):
            st["x1_only"] += 1
        if res[1] == "yes" and res[0] != "yes":
            st["x2_only"] += 1
        if res == ["yes", "yes"]:
            st["hits"] += 1; hits.append(prog)
    st["seconds"] = round(time.time() - t0, 1)
    print(f"ST m={m} c={c} moves={m+1+c} n={n} seed={seed}: {st}", flush=True)
    if hits:
        print("hits:", hits[:5])
    if out:
        json.dump({"args": sys.argv[1:], "stats": st, "hits": hits[:50]}, open(out, "w"))

if __name__ == "__main__":
    main()
