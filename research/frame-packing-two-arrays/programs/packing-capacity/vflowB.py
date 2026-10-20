#!/usr/bin/env python3
"""vflowB.py (fp-capacity): Theorem A with coarser classes, raw gauge (PROOF.md section 3a).

Classes modulo a subgroup H of the free group B (term g y_l  ->  class g H):
  mode "img":  H = B. The class of g is its IMAGE FRAME F = gbar L0 (G/B = Sp/P = Lagrangians). The quotient graph
               B\\G/B is the level path 0 - 1 - ... - m.  Claim: for every level step i-1 -> i, the number of
               (paid move, live image frame of the moved array at level i-1 that the move ascends) is >= W.
               Hence Delta <= P_img_asc := sum_moves (live image frames ascended - 1).
  mode "coimg": gauge h_a = K^-1, H = B: the class of g is its CO-IMAGE FRAME gbar L1; co-level d(gbar L1, L0) runs
               from m to 0; every co-level step i -> i-1 must be crossed >= W times by live co-image frames.
  mode "opp":  H = Paulis x GL(m,2) (address maps).  The class of g is the opposite pair (gbar L0, gbar L1); the
               weight of a virtual move is max(0, Delta f)/2 with f = d(F, L0) - d(F', L0).
               Claim: sum of weights over live virtual moves >= W m; Delta <= P_opp_w.
Also reports the frame-label (Sp) live packing for comparison.  usage: python3 vflowB.py [seed] [n_random]"""
import sys, os, random, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "ledger"))
import vflow
from vflow import *          # Engine, run, cases_for, random_order3, ...
from vflow import _step, analyse, setup, cases_for, random_order3

def lag_img(E, S, base):
    """bitset of the Lagrangian S . base, base given as a list of 2m-vectors (columns)"""
    m = E.m
    cols = [(S @ b) % 2 for b in base]
    out = 0
    for mask in range(1 << m):
        v = np.zeros(2 * m, dtype=np.int64)
        for j in range(m):
            if (mask >> j) & 1: v ^= cols[j]
        idx = 0
        for k in range(2 * m):
            if v[k] & 1: idx |= 1 << k
        out |= 1 << idx
    return out

def run_mode(E, W, R, word, mode, labs_index, Dm):
    m = E.m
    L0b = [np.array([0] * m + [1 if i == j else 0 for i in range(m)], dtype=np.int64) for j in range(m)]   # (0|e_j)
    L1b = [(E.Ksym @ b) % 2 for b in L0b]
    i0 = labs_index[lag_img(E, np.eye(2 * m, dtype=np.int64), L0b)]
    cache = {}
    def key(gid):
        if gid not in cache:
            S = E.syms[gid]
            F = labs_index[lag_img(E, S, L0b)]
            if mode == "img": cache[gid] = (F,)
            elif mode == "coimg": cache[gid] = (labs_index[lag_img(E, S, L1b)],)
            else: cache[gid] = (F, labs_index[lag_img(E, S, L1b)])
        return cache[gid]
    def level(k): return int(Dm[i0, k[0]])
    def fval(k): return int(Dm[i0, k[0]] - Dm[i0, k[1]])
    start = {r: r for r in range(W)}
    C = [dict() for _ in range(R)]
    for r, l in start.items(): C[r] = {(l, E.ID): ONE}
    snaps = [C]
    for st in word:
        C = _step(E, C, st); snaps.append(C)
    res = run(E, W, R, word)
    if not res["unit"]: return dict(unit=False)
    def classes(c):
        out = {}
        for (l, gid) in c: out.setdefault(key(gid), []).append((l, gid))
        return out
    def img_after(k, name, gid_example):
        nid, _ = E.lmul(name, gid_example); return key(nid)
    T = len(word)
    live = [None] * (T + 1); live[T] = set()
    for r in range(W):
        for k in classes(snaps[T][r]): live[T].add((r, k))
    for t in range(T - 1, -1, -1):
        st = word[t]; cur = [classes(c) for c in snaps[t]]
        L = set(); touched = set(); kind = st[0]
        if kind in ("move", "free"):
            r = st[1]; touched.add(r); name = ("D%d" % st[2]) if kind == "move" else st[2]
            for k, terms in cur[r].items():
                if (r, img_after(k, name, terms[0][1])) in live[t + 1]: L.add((r, k))
        elif kind == "mix":
            r, terms_ = st[1], st[2]; touched.add(r)
            for (sreg, mu, name) in terms_:
                for k, terms in cur[sreg].items():
                    k2 = k if name is None else img_after(k, name, terms[0][1])
                    if (r, k2) in live[t + 1]: L.add((sreg, k))
        elif kind == "gate":
            regs, A = st[1], st[2]
            for i, r in enumerate(regs):
                touched.add(r)
                for j, s in enumerate(regs):
                    if gq(A[i][j]).iszero(): continue
                    for k in cur[s]:
                        if (r, k) in live[t + 1]: L.add((s, k))
        for r in range(R):
            if r in touched: continue
            for k in cur[r]:
                if (r, k) in live[t + 1]: L.add((r, k))
        live[t] = L
    M = res["M"]; Delta = W * m - M
    cross = [0] * (m + 1)        # cross[i] = live ascents from level i-1 to i
    P_live = 0; P_asc = 0; wsum = 0.0; P_w = 0.0; maxk = 0
    for t, st in enumerate(word):
        if st[0] != "move": continue
        r = st[1]; cur = classes(snaps[t][r]); name = "D%d" % st[2]
        kl = [(k, terms) for k, terms in cur.items() if (r, k) in live[t]]
        maxk = max(maxk, len(kl))
        P_live += max(0, len(kl) - 1)
        a = 0; ww = 0.0
        for k, terms in kl:
            k2 = img_after(k, name, terms[0][1])
            if mode == "coimg":
                if level(k2) < level(k): a += 1; cross[level(k)] += 1     # co-level step level(k) -> level(k)-1
            elif level(k2) > level(k): a += 1; cross[level(k2)] += 1
            if mode == "opp": ww += max(0, fval(k2) - fval(k)) / 2.0
        P_asc += a - 1; wsum += ww; P_w += ww - 1
    out = dict(unit=True, mode=mode, M=M, Wm=W * m, Delta=Delta, P_live=P_live, P_asc=P_asc, maxk=maxk,
               cross=cross[1:], ok_cross=all(c >= W for c in cross[1:]), ok_asc=(Delta <= P_asc))
    if mode == "opp": out.update(wsum=wsum, P_w=P_w, ok_w=(wsum >= W * m - 1e-9 and Delta <= P_w + 1e-9))
    return out

def crossing_unit(E, rng, W=None):
    """random 'crossing' unit (no saving, live image-frame packing): input a moves along its own prefix set P_a,
    a random invertible Gaussian-rational gate G mixes all W registers, every register moves along the common set C,
    the inverse gate G^-1 unmixes, input a moves along its suffix set Q_a; P_a, C, Q_a partition the m directions."""
    m = E.m
    if W is None: W = rng.choice([2, 3])
    dirs = list(range(m)); rng.shuffle(dirs)
    c = rng.randint(1, m - 1); Cset = dirs[:c]; rest = dirs[c:]
    w = []
    Q = []
    for a in range(W):
        Pa = [d for d in rest if rng.random() < 0.5]; Qa = [d for d in rest if d not in Pa]
        w += [("move", a, 1 << d) for d in Pa]; Q.append(Qa)
    import t_random
    while True:
        A = [[GQ(rng.randint(-2, 2), rng.randint(-1, 1)) for _ in range(W)] for _ in range(W)]
        Ai = t_random.inv_gq_matrix(A)
        if Ai is not None: break
    w.append(("gate", list(range(W)), A))
    mv = [("move", r, 1 << d) for r in range(W) for d in Cset]; rng.shuffle(mv); w += mv
    w.append(("gate", list(range(W)), Ai))
    for a in range(W): w += [("move", a, 1 << d) for d in Q[a]]
    return ("random crossing unit W=%d |C|=%d" % (W, c), W, W, w)

def crossing_pair(E):
    """r0 = D1 x0 + x1, r1 = D1 x0 - x1; both moved along D2..Dm; (r0 + r1)/2 = K x0; (r0 - r1)/2 then D1 = K x1"""
    m = E.m; h = Fr(1, 2)
    w = [("move", 0, 1), ("mix", 2, [(0, 1, None)]),
         ("mix", 0, [(2, 1, None), (1, 1, None)]), ("mix", 1, [(2, 1, None), (1, -1, None)])]
    for j in range(1, m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
    w += [("gate", [0, 1], [[h, h], [h, -h]]), ("move", 1, 1)]
    return ("crossing pair", 2, 3, w)

def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    nrand = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    rng = random.Random(seed)
    rows = []; bad = []
    for m in (2, 3, 4):
        E = Engine(m)
        labs, L0bits, Dm = setup(m)
        labs_index = {x: i for i, x in enumerate(labs)}
        cs = cases_for(E, rng) + [crossing_pair(E)] + [crossing_unit(E, rng) for _ in range(nrand)]
        if m >= 3: cs += [random_order3(E, rng) for _ in range(nrand)]
        for (name, W, R, w) in cs:
            fr = analyse(E, W, R, w, gauge_trivial(E, W), labs_index, Dm, raw=True)
            for mode in ("img", "coimg", "opp"):
                o = run_mode(E, W, R, w, mode, labs_index, Dm)
                o.update(case=name, m=m, W=W, P_live_frame=fr.get("P_live"), P_asc_frame=fr.get("P_asc"))
                rows.append(o)
                okk = o["unit"] and o["ok_cross"] and o["ok_asc"] and o.get("ok_w", True)
                if not okk: bad.append(o)
                print("m=%d %-28s %-3s M=%d Delta=%d  P_live: img/opp=%d frame=%s | P_asc=%d (frame %s) | crossings per level %s %s"
                      % (m, name, mode, o["M"], o["Delta"], o["P_live"], fr.get("P_live"), o["P_asc"], fr.get("P_asc"),
                         o["cross"], ("P_w=%.1f" % o["P_w"]) if mode == "opp" else ""), flush=True)
    print("\nSUMMARY: %d rows, %d violations" % (len(rows), len(bad)))
    with open(os.path.join(HERE, "vflowB_s%d.json" % seed), "w") as fh: json.dump(dict(rows=rows, bad=bad), fh, indent=1, default=str)

if __name__ == "__main__":
    main()
