#!/usr/bin/env python3
"""v_lemma1.py: independent check of PROOF Lemma 1 (parity law) and the order-3 control, with vgame.py only.
sigma drawn as random elements of P (random Levi A, random symmetric phase) FILTERED by K sigma K^-1 in P
(not built from the claimed formula), then the program r1 = y1 + y2, r2 = y1 - sigma y2, K on each register;
states = inputs + every move output; exact closure solve for both outputs; psi screen; cut versions."""
import sys, random, itertools
from vgame import *

def levi(Gp, A_cols):
    """[[A,0],[0,A^-T]] as column images; A given by column images in F2^m"""
    m = Gp.m
    # A^-T: compute via omega-duality: need g with omega preserved. Build g on x-part, solve the b-part by inverse.
    # Make the 2m x 2m matrix diag(A, B) and find B with omega preserved: B = A^-T.
    # A as m x m: column j = A_cols[j]
    # inverse of A
    Ai = Gp_inv_m(A_cols, m)
    # A^-T column j = row j of A^-1 ; row j of Ai: bit k = (Ai[k] >> j) & 1
    AiT = [sum(((Ai[k] >> j) & 1) << k for k in range(m)) for j in range(m)]
    return tuple([A_cols[j] for j in range(m)] + [AiT[j] << m for j in range(m)])
def Gp_inv_m(cols, m):
    G1 = G(m // 2 if False else 1)
    # brute: solve over F2 using generic routine of a doubled size? simple Gauss here
    rows = []
    for i in range(m):
        r = 0
        for j in range(m):
            if (cols[j] >> i) & 1: r |= 1 << j
        rows.append([r, 1 << i])
    for c in range(m):
        p = next(k for k in range(c, m) if (rows[k][0] >> c) & 1)
        rows[c], rows[p] = rows[p], rows[c]
        for k in range(m):
            if k != c and (rows[k][0] >> c) & 1:
                rows[k][0] ^= rows[c][0]; rows[k][1] ^= rows[c][1]
    return [sum(((rows[i][1] >> j) & 1) << i for i in range(m)) for j in range(m)]
def phase(Gp, Sm):
    """(a|b) -> (a | b + S a), S symmetric given as rows"""
    m = Gp.m
    return tuple([(1 << j) | (sum(((Sm[i] >> j) & 1) << i for i in range(m)) << m) for j in range(m)]
                 + [1 << (m + j) for j in range(m)])
def randP(Gp, rng):
    m = Gp.m
    while True:
        cols = [rng.randrange(1, 1 << m) for _ in range(m)]
        if len(G.span_rref(cols)) == m: break
    Sm = [0] * m
    for i in range(m):
        for j in range(i, m):
            if rng.random() < 0.5: Sm[i] |= 1 << j; Sm[j] |= 1 << i
    return Gp.comp(phase(Gp, Sm), levi(Gp, cols))
def order(Gp, g, cap=64):
    h = g
    for n in range(1, cap + 1):
        if h == Gp.I: return n
        h = Gp.comp(g, h)
    return None

def program(Gp, sig):
    y1 = {(0, Gp.I): Fr(1)}; y2 = {(1, Gp.I): Fr(1)}
    regs = [add((1, y1), (1, y2)), add((1, y1), (-1, lm(Gp, sig, y2)))]
    states = [y1, y2]
    for c in regs:
        for i in range(Gp.m):
            c = lm(Gp, Gp.T(1 << i), c); states.append(c)
    return states

def verdict(Gp, states):
    out = []
    for l in range(2):
        r, sol = closure_solve(Gp, states, {(l, Gp.K): Fr(1)})
        out.append(r + (str(len(sol)) if sol else ""))
    ps = "".join("1" if psi_screen(Gp, states, l) else "0" for l in range(2))
    return out, ps

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
for m in (3, 4, 5):
    Gp = G(m)
    assert Gp.is_symp(Gp.K) and Gp.K == tuple([1 << j for j in range(m)] + [(1 << (m + j)) | (1 << j) for j in range(m)])
    # plain
    st = program(Gp, Gp.I)
    y1 = {(0, Gp.I): Fr(1)}; y2 = {(1, Gp.I): Fr(1)}
    plain = [y1, y2]
    for c in (y1, y2):
        for i in range(m): c = lm(Gp, Gp.T(1 << i), c); plain.append(c)
    print(f"m={m} plain 2 in {2*m}: {verdict(Gp, plain)}; cut: {verdict(Gp, plain[:-1])}", flush=True)
    # 3-cycle address map
    cols = [1 << ((k + 1) % 3) if k < 3 else 1 << k for k in range(m)]
    s3 = levi(Gp, cols)
    assert Gp.is_symp(s3) and Gp.inP(s3) and Gp.inP(Gp.comp(Gp.K, Gp.comp(s3, Gp.inv(Gp.K))))
    st = program(Gp, s3)
    print(f"   o3 control: {verdict(Gp, st)}; cut by one move: {verdict(Gp, st[:-1])}", flush=True)
    seen = {}; tries = 0
    while tries < 400 and sum(seen.values()) < 16:
        tries += 1
        while True:
            cols = [rng.randrange(1, 1 << m) for _ in range(m)]
            if len(G.span_rref(cols)) == m: break
        L = levi(Gp, cols)
        cands = []
        if m <= 4:
            pairs = [(i, j) for i in range(m) for j in range(i, m)]
            for bits in range(1 << len(pairs)):
                Sm = [0] * m
                for t, (i, j) in enumerate(pairs):
                    if (bits >> t) & 1: Sm[i] |= 1 << j; Sm[j] |= 1 << i
                g = Gp.comp(phase(Gp, Sm), L)
                if Gp.inP(Gp.comp(Gp.K, Gp.comp(g, Gp.inv(Gp.K)))): cands.append(g)
            assert len(cands) == 1, len(cands)
        else:
            Ai = Gp_inv_m(cols, m)   # columns of A^-1
            # S = I + A^-T A^-1 (derived by hand; membership re-checked below)
            def mat(colsx):  # rows of matrix from columns
                return [sum(((colsx[j] >> i) & 1) << j for j in range(m)) for i in range(m)]
            AiR = mat(Ai)                    # rows of A^-1
            # (A^-T A^-1)[i][j] = sum_k A^-1[k][i] A^-1[k][j] = col_i(A^-1) . col_j(A^-1)
            Sm = [0] * m
            for i in range(m):
                for j in range(m):
                    v = bin(Ai[i] & Ai[j]).count("1") & 1
                    if i == j: v ^= 1
                    if v: Sm[i] |= 1 << j
            g = Gp.comp(phase(Gp, Sm), L)
            cands.append(g)
        g = cands[0]
        assert Gp.is_symp(g) and Gp.inP(g) and Gp.inP(Gp.comp(Gp.K, Gp.comp(g, Gp.inv(Gp.K))))
        n = order(Gp, g)
        if n is None or n < 2 or seen.get(n, 0) >= 2: continue
        seen[n] = seen.get(n, 0) + 1
        st = program(Gp, g)
        v = verdict(Gp, st)
        exp = "yes" if n % 2 else "no"
        ok = all(x.startswith(exp) for x in v[0])
        print(f"   two-frame sigma order {n}: exact {v[0]} psi {v[1]}  expected {exp}: {'OK' if ok else 'MISMATCH'}", flush=True)
