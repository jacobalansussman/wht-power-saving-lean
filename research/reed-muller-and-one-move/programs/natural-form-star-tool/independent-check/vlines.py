"""vlines.py: verify scanrc lines (D, X ten-piece, Y one-level star) from scratch with vf.py tables.
Every bound is recomputed: min over ALL nonzero forms G with D_G >= piece of s_G = |D_G| - r_G (and the piece size)."""
import sys, numpy as np, vf
MASK, WT, PAR = vf.MASK, vf.WT, vf.PAR
def lin(c): return sum(1 << (z - 1) for z in range(1, 32) if bin(c & z).count('1') & 1)
LIN = [0] + [lin(c) for c in range(1, 32)]
FULL = (1 << 31) - 1
def pc(x): return bin(x).count('1')
def minslack(S, s):
    sel = (MASK[1:] & S) == S
    return int(s[1:][sel].min()) if sel.any() else 10**6
def check(line):
    t = line.split(); ks = [int(x) for x in t[1].split(',')]; kind = t[2]
    r = PAR[:, ks].sum(axis=1).astype(np.int64); s = WT - r
    if kind == 'D': return bool((s[1:] < 0).any())
    if kind == 'P': return None
    F, ph, ps = int(t[3]), int(t[4]), int(t[5]); A = int(MASK[F]); sF = int(s[F])
    if kind == 'X':
        P0, P1 = A & ~LIN[ph] & FULL, A & LIN[ph]
        def blk(R, C): return min(pc(R), pc(C), minslack(R | C, s))
        q = lambda P, c: P & (LIN[ps] if c else ~LIN[ps] & FULL)
        bound = blk(P0, P0) + blk(P1, P1)
        for a in (0, 1):
            for b in (0, 1):
                bound += 2 * blk(q(P0, a), q(P1, b))
        # sanity: the ten blocks partition D_F x D_F (cells counted)
        cells = pc(P0) ** 2 + pc(P1) ** 2 + 2 * sum(pc(q(P0, a)) * pc(q(P1, b)) for a in (0, 1) for b in (0, 1))
        assert cells == pc(A) ** 2
        return bound < sF and sF > 0
    if kind == 'Y':
        c1, c2 = int(t[6]), int(t[7])
        if r[F] != 0: return False
        h = lambda c, v: A & (LIN[c] if v else ~LIN[c] & FULL)
        Psi = h(ph, c1) & h(ps, c2)
        H = [h(ph, c1), h(ps, c2), h(ph ^ ps, c1 ^ c2)]
        pieces = [x & ~Psi for x in H]
        assert all(Psi & ~x == 0 for x in H)
        assert pieces[0] | pieces[1] | pieces[2] == A & ~Psi and sum(pc(p) for p in pieces) == pc(A & ~Psi)
        tot = sum(min(pc(x), minslack(x, s)) for x in H)
        return 2 * pc(Psi) > tot
if __name__ == "__main__":
    import random; random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    L = [l for l in open(sys.argv[1]) if l.startswith('S ')]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else len(L)
    L = random.sample(L, min(n, len(L)))
    from collections import Counter; ok = Counter(); bad = []
    for l in L:
        v = check(l); k = l.split()[2]
        ok[(k, v)] += 1
        if v is False: bad.append(l)
    print(sys.argv[1], dict(ok)); [print("BAD", b.strip()) for b in bad[:5]]
