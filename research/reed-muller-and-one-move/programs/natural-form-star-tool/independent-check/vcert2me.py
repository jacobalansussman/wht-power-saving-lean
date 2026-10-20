"""vcert2me.py (v-rm2-exact): own line-by-line check of starrec2 CERT=1 output lines (C, M, Z; round 1).
Each cited zero cell x,y:F/S/round/d/i: r_F = 0; S <= D_F; if d >= 0: S = D_F n A for a d-flat A on which OWN round-1
recursion of F alone gives u(A) = 0; if d = -1 (quadratic star, i = G1*65536 + G2): S = D_F n D_G1 n D_G2 and
s_G1 + s_G2 + min(|D_F n D_G3|, s_G3) - 2|S| <= 0 with G3 = F + G1 + G2 (own slacks). Then:
C: every cell of Y x Y cited and |Y| > W.  M: matching of D_Q x D_Q avoiding the cited cells < s_Q (own s_Q).
Z: own round-1 recursion of F reaches a negative right side."""
import sys, numpy as np, vf, vrec
from collections import Counter
def rec_u(F, w, nu):
    P = vrec.pre(F); u = []
    for dd in range(5):
        mw, sz, S = P[dd]; u.append(np.minimum(sz, (mw - w[None, :]).min(axis=1)).astype(np.int64))
    for dd in range(3, -1, -1):
        g = (1 << (5 - dd)) - 1; sz = P[dd][1]; Bu = u[dd + 1][vrec.UP[dd]]
        val = Bu.sum(axis=1) - 2 * sz + 2 * nu
        if ((val < 0) & (sz > 0)).any(): return None, P, dd
        u[dd] = np.where(sz > 0, np.minimum(u[dd], np.minimum(val // (g - 2), Bu.min(axis=1))), u[dd])
    return u, P, None
pc = lambda x: bin(x).count('1')
def check(line):
    t = line.split(); ks = [int(x) for x in t[1].split(',')]; kind = t[2]; W = len(ks)
    w = vrec.CLSPAR[:, ks].sum(axis=1); r = w[vrec.VID]; s = vf.WT - r
    if kind == 'Z':
        F, d = int(t[4]), int(t[5]); assert r[F] == 0 and int(t[3]) == 1
        u, P, dd = rec_u(F, w, 0); return dd is not None
    cache = {}; Z = set(); uc = {}
    if kind == 'C': cells = t[7:]
    elif kind == 'M': cells = t[7:]
    else: return None
    for c in cells:
        xy, rest = c.split(':'); x, y = map(int, xy.split(',')); fs = list(map(int, rest.split('/')))
        if len(fs) == 3: fs += [None, None]
        F, S, rnd, d, i = fs
        if rnd != 1: return 'round>1'
        key = (F, S, d, i)
        if key not in cache:
            if r[F] != 0 or S & ~int(vf.MASK[F]): cache[key] = False
            elif d is None or d >= 0:
                if F not in uc: uc[F] = rec_u(F, w, 0)
                u, P, dd = uc[F]
                dims = range(5) if d is None else [d]
                cache[key] = u is not None and any(int(P[e][2][j]) == S and u[e][j] == 0 for e in dims for j in range(len(P[e][2])))
            else:
                G1, G2 = i // 65536, i % 65536; G3 = F ^ G1 ^ G2; A = int(vf.MASK[F])
                Sq = A & int(vf.MASK[G1]) & int(vf.MASK[G2])
                S3 = (A & int(vf.MASK[G3])) if G3 else 0
                s3 = pc(S3) if not G3 or G3 == F else min(pc(S3), int(s[G3]))
                cache[key] = (Sq == S) and (int(s[G1]) + int(s[G2]) + s3 - 2 * pc(S) <= 0)
        if not cache[key]: return 'BLOCK-FAIL %s' % (key,)
        if not (S >> (x - 1) & 1 and S >> (y - 1) & 1): return 'CELL-FAIL'
        Z.add((x, y))
    if kind == 'C':
        q = int(t[4]); Y = int(t[6][2:]); Yd = [z for z in range(1, 32) if Y >> (z - 1) & 1]
        return len(Yd) == q > W and all((a, b) in Z for a in Yd for b in Yd)
    if kind == 'M':
        Q = int(t[4]); Dq = [z for z in range(1, 32) if int(vf.MASK[Q]) >> (z - 1) & 1]
        mm = vrec.matching(Dq, Z)
        return mm < int(s[Q]) and int(s[Q]) == int(t[6])
if __name__ == "__main__":
    cnt = Counter(); bad = []
    for l in open(sys.argv[1]):
        if not l.startswith('S '): continue
        v = check(l); k = l.split()[2]; cnt[(k, str(v))] += 1
        if v is not True: bad.append((v, l[:120]))
    print(sys.argv[1], dict(cnt)); [print("NOT VERIFIED", b) for b in bad[:10]]
