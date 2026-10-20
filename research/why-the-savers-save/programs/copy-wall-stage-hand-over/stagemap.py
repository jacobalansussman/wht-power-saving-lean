"""stagemap.py (d-g-crossstage) -- the stage map of one bank pair in the hub's B_K compressed-window format, and the
hand-over test.  Standard library + our own u-floor/copyfloor.py (E8 lines, F2 labels).  usage: python3 -I stagemap.py
(1) scalar: replay the B_K recipe of hub global_k.py (re-implemented, K = 4, s = 9) on symbolic bank vectors and check,
    for every pair of consecutive active stages of one bank pair, that the stage-(i+1) SOURCE bank holds
    (stage-i target before stage i) + (stage-i source) [+ a bridge term], exactly.
(2) frames: one port's registers at the boundary i -> i+1, written in window-i (old) and window-(i+1) (new)
    identification, from the hub's boundary identities H_i = P_(i+1), K_i = 0_(i+1) (five-stage geometry527.py asserts).
    Model: old coordinates U (+) B relative to the common base E_i (U = F2^9 dot form, B = the next 8-block);
    new identification of port T: t_T(a_T^perp) = Enew, t_T(a_T) = a_T', t_T(B) = a_T'^perp in U'.
(3) the transitions t_T are port dependent; in the OLD identification no stage-(i+1) target knows any stage-(i+1)
    source (all 112 x 112 pairs, diagonal included).
(4) frame distances of the hand-over routes on the complex P112 frame decoder, against the star drop w_in."""
import sys, os, itertools
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-floor-checker')))
import copyfloor as cf

# ---------------- (1) scalar replay of the B_K recipe ----------------
def recipe(K):
    out = [dict(kind='local', stage=0, pair=0, reverse=False)]; stage = 1
    for j in range(K):
        out.append(dict(kind='local', stage=stage, pair=j, reverse=True)); stage += 1
        if j < K - 1: out.append(dict(kind='bridge', dst=j, src=j + 1, c=1))
        out.append(dict(kind='local', stage=stage, pair=j, reverse=False)); stage += 1
        if j < K - 1: out.append(dict(kind='bridge', dst=j, src=j + 1, c=-1))
    for j in reversed(range(K - 1)): out.append(dict(kind='bridge', dst=j + 1, src=j, c=1))
    out.append(dict(kind='exchange'))
    return out

def scalar(K=4):
    n = 2 * K
    z = [[int(i == j) for j in range(n)] for i in range(n)]      # bank b holds sum_j z[b][j] * (initial bank j)
    def add(d, s, c): z[d] = [x + c * y for x, y in zip(z[d], z[s])]
    log = []; last = {}
    for e in recipe(K):
        if e['kind'] == 'local':
            p, rev = e['pair'], e['reverse']
            src, tgt = (2 * p + 1, 2 * p) if rev else (2 * p, 2 * p + 1)   # P-bank = 2p, 0-bank = 2p+1
            before_tgt, before_src = z[tgt][:], z[src][:]
            add(tgt, src, -1 if rev else 1)
            if p in last:   # previous active stage of this pair: check source(i+1) = target(i) after stage i (+ bridge)
                pi = last[p]
                ok_same = z[src] == before_src and before_src == pi['tgt_after']
                diff = [a - b for a, b in zip(before_src, pi['tgt_after'])]
                log.append((pi['stage'], e['stage'], p, 'rev' if pi['rev'] else 'fwd', 'rev' if rev else 'fwd',
                            pi['tgt'] == src, before_src == pi['tgt_after'], [k for k, x in enumerate(diff) if x]))
            last[p] = dict(stage=e['stage'], rev=rev, tgt=tgt, tgt_after=z[tgt][:], tgt_before=before_tgt, src_content=before_src)
        elif e['kind'] == 'bridge':
            d, s, c = e['dst'], e['src'], e['c']
            add(2 * d, 2 * s, c); add(2 * s + 1, 2 * d + 1, -c)
        else:
            for j in range(K): z[2 * j], z[2 * j + 1] = z[2 * j + 1], [-x for x in z[2 * j]]
    ident = z == [[int(i == j) for j in range(n)] for i in range(n)]
    return log, ident

# ---------------- F2 helpers ----------------
def rank(vs):
    basis = {}
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in basis: v ^= basis[h]
            else: basis[h] = v; break
    return len(basis)
def dot(a, b): return bin(a & b).count('1') & 1
def perp(vs, n):   # orthogonal complement in F2^n (dot form), as a list of basis ints
    cons = [v for v in vs]
    sol = []
    # brute force for n <= 9 over the space: fine (512)
    allv = [x for x in range(1, 1 << n) if all(dot(x, c) == 0 for c in cons)]
    b = []
    for x in allv:
        if rank(b + [x]) > len(b): b.append(x)
    return b
def dim(vs): return rank(vs)
def meet_dim(A, B): return dim(A) + dim(B) - dim(A + B)

def main():
    print('(1) scalar replay of the B_K recipe (hub global_k.py), K = 4, s = 9; P-bank = 2p, 0-bank = 2p+1')
    log, ident = scalar(4)
    print('    whole word returns to the identity (exchange included):', ident)
    for si, sj, p, oi, oj, same_bank, exact, extra in log:
        print('    pair %d stage %d (%s) -> %d (%s): source(i+1) bank = target(i) bank: %s; content = target(i) after stage i '
              '[= target before + source of stage i]: %s; extra banks in between (bridge): %s'
              % (p, si, oi, sj, oj, same_bank, exact, extra or '-'))
    # ---------------- (2) the frame table ----------------
    print('\n(2) frames of ONE port at the boundary i -> i+1 (dims relative to the common base E_i; h = 9, a = h - k = 8)')
    rows = [
     ('fwd i -> rev i+1', [
      ('P-bank (src i, tgt i+1)', 'i: <a_T> -> U', 'U  = <a_T>+E_(i+1)  [new: line]', 'i+1: U -> U+B  [new: <a_T> -> U\']'),
      ('0-bank (tgt i, src i+1)', 'i: 0 -> a_T^perp', 'a_T^perp = E_(i+1)  [new: frame 0]', 'i+1: a_T^perp -> a_T^perp+B  [new: 0 -> a_T\'^perp]'),
      ('stage-i total', 'formed at P_j (dim 8)', 'not a port register: no new identification', ''),
      ('stage-i copy (star)', 'P_j -> 0 (paid 8)', 'frame 0 of window i = E_i, 8 below E_(i+1)', ''),
      ('stage-(i+1) total', '', 'formed at <= P_j\'^perp (new)', 'copy climbs 8 to U\' (new)')]),
     ('rev i -> fwd i+1', [
      ('0-bank (src i, tgt i+1)', 'i: 0 -> a_T^perp', 'a_T^perp = E_(i+1)  [new: frame 0]', 'i+1: [new: 0 -> a_T\'^perp]'),
      ('P-bank (tgt i, src i+1)', 'i: <a_T> -> U', 'U = <a_T>+E_(i+1)  [new: line]', 'i+1: [new: <a_T> -> U\']'),
      ('stage-i total (reverse)', 'formed at <= P_j^perp', '', ''),
      ('stage-i copy (reverse)', 'climbs 8 to U (paid)', 'U = <a_T>+E_(i+1) for EVERY T in old coords', ''),
      ('stage-(i+1) total', '', 'formed at P_j\' (new)', 'copy drops 8 to E_(i+1) (new)')])]
    for title, rs in rows:
        print('   ', title)
        for r in rs: print('      %-26s | in window i: %-24s | at the boundary: %-44s | %s' % r)
    # ---------------- (3) port dependence and the old-identification pattern ----------------
    L, G = cf.e8_lines()
    B0 = [[(1 + G[i][j]) & 1 for j in range(120)] for i in range(120)]
    U = cf.f2_labels(B0)
    u = [sum(U[S][c] << c for c in range(len(U[S]))) for S in range(120)]
    h = len(U[0]); assert h == 9
    live = list(range(8, 120))
    assert all(dot(u[S], u[S]) == 1 for S in live)
    # t_T(a_T^perp) = Enew is a statement about the PRE-image: two ports with a_T != a_T' have different hyperplanes
    # mapped onto the same Enew, so t_T != t_T' (both injective).  Count the distinct preimages:
    hyp = {tuple(sorted(perp([u[S]], 9))) for S in live}
    canon = set()
    for S in live:
        b = perp([u[S]], 9)
        canon.add(frozenset(x for x in range(1, 512) if dot(x, u[S]) == 0))
    print('\n(3) distinct labels on P112: %d; distinct hyperplanes a_T^perp (= preimages of the common E_(i+1)): %d'
          % (len(set(u[S] for S in live)), len(canon)))
    print('    so the 112 transitions t_T are pairwise different maps: a register holding content of two ports cannot be')
    print('    re-identified into window i+1 by one map.')
    # knows in OLD identification at stage i+1: need a frame G with  (target side) and (source side) both satisfied.
    # fwd->rev: target = P-bank S: U <= G <= U+B ; source = 0-bank T: a_T^perp <= G <= a_T^perp + B.
    #   U <= a_T^perp + B  iff  U cap (a_T^perp+B) = U  iff a_T in a_T^perp  (never: N(a,a) = 1).
    # rev->fwd: target = 0-bank S: a_S^perp <= G <= a_S^perp + B ; source = P-bank T: U <= G: same obstruction with S.
    known_old = sum(1 for S in live for T in live if dot(u[T], u[T]) == 0)
    known_old2 = sum(1 for S in live for T in live if dot(u[S], u[S]) == 0)
    known_new = sum(1 for S in live for T in live if S != T and dot(u[S], u[T]) == 0)
    print('    pairs (S,T) with a common legal frame at stage i+1 in the OLD identification: fwd->rev %d, rev->fwd %d (of %d)'
          % (known_old, known_old2, 112 * 112))
    print('    same count in the NEW identification (the stage-(i+1) side information, off-diagonal): %d' % known_new)
    # ---------------- (4) distances on the complex P112 frame decoder ----------------
    print('\n(4) hand-over routes for the 8 totals A_a of the complex P112 frame decoder (served = read set, as cxp.py)')
    print('    a | w_in=dim P | R2 route: E_i -> U | mid-climb read at P^perp -> U | rev->fwd drop U -> common cut | route N: E_i -> E_(i+1) | star copy (new)')
    tot = [0] * 5
    for a in range(8):
        A = [S for S in live if G[S][a] != 0]
        P = [u[S] for S in A]
        wP = dim(P)
        Pp = perp(P, 9)
        cut = perp([u[S] for S in A], 9)        # cap of u_S^perp over served S = (span served)^perp
        r2 = 9 - 0
        mid = 9 - dim(Pp)
        drop = 9 - dim(cut)
        routeN = 8
        star = wP
        vals = [r2, mid, drop, routeN, star]
        tot = [x + y for x, y in zip(tot, vals)]
        print('    %d | %d | %d | %d | %d | %d | %d' % (a, wP, *vals))
    print('    per call (8 totals): R2 route %d, mid-climb %d, rev->fwd %d, route N %d, star copies %d' % tuple(tot))

if __name__ == '__main__':
    main()
