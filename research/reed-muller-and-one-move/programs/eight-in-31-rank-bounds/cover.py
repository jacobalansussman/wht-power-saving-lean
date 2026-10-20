"""cover.py: forced-full-rank test from matching equality (r3-classes).
For a form Q with nu(D_Q) = s_Q (term rank of E[D_Q,D_Q] avoiding derived zeros equals the required rank), every
minimum vertex cover (R', C') (rows R', cols C', |R'|+|C'| = s_Q) has E[D\\R', D\\C'] = 0, so rank = s_Q forces
E[R', D\\C'] to have full row rank |R'| and E[D\\R', C'] full column rank |C'|.  Any odd form G whose support holds
R' u (D\\C') then needs s_G >= |R'| (and likewise s_G >= |C'| for (D\\R') u C').  Here: one cover per form (from the
maximum matching, Koenig), checked against all forms.  Input: tables5.bin, targets, linkrec/starrec Z line (hex masks)."""
import struct, sys
raw = open(sys.argv[1], 'rb').read(); nf = struct.unpack_from('<i', raw, 0)[0]
mask = struct.unpack_from('<%dI' % nf, raw, 4); wt = struct.unpack_from('<%di' % nf, raw, 4 + 4 * nf)
vid = struct.unpack_from('<%di' % nf, raw, 4 + 8 * nf); o = 4 + 12 * nf
nv = struct.unpack_from('<i', raw, o)[0]; pk = struct.unpack_from('<%dI' % nv, raw, o + 4)
ks = [int(x) for x in sys.argv[2].split(',')]
Z = [0] + [int(h, 16) for h in sys.argv[3].split()]       # Z[x]: bit y-1 set <-> E_xy derived 0 (x,y in 1..31)
w = [sum(pk[v] >> k & 1 for k in ks) for v in range(nv)]
s = [wt[q] - w[vid[q]] for q in range(nf)]
def bits(m): return [i + 1 for i in range(31) if m >> i & 1]
def match(R, C):
    mt = {}
    def aug(x, seen):
        for y in C:
            if (Z[x] >> (y - 1) & 1) or y in seen: continue
            seen.add(y)
            if y not in mt or aug(mt[y], seen): mt[y] = x; return True
        return False
    n = 0
    for x in R:
        if aug(x, set()): n += 1
    return n, mt
def koenig(R, C, mt):
    # Koenig cover from a maximum matching: Zr = vertices reachable from unmatched rows by alternating paths
    rowm = {x: y for y, x in mt.items()}
    visR, visC = set(), set()
    stack = [x for x in R if x not in rowm]
    visR.update(stack)
    while stack:
        x = stack.pop()
        for y in C:
            if (Z[x] >> (y - 1) & 1) or y in visC: continue
            visC.add(y)
            if y in mt and mt[y] not in visR: visR.add(mt[y]); stack.append(mt[y])
    Rc = [x for x in R if x not in visR]; Cc = [y for y in C if y in visC]
    return Rc, Cc
# super-set minimum slack table: smallest s_G over forms G whose support contains a given set (brute force on demand)
supp = [(mask[q], s[q]) for q in range(1, nf)]
def minslack(S):
    b = 99
    for m_, sg in supp:
        if m_ & S == S and sg < b: b = sg
    return b
neq = 0; viol = 0; checked = 0
for q in range(1, nf):
    if s[q] <= 0 or w[vid[q]] == 0: continue
    D = bits(mask[q]); n, mt = match(D, D)
    if n < s[q]: print('MATCHING VIOLATION', q, n, s[q]); viol += 1; continue
    if n > s[q]: continue
    neq += 1
    Rc, Cc = koenig(D, D, mt); assert len(Rc) + len(Cc) == s[q]
    for (Rr, Cs, lab) in ((Rc, [y for y in D if y not in Cc], 'row'), ([x for x in D if x not in Rc], Cc, 'col')):
        need = len(Rr) if lab == 'row' else len(Cs)
        if need == 0: continue
        S = 0
        for x in set(Rr) | set(Cs): S |= 1 << (x - 1)
        checked += 1
        b = minslack(S)
        if b < need:
            print('FORCED-RANK VIOLATION form %d (s=%d) %s cover |R'"'"'|=%d |C'"'"'|=%d: set %s inside a form of slack %d < %d' % (q, s[q], lab, len(Rc), len(Cc), bits(S), b, need)); viol += 1
print('forms at matching equality', neq, 'rectangles checked', checked, 'violations', viol)
