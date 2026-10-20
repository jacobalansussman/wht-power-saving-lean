"""chkL.py: sanity check of one LB2 contradiction with plain Python (tables5.bin data, final zero dump)."""
import struct, sys
raw = open('tables5.bin', 'rb').read(); nf = struct.unpack_from('<i', raw, 0)[0]
mask = struct.unpack_from('<%dI' % nf, raw, 4); wt = struct.unpack_from('<%di' % nf, raw, 4 + 4 * nf)
vid = struct.unpack_from('<%di' % nf, raw, 4 + 8 * nf); o = 4 + 12 * nf
nv = struct.unpack_from('<i', raw, o)[0]; pk = struct.unpack_from('<%dI' % nv, raw, o + 4)
ks = [int(x) for x in sys.argv[1].split(',')]; F = int(sys.argv[2]); S = int(sys.argv[3])
Z = [0] + [int(h, 16) for h in sys.argv[4].split()]
w = [sum(pk[v] >> k & 1 for k in ks) for v in range(nv)]
s = [wt[q] - w[vid[q]] for q in range(nf)]
bits = lambda m: [i + 1 for i in range(31) if m >> i & 1]
D = mask[F]; print('F', F, 'weight', wt[F], 'odd arrays', w[vid[F]], '| S inside D_F:', S & ~D == 0, '|S| =', bin(S).count('1'), bits(S))
def match(rows, colsfor):
    mt = {}
    def aug(x, seen):
        for y in colsfor(x):
            if y in seen: continue
            seen.add(y)
            if y not in mt or aug(mt[y], seen): mt[y] = x; return True
        return False
    return sum(1 for x in rows if aug(x, set()))
nz = lambda x, y: not (Z[x] >> (y - 1) & 1)
# upper bound on rank E[S,S]: forms containing S, and term rank of S x S
ub_forms = min(s[q] for q in range(1, nf) if mask[q] & S == S)
ub_match = match(bits(S), lambda x: [y for y in bits(S) if nz(x, y)])
print('upper bounds on rank E[S,S]: min slack of a form containing S =', ub_forms, '; term rank of S x S =', ub_match)
# lower bound: rank E[D_F] = |D_F| - nu <= rank E[S,S] + term rank of (D_F^2 minus S^2)
Sset = set(bits(S))
nuL = match(bits(D), lambda x: [y for y in bits(D) if nz(x, y) and not (x in Sset and y in Sset)])
nu = w[vid[F]]
print('rank E[D_F] =', bin(D).count('1') - nu, '; term rank of D_F^2 minus S^2 =', nuL, '-> rank E[S,S] >=', bin(D).count('1') - nu - nuL)
