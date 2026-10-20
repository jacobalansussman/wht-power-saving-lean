"""v-g-crossstage-model: what stagemap.py's 'computed' columns actually compute.
(a) its old-identification count tests dot(u_T,u_T)==0 only (no frames, no S); (b) in the stated frame model
    (target side contains U, source side inside a^perp + B, B a fresh 8-block) the obstruction holds for EVERY
    nonzero a, isotropic or not, so N(a,a)=1 is not the cause; (c) its 'mid-climb' and 'rev->fwd drop' columns are
    the same perp computed twice; R2 = 9 and route N = 8 are literals.   usage: python3 -I check_step2.py"""
import itertools
def rank(vs):
    b = {}
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in b: v ^= b[h]
            else: b[h] = v; break
    return len(b)
def dot(a, b): return bin(a & b).count('1') & 1
n = 9; Bbits = [1 << (n + j) for j in range(8)]          # B = 8 fresh coordinates above U = F2^9
U = [1 << j for j in range(n)]
def perp(a): 
    out = []
    for x in range(1, 1 << n):
        if dot(x, a) == 0 and rank(out + [x]) > len(out): out.append(x)
    return out
bad = 0; iso = 0
for a in range(1, 1 << n):
    if dot(a, a) == 0: iso += 1
    src = perp(a) + Bbits                                  # largest source frame a^perp + B
    if rank(src + U) == rank(src): bad += 1                # U inside a^perp + B ?
print('nonzero a: 511, isotropic a: %d; a with U inside a^perp+B: %d  -> obstruction holds for all a, N(a,a) irrelevant' % (iso, bad))
