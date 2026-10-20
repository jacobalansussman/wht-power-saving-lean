# v_plain4.py (v-e3-slice6): validity calibration of star rows at block 4: the plain 1 in 4 (c = 1 on e_1..e_4) must
# satisfy every "point + its 3 lines inside an even plane support" row  2 c_x <= sum_L w_L, with w_L at its largest
# allowed value min(c_x + c_y, min over odd supports G containing L of (sum_G c - 1)).
import itertools
m = 4; pts = list(range(1, 16)); pc = lambda x: bin(x).count('1')
mons = [(j, k) for j in range(m) for k in range(j, m)]
monmask = [sum(1 << z for z in pts if (z >> j) & 1 and (z >> k) & 1) for (j, k) in mons]
odd, even = set(), set(); mask = 0; tr = 0
for g in range(1 << len(mons)):
    if g:
        b = (g & -g).bit_length() - 1; mask ^= monmask[b]; tr ^= (mons[b][0] == mons[b][1])
    if mask: (odd if tr else even).add(mask)
c = {z: (1 if pc(z) == 1 else 0) for z in pts}
def wL(L):
    Lm = sum(1 << z for z in L)
    return min([sum(c[z] for z in L)] + [sum(c[z] for z in pts if (G >> z) & 1) - 1 for G in odd if G & Lm == Lm])
worst = None; n = 0
for F in even:
    if pc(F) != 4: continue
    Fp = [z for z in pts if (F >> z) & 1]
    if any(a ^ b ^ cc not in Fp for a, b, cc in itertools.combinations(Fp, 3)): continue   # F must be a plane
    for x in Fp:
        lines = [(x, y) for y in Fp if y != x]
        sl = sum(wL(L) for L in lines) - 2 * c[x]; n += 1
        worst = sl if worst is None else min(worst, sl)
print("block 4: even plane supports used; rows", n, "; plain 1 in 4 min slack", worst, "(>= 0: the rows accept the plain schedule)")
