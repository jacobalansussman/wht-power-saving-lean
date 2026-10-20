# rl-plan: #1042 under the separated variants (same formula as rescue.py). Counts from rl-vprice/out/vprice.json.
import sys; sys.path.insert(0, __file__.rsplit('/', 1)[0])
from rescue import root, add, sep, D
CALLS = {1: 347328, 2: 60192, 3: 62208, 4: 20448, 5: 11808, 8: 1152}; W = 14336
NEITHER = {1: 901, 2: 1, 3: 16, 16: 29, 48: 2539}; NET = {1: 1348, 3: 384, 7: 3}; LIT48 = 2311; LIT16 = 18
HUB = D('1.6797549561255445e-3')
def run(name, ch, tg):
    b = root(ch, tg); print(f"{name:66s} b={b:.10e} vs hub {float((b-HUB)/HUB)*100:+.3f}% separated={sep(ch, tg)}")
base = dict(CALLS)
for r, n in NEITHER.items(): base = add(base, r, n)
for r, n in NET.items(): base = add(base, r, n)
run("#1042 uniform (expect 1.6749890159e-3)", add(add(base, 48, LIT48), 16, LIT16), {73: W})
# A: literal heads fully cut (48->target 25, 16->57), neither rank 48 split 24+24
ch = dict(CALLS)
for r, n in NEITHER.items():
    if r == 48: ch = add(ch, 24, 2 * n)
    else: ch = add(ch, r, n)
for r, n in NET.items(): ch = add(ch, r, n)
run("#1042 literal heads cut fully, rank-48 neither split 24+24", ch, {73: W - LIT48 - LIT16, 25: LIT48, 57: LIT16})
# B: literal rank-48 heads cut by 24 (target 49), rank-16 heads cut fully (57), all else paid as is
ch = add(base, 24, LIT48)
run("#1042 literal 48-heads cut by 24 (49), 16-heads cut (57), rest paid", ch, {73: W - LIT48 - LIT16, 49: LIT48, 57: LIT16})
