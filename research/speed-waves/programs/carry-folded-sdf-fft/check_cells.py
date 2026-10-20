# check_cells.py (as-5): exhaustive-ish check of the fused CCU2C cells against the intended arithmetic (Python model of CCU2C)
import random, sys
sys.path.insert(0, '.')
from gen import G
g = G(8, 'sdfx')
def chain(kind, A, B, C, D, cin):
    init = g.cell_init(kind); c = cin; out = 0
    for i in range(16):
        idx = A[i] | (B[i] << 1) | (C[i] << 2) | (D[i] << 3)
        p = (init >> idx) & 1; gg = (init >> (A[i] | (B[i] << 1))) & 1
        out |= (p ^ c) << i
        c = c if p else gg
    return out - 65536 if out & 0x8000 else out
def bits(v): return [(v >> i) & 1 for i in range(16)]
def hb(v): return [((v >> 1) >> i) & 1 for i in range(16)]   # bits of v>>1 (arith)
bad = 0
for t in range(20000):
    m = random.randint(-32768, 32767); x = random.randint(-32768, 32767); ph = random.randint(0, 1); jf = random.randint(0, 1) & ph
    xh = x >> 1; x0 = x & 1
    o = chain('OUT', bits(m), hb(x), [ph]*16, [1]*16, x0 & ph)
    e = s16 = None
    exp_o = ((m + xh + x0) if ph else m); exp_o = ((exp_o + 32768) % 65536) - 32768
    w = chain('MEM', hb(x), [ph]*16, bits(m), [1]*16, ph)
    exp_w = ((m - xh) if ph else xh); exp_w = ((exp_w + 32768) % 65536) - 32768
    oj = chain('OUTJ', bits(m), hb(x), [ph]*16, [jf]*16, ph & (jf | x0))
    exp_oj = (m - xh if jf else m + xh + x0) if ph else m; exp_oj = ((exp_oj + 32768) % 65536) - 32768
    wj = chain('MEMJ', bits(m), [ph]*16, hb(x), [1 - jf]*16, ph & ((1 - jf) | x0))
    exp_wj = (m + xh + x0 if jf else m - xh) if ph else xh; exp_wj = ((exp_wj + 32768) % 65536) - 32768
    for nm, a, b in (('OUT', o, exp_o), ('MEM', w, exp_w), ('OUTJ', oj, exp_oj), ('MEMJ', wj, exp_wj)):
        if a != b:
            bad += 1
            if bad < 6: print(nm, m, x, ph, jf, a, b)
print('bad', bad)
# Nexus cells: CCU2 carry = P ? cin : LUT3(A,B,C) = INIT[7:0]
def chain_n(kind, A, B, C, D, cin):
    init = g.cell_init_n(kind); c = cin; out = 0
    for i in range(16):
        idx = A[i] | (B[i] << 1) | (C[i] << 2) | (D[i] << 3)
        p = (init >> idx) & 1; gg = (init >> (A[i] | (B[i] << 1) | (C[i] << 2))) & 1
        out |= (p ^ c) << i
        c = c if p else gg
    return out - 65536 if out & 0x8000 else out
badn = 0
for t in range(20000):
    m = random.randint(-32768, 32767); x = random.randint(-32768, 32767); ph = random.randint(0, 1); jf = random.randint(0, 1) & ph
    xh = x >> 1; x0 = x & 1
    w16 = lambda v: ((v + 32768) % 65536) - 32768
    r = [(chain_n('OUT', bits(m), hb(x), [ph]*16, [1]*16, x0 & ph), w16(m + xh + x0 if ph else m)),
         (chain_n('MEM', hb(x), [ph]*16, bits(m), [1]*16, ph), w16(m - xh if ph else xh)),
         (chain_n('OUTJN', bits(m), hb(x), [jf]*16, [ph]*16, ph & (jf | x0)), w16((m - xh if jf else m + xh + x0) if ph else m)),
         (chain_n('MEMJN', [ph]*16, hb(x), [1 - jf]*16, bits(m), ph & ((1 - jf) | x0)), w16((m + xh + x0 if jf else m - xh) if ph else xh))]
    for k, (a_, b_) in enumerate(r):
        if a_ != b_:
            badn += 1
            if badn < 5: print('nexus cell', k, m, x, ph, jf, a_, b_)
print('nexus bad', badn)
