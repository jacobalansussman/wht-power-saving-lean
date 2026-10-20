import sys
from fractions import Fraction as F
sys.path.insert(0, __import__('os').path.dirname(__file__))
from wordlib import *
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
prog = program(w); chain, bad = frame_replay(w, prog)
rows = [{k: F(1)} for k in range(n)]; readval = {}
for T, (ph, e) in enumerate(prog):
    if e[0] == 'add': readval[T] = dict(rows[e[2]])
    apply(e, rows)
nm = lambda k: ('x%d' % k if k < v else 'y%d' % k if k < 2*v else 'h%d' % k)
def fmt(d): return ' '.join('%s%s' % (str(x) if x != 1 else '+', nm(k)) for k, x in sorted(d.items()))
for q in map(int, sys.argv[2:]):
    print('REG', nm(q), 'chain', [(T, len(U), U) for T, U in chain[q]])
    for T, (ph, e) in enumerate(prog):
        if q in regs(e) and ph != 'compensation':
            if e[0] == 'add':
                other = e[2] if e[1] == q else e[1]
                d = 'READS %s' % nm(e[2]) if e[1] == q else 'read by %s' % nm(e[1])
                print('  %6d %-8s %-14s q=%-5s fr=%s  val[b]=%s' % (T, ph, d, F(e[3], e[4]), (len(basis(e[-1])) if e[-1] is not None else '-'), fmt(readval[T])[:160]))
            else: print('  ', T, ph, e)
    nc = sum(1 for ph, e in prog if ph == 'compensation' and q in regs(e)); print('  compensation gates', nc)
