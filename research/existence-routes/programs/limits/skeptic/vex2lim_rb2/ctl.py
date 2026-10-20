"""control: my letters (myxw.py, written from the spec) vs xw2lib (reference) on every letter and some words."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, myxw
from xw2lib_ref import Game
names = ['P','C01','C10','X1','X2','X3','Z1','Z2','Z3','S1','S2','S3','CZ','CS','CSi','Om','Neg','A','Ai','a0','a1',
         'a0i','a1i','E01','E10','B','Bi','Pr0','Pr1','lX1','lX2','lZ1','lZ2','lS1','lS2','lCZ','lCS','lC01','fX1','fX2',
         'fZ1','fZ2','fS1','fS2','fCZ','fCS','A.CS.Ai','Ai.lS1.A','B.A','Pr0.X1','X1.Pr0','A.B','E01.Pr1.lCS','Bi.fCS.B']
G = Game([2, 3, 4, 5]); bad = 0
for nm in names:
    for f in (2, 3, 4, 5):
        t, p = myxw.word(nm, f); g = G.mono(nm)[f]
        if not (np.array_equal(t, g.sig) and np.array_equal(p % 4, g.ph % 4)):
            bad += 1; print('MISMATCH', nm, f)
# group-theoretic sanity of my code alone: A.Ai = 1, B.Bi = 1, Pr0^4 = 1, P.P = 1
for w in ['A.Ai', 'Ai.A', 'B.Bi', 'Pr0.Pr0.Pr0.Pr0', 'P.P', 'a0.a0i', 'CS.CSi']:
    for f in (3, 5):
        t, p = myxw.word(w, f)
        assert np.array_equal(t, np.arange(4 ** f)) and not p.any(), (w, f)
# Pr0 conjugated by X1: X1.Pr0 vs Pr0.X1 differ by (-1)^f? (claimed in PROOF section 3 E1)
for f in (3, 4, 5, 6):
    t1, p1 = myxw.word('X1.Pr0', f); t2, p2 = myxw.word('Pr0.X1', f)
    d = (p1 - p2) % 4
    print('f=%d: X1.Pr0 vs Pr0.X1: same targets %s, phase difference values %s' % (f, np.array_equal(t1, t2), np.unique(d)))
print('ctl: %d letters/words x 4 f compared with xw2lib: %d mismatches' % (len(names), bad))
