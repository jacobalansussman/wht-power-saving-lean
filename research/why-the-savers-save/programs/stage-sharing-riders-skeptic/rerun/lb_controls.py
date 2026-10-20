"""r-riders / lb_controls.py -- controls for lb.py's bound (does it move when the target moves?).
  K  : Lemma A's own case, 6 live arrays each owing the kernel (expect 18 = 3W).
  S  : the stage-i end state alone (X_T = x_T@F, Y_T = (y_T + x_T)@psi_T), expect 15 if the scatter is forced.
  P  : planted cheap target, the data climbs only (X_T = x_T@F, Y_T = y_T@psi_T): expect 12.
  Bd : version B WITHOUT the cleanup (X_T keeps v-s-wide's rider: X_1 = (x_1+y_4)@F, X_2 = (x_2+y_1)@F,
       X_4 = (x_4+y_2)@F), Y_T and Z as in B: the 'what-if no cleanup' number.
  Zd : riders' best case: only the dirty X_T and Z = T2 owed (no scatter, Y_T unconstrained).
usage: python3 -B lb_controls.py K|S|P|Bd|Zd"""
import sys
sys.dont_write_bytecode = True
import toy, lb
F, PSI, E, phys = toy.F, toy.PSI, toy.E, toy.phys
v = sys.argv[1]
T = {}
RIDER = {1: ('y', 4), 2: ('y', 1), 4: ('y', 2)}
if v == 'K':
    for inp in toy.INPUTS: T['%s%d' % inp] = {(inp, F): 1}
else:
    for a in E:
        x, y = ('x', a), ('y', a)
        T['X%d' % a] = {(x, phys(x, F)): 1}
        if v in ('Bd', 'Zd'): T['X%d' % a][(RIDER[a], phys(RIDER[a], F))] = 1
        if v in ('S', 'Bd'): T['Y%d' % a] = {(y, phys(y, PSI[a])): 1, (x, phys(x, PSI[a])): 1}
        if v == 'P': T['Y%d' % a] = {(y, phys(y, PSI[a])): 1}
    if v in ('Bd', 'Zd'):
        T['Z'] = {}
        for a in E:
            T['Z'][(('y', a), phys(('y', a), F))] = 1; T['Z'][(('x', a), phys(('x', a), F))] = 1
lb.bound(T, v)
