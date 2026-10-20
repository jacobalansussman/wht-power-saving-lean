"""v-s-wide / rider_check.py -- escape (b) inside s-wide's own toy (handover_toy.py geometry), exact frames.
Frames are q: F2^3 -> Z4 (Fourier side), dir(z) adds chi_z EXACTLY (no shift quotient used anywhere below).
(1) every free 'ride': an input term at one of its free frames joins a stage-i DATA climb (X_T or Y_T, paid anyway)
    at its start or mid point and is carried to the climb's end frame;
(2) twin paths: two different data climbs that pass the same frame and then take the same remaining moves;
(3) a concrete rider assignment that makes X_1 + X_2 + X_4 equal s-wide's coded total T2 = sum_T (y_T + x_T)@F as a
    formal sum, with 0 paid moves beyond the data climbs (s-wide's toy reports 3 for T2 in both N and H), and the
    pollution it leaves in each live X_T (the cost s-wide's toy does not model).   usage: python3 -I rider_check.py"""
import sys
sys.dont_write_bytecode = True
Y = range(8)
def chi(z): return tuple(bin(y & z).count('1') & 1 for y in Y)
def add(q, r): return tuple((a + b) & 3 for a, b in zip(q, r))
Z0 = tuple(0 for _ in Y)
E = [1, 2, 4]
ports = {a: [d for d in E if d != a] for a in E}            # same climb order as handover_toy.py
F = add(add(chi(1), chi(2)), chi(4))
def climb(start, dirs):
    fr = [start]
    for z in dirs: fr.append(add(fr[-1], chi(z)))
    return fr
X = {a: climb(chi(a), ports[a]) for a in E}                 # chi_a -> chi_a+chi_d1 -> F
Yc = {a: climb(Z0, ports[a]) for a in E}                    # 0 -> chi_d1 -> psi
psi = {a: Yc[a][2] for a in E}
assert all(X[a][2] == F for a in E)
free = {}
for a in E:
    free[('x', a)] = {'start': X[a][0], 'mid': X[a][1], 'end(F)': X[a][2]}
    free[('y', a)] = {'start(0)': Yc[a][0], 'mid': Yc[a][1], 'end(psi)': Yc[a][2]}
name = {F: 'F'}; name.update({psi[a]: 'psi_%d' % a for a in E})
print('(1) free rides on stage-i data climbs that land on a frame some output needs (F or a psi):')
rides = []
for (kind, s), frs in free.items():
    for lab, f in frs.items():
        for host, path in [('X', X), ('Y', Yc)]:
            for a in E:
                for j, jp in [(0, 'start'), (1, 'mid')]:
                    if path[a][j] == f:
                        end = path[a][2]
                        if end in name and not (host == ('X' if kind == 'x' else 'Y') and a == s):
                            rides.append((kind, s, lab, host, a, jp, name[end]))
                            print('   %s_%d at its %-8s rides %s_%d from its %-5s -> %s_%d @ %s' % (kind, s, lab, host, a, jp, kind, s, name[end]))
print('\n(2) twin paths (different data climbs, same frame, same remaining moves):')
paths = [('X_%d' % a, X[a], ports[a]) for a in E] + [('Y_%d' % a, Yc[a], ports[a]) for a in E]
for i in range(len(paths)):
    for j in range(i + 1, len(paths)):
        (n1, f1, d1), (n2, f2, d2) = paths[i], paths[j]
        for k in range(2):
            if f1[k] == f2[k] and d1[k:] == d2[k:]:
                print('   %s and %s: same frame after %d move(s), same remaining moves %s, same end frame %s'
                      % (n1, n2, k, d1[k:], name.get(f1[2], f1[2])))
print('\n(3) riders: y_4 into X_1 at start (from Y_4 mid), y_1 into X_2 at start (from Y_1 mid), y_2 into X_4 at mid')
print('    (from Y_2 end, before Y_2\'s stage-i scatter).  Formal sums after the X climbs (input -> exact frame):')
reg = {a: {('x', a): [X[a][0]]} for a in E}
def ride(a, inp, f, at):
    assert f == X[a][at], 'rider frame does not match host'
    reg[a].setdefault(inp, []).append(('at', at, f))
assert Yc[4][1] == X[1][0]; assert Yc[1][1] == X[2][0]; assert Yc[2][2] == X[4][1]
# carry terms exactly
final = {}
for a, (rider, f, at) in {1: (('y', 4), Yc[4][1], 0), 2: (('y', 1), Yc[1][1], 0), 4: (('y', 2), Yc[2][2], 1)}.items():
    terms = {('x', a): X[a][0], rider: f}
    for k, z in enumerate(ports[a]):
        for inp in terms:
            if inp == rider and k < at: continue
            terms[inp] = add(terms[inp], chi(z))
    final[a] = terms
    print('   X_%d = %s' % (a, ' + '.join('%s_%d@%s' % (i[0], i[1], name.get(g, g)) for i, g in sorted(terms.items()))))
tot = {}
for a in E:
    for inp, g in final[a].items(): tot[(inp, g)] = tot.get((inp, g), 0) + 1
T2 = {}
for a in E: T2[(('y', a), F)] = 1; T2[(('x', a), F)] = 1
print('   X_1 + X_2 + X_4 == T2 (coded total sum_T (y_T + x_T)@F, exact frames):', tot == T2)
print('   paid moves beyond the data climbs: 0 (s-wide toy: 3 for T2 with N and with H)')
print('   left to pay for (outside s-wide\'s toy): each live X_T now holds a rider y_S@F; the stage-i copy of X_T')
print('   (F -> psi_T) carries it into Y_T as junk unless cleaned.  Whether cleaning costs less than the')
print('   3 moves saved is the open question escape (b) poses; neither s-wide\'s search nor this script decides it.')
