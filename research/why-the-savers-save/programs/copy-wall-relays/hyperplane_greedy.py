"""hyperplane_greedy.py (d-relay) -- decoder-level (NO closure) witness search, restricted family: rank-8 row space
V_f = {r_q : q.f = 0} inside the row space of B0 = U U^T (r_q[T] = q.u_T), any total order of the 120 targets.
Target S, with T 'later' than S, needs q with q.f = 0, q.u_S = 1, q.u_T = 0 for every later T not orthogonal to S:
possible iff u_S is not in span({f} u {u_T : T later, not orthogonal to S}).  Backward greedy (place latest first),
random restarts; ports with u_S.f = 0 are served statically by q = u_S and are placed earliest.  A success = a rank-8 relay decoder at the decoder level (closure ignored).  Standard library.
usage: python3 -I hyperplane_greedy.py [P120|P111ours] [tries per f]"""
import sys, random
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-floor-checker'))
import copyfloor as cf
L, G = cf.e8_lines()
B0 = [[(1 + G[i][j]) % 2 for j in range(120)] for i in range(120)]
U = cf.f2_labels(B0)
UI = [sum(U[i][c] << c for c in range(len(U[i]))) for i in range(120)]
inst = sys.argv[1] if len(sys.argv) > 1 else 'P120'
tries = int(sys.argv[2]) if len(sys.argv) > 2 else 20
live = list(range(120)) if inst == 'P120' else [i for i in range(120) if i not in (8, 9, 10, 12, 16, 17, 24, 32, 64)]
assert all(((bin(UI[a] & UI[b]).count('1')) & 1) == B0[a][b] for a in range(120) for b in range(120))
def inspan(x, vecs):
    basis = []
    for y in vecs:
        for b in basis: y = min(y, y ^ b)
        if y: basis.append(y); basis.sort(reverse=True)
    for b in basis: x = min(x, x ^ b)
    return x == 0
rng = random.Random(1)
best = (len(live) + 1, None)
for f in range(1, 1 << 9):
    for t in range(tries):
        placed = []; rest = list(live); rng.shuffle(rest)
        while rest:
            ok = [S for S in rest if not inspan(UI[S], [f] + [UI[T] for T in placed if G[S][T] != 0])]
            if not ok: break
            ok1 = [S for S in ok if bin(UI[S] & f).count('1') & 1]   # u_S.f = 1 first (latest); u_S.f = 0 ports are
            S = rng.choice(ok1 or ok); placed.append(S); rest.remove(S)  # served statically and go earliest
        if len(rest) < best[0]: best = (len(rest), f)
        if not rest:
            print('FOUND f=%d order(latest first)=%s' % (f, placed)); sys.exit(0)
print('%s: no rank-8 order found in the V_f family (%d tries per f); fewest stuck targets %d (f=%d)' % (inst, tries, best[0], best[1]))
