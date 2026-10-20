# write the pi classes (mod S3 bit conjugation) whose cut-LP value is <= 11 (others need >= 12 calls) [uses cutip]
import itertools, sys
sys.path.insert(0, '.')
import cutip
N = 8
bitperms = list(itertools.permutations(range(3)))
def act(g, S): return sum(((S >> k) & 1) << g[k] for k in range(3))
seen = set(); reps = []
for pi in itertools.permutations(range(N)):
    if pi in seen: continue
    orb = set()
    for g in bitperms:
        gi = [act(g, S) for S in range(N)]; ginv = [0]*N
        for S in range(N): ginv[gi[S]] = S
        orb.add(tuple(gi[pi[ginv[c]]] for c in range(N)))
    seen |= orb; reps.append(pi)
out = []
for pi in reps:
    v, r = cutip.lp(pi)
    if v <= 11 + 1e-7: out.append((v, pi))
out.sort()
with open('pis_b3_lp11.txt', 'w') as f:
    for v, pi in out: f.write(' '.join(map(str, pi)) + f'   # LP {v:.6f}\n')
print(len(reps), 'classes;', len(out), 'with LP <= 11 written to pis_b3_lp11.txt')
