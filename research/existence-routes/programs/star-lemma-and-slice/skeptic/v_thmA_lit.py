# breakdown of the literal program's failures by linear part lambda (diagnostic for the '448' statement)
exec(open('v_thmA.py').read().split('# --- the whole character group')[0].replace('for f in range(1, 5):', 'for f in []:'))
gens = [tuple(bin(z & y).count('1') % 2 for z in range(1, 8)) for y in range(8)]
S = {tuple([0] * 7)}
while True:
    S2 = S | {tuple((a + b) % 4 for a, b in zip(s, g)) for s in S for g in gens}
    if S2 == S: break
    S = S2
from collections import Counter
tot = Counter(); bad = Counter()
for ch in S:
    psi = {z: ch[z - 1] for z in range(1, 8)}
    lam = tuple(psi[1 << j] % 2 for j in range(m))
    target = 1j ** (sum(psi[1 << j] for j in range(m)) % 4)
    e = np.abs(run(psi, False) - target * np.eye(W)).max()
    tot[lam] += 1; bad[lam] += e > 1e-8
for lam in sorted(tot): print("lambda", lam, "characters", tot[lam], "literal wrong", bad[lam])
