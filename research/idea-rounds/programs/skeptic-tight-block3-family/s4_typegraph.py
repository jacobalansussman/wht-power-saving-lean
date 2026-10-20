# skeptic-1 (wht44/ideas), 16 Oct 2026. The "type graph" behind s1's depth numbers, for multiplicity-free irreducibles:
# P0-types are vertices; h joins tau1 -> tau2 when its (tau2 <- tau1) block is non-zero. In a multiplicity-free module
# rho(C_d) is the sum of Hom(tau1, tau2) over pairs joined by a walk of length <= d (length 0: tau1 = tau2), so the
# depth a module sees is the largest h-graph distance between two types that K's blocks join.
import numpy as np
np.seterr(all='ignore')
import importlib.util
spec = importlib.util.spec_from_file_location('s1', __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 's1_depth_cuts.py'))
s1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s1)

for j, lab in ((8, 'V56'), (11, 'V105a'), (10, 'V84'), (12, 'V105 (irr12)'), (13, 'V105 (irr13)')):
    n, h, k, k2, T = s1.load([j])
    if any(len(t['slots']) > 1 for t in T):
        print(lab, 'has repeated types, skipped'); continue
    nt = len(T); names = [f"{t['d']}{'abcdefg'[sum(1 for u in T[:i] if u['d'] == t['d'])]}" for i, t in enumerate(T)]
    def blk(X, a, b): return np.linalg.norm(X[np.ix_(T[a]['slots'][0], T[b]['slots'][0])])
    H = np.array([[blk(h, a, b) > 1e-9 for b in range(nt)] for a in range(nt)])
    Kb = np.array([[blk(k, a, b) > 1e-9 for b in range(nt)] for a in range(nt)])
    K2b = np.array([[blk(k2, a, b) > 1e-9 for b in range(nt)] for a in range(nt)])
    # distances in the h-graph (directed tau1 -> tau2 via H[tau2, tau1])
    INF = 99; Dm = np.full((nt, nt), INF)
    for a in range(nt): Dm[a, a] = 0
    for a in range(nt):
        for b in range(nt):
            if H[a, b]: Dm[a, b] = min(Dm[a, b], 1)
    for c in range(nt):
        for a in range(nt):
            for b in range(nt):
                Dm[a, b] = min(Dm[a, b], Dm[a, c] + Dm[c, b])
    print(f'== {lab}: types {names}')
    print('   h joins (row <- col):'); [print('     ', names[a], ''.join('x' if H[a, b] else '.' for b in range(nt))) for a in range(nt)]
    print('   K blocks (row <- col):'); [print('     ', names[a], ''.join('x' if Kb[a, b] else '.' for b in range(nt))) for a in range(nt)]
    need = max(Dm[a, b] for a in range(nt) for b in range(nt) if Kb[a, b])
    need2 = max(Dm[a, b] for a in range(nt) for b in range(nt) if K2b[a, b])
    far = [(names[b], names[a], int(Dm[a, b])) for a in range(nt) for b in range(nt) if Kb[a, b] and Dm[a, b] == need]
    print(f'   depth seen for K = {need} (pairs at that distance: {far}); for K2 = {need2}')
