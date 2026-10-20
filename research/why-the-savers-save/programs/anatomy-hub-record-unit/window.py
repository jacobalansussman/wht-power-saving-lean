"""an-hub: cost share by child rank for T74, and the grouping ceiling of the window architecture
(local children at most the window rank 8 or 9, idle at most R = 51). Stdlib only."""
import importlib.util, sys, os
spec = importlib.util.spec_from_file_location('an', os.path.join(os.path.dirname(__file__), 'anatomy.py'))
an = importlib.util.module_from_spec(spec); spec.loader.exec_module(an)
m, W = 73, 66304
t74 = {1: 1559100, 2: 339000, 3: 277392, 4: 84912, 5: 54612, 8: 5328, 48: 22748}
a = an.saving(t74, m, W)
tot = sum(n * (r / m) ** (1 - a) for r, n in t74.items())
print('T74 a =', a)
for r, n in sorted(t74.items()):
    print(f'  rank {r:2d}: mass share {r*n/4816512:.3f}  cost share {n*(r/m)**(1-a)/tot:.3f}')
local = {k: v for k, v in t74.items() if k < 48}
lm = sum(r * n for r, n in local.items()); im = 48 * 22748
for L in (8, 9):
    for I in (48, 51):
        p = {L: lm / L}; p[I] = p.get(I, 0) + im / I
        print(f'local mass all at rank {L}, idle at {I}: a = {an.saving(p, m, W):.6e}')
# centre path at rank 8 is part of local; keep it, only regroup the rest
