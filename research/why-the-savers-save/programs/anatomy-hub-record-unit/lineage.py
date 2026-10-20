"""an-hub: D/E split across hub checkpoints #588..#1042 (+T74). W is implied by bill(a)=1 at the hub's complex saving.
Stdlib only."""
import json, math, sys
m = 73
rows = [json.loads(l) for l in open(sys.argv[1])]
rows.append({'name': 'T74 (PR1424 era)', 'complex': '0.0016796717955047376', 'kappa': '0.001676855229071321',
             'profile': {'1': 1559100, '2': 339000, '3': 277392, '4': 84912, '5': 54612, '8': 5328, '48': 22748}})
print(f"{'checkpoint':18s} {'kappa':>12s} {'a (hub)':>12s} {'implied W':>12s} {'W/896':>8s} {'mass':>10s} {'1-bill0':>10s} {'a rank-1':>10s} {'r_eff':>7s} {'idle48 share':>12s}")
for r in rows:
    p = {int(k): int(v) for k, v in r['profile'].items()}
    a = float(r['complex'])
    W = sum(n * (k / m) ** (1 - a) for k, n in p.items())
    mass = sum(k * n for k, n in p.items())
    x = 1 - mass / (m * W)
    reff = math.exp(sum(k * n * math.log(k) for k, n in p.items()) / mass)
    big = sum(k * n for k, n in p.items() if k >= 16) / mass
    print(f"{r['name']:18s} {float(r['kappa']):12.6e} {a:12.6e} {W:12.2f} {W/896:8.3f} {mass:10d} {x:10.4e} {-math.log(1-x)/math.log(m):10.4e} {reff:7.3f} {big:12.3f}")
