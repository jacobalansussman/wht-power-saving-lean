#!/usr/bin/env python3
"""d-hpe-multigate STEP 3: price the rank-7 end words with the m-anatomy pricer (exact moment engine + outer recipe).
Global = calls x local + idle residual (hub rule).  Packing brackets: 'perfect' = the freed width re-packed into full
rank-48 bins (needs a literal owner map, not built here); 'holes' = the predecessor owner map kept, each narrowed helper
leaves a rank-1 hole in its bin (a legal but bad assignment).  Usage: python3 -I price_variants.py"""
import subprocess, json, sys
P = __import__('os').path.join(__import__('os').environ.get('PRICER_DIR', '.'), 'pricer.py')
def price(local, calls, idle, W, D):
    G = {r: calls*x for r, x in local.items()}
    for r, x in idle.items(): G[r] = G.get(r, 0)+x
    out = subprocess.run([sys.executable, '-I', P, '--profile', json.dumps({str(k): v for k, v in sorted(G.items())}),
                          '--m', '73', '--W', str(W), '--deficit', str(D)], capture_output=True, text=True, check=True).stdout
    d = json.loads(out); return dict(profile=G, b=d['b_decimal'], kappa=d['kappa_decimal'], binding=d['binding'])
U8 = {1: 2548, 2: 338, 3: 446, 4: 162, 5: 72, 8: 8}
G17 = {1: 2542, 2: 352, 3: 438, 4: 162, 5: 72, 8: 8}
G1718 = {1: 2536, 2: 366, 3: 430, 4: 162, 5: 72, 8: 8}
rows = {
 'control #942 (333 calls, T37 idle)': price(U8, 333, {48: 11086, 1: 6}, 33152, 11840),
 'PR1030 pair (group 17) in 72-call geometry, its map {48:2400}': price(G17, 72, {48: 2400}, 7168, 2560),
 'groups 17+18, 333 calls, perfect re-pack {48:11113,42:1}': price(G1718, 333, {48: 11113, 42: 1}, 33152, 11840),
 'groups 17+18, 333 calls, T37 map kept + 1332 holes': price(G1718, 333, {48: 11086, 1: 1338}, 33152, 11840),
 'groups 17+18, 72 calls, perfect re-pack {48:2403}': price(G1718, 72, {48: 2403}, 7168, 2560),
 'groups 17+18, 72 calls, PR1030 map kept + 144 holes': price(G1718, 72, {48: 2400, 1: 144}, 7168, 2560),
}
for k, r in rows.items(): print(f"{k:70s} b {r['b']}  kappa {r['kappa']}  {r['binding']}")
json.dump(rows, open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'price_variants.json'), 'w'), indent=1)
