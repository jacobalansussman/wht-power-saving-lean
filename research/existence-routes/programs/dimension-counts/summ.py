import json, sys
for fn in sys.argv[1:]:
    for l in open(fn):
        o = json.loads(l)
        if 'error' in o: print(o['label'], o['error']); continue
        g = o.get('gap') or [0, 0]
        print(f"{o['label']:34s} P={o['P']:7d} Eeff={o['Eeff']:6d} naive={o['naive_eff']:7d} r={o.get('r',-1):6d} dom={str(o.get('dominant')):5s} "
              f"fib={o.get('fibre',-1):7d} orb={o['orbit']:6d} ok={str(o.get('orbit_ok')):5s} dimR={o['dimR']}/{o['Eeff']//o['W']} gap={g[0]:.0e}/{g[1]:.0e} {o['sec']}s")
