# k8-twin: exact prices with m-anatomy/pricer.py (reproduces the hub's moment engine). Control = #1042.
import subprocess, json, sys
P = __import__('os').path.join(__import__('os').environ.get('PRICER_DIR', '.'), 'pricer.py')
def price(G, W=14336, D=5120):
    out = subprocess.run([sys.executable, '-I', P, '--profile', json.dumps({str(k): v for k, v in sorted(G.items())}), '--m', '73', '--W', str(W), '--deficit', str(D)], capture_output=True, text=True, check=True).stdout
    d = json.loads(out); return d['b_decimal'], d['kappa_decimal'], d['binding']
IDLE1042 = {1: 2249, 2: 1, 3: 400, 7: 3, 16: 47, 48: 4850}
K8 = {1:2412,2:418,3:432,4:142,5:82,8:8}; G14 = {1:2406,2:432,3:424,4:142,5:82,8:8}; G1314 = {1:2400,2:446,3:416,4:142,5:82,8:8}
def glob(local, idle):
    G = {r: 144*x for r, x in local.items()}
    for r, x in idle.items(): G[r] = G.get(r, 0) + x
    return G
def add(d, e):
    d = dict(d)
    for k, v in e.items(): d[k] = d.get(k, 0) + v
    return {k: v for k, v in d.items() if v}
if __name__ == '__main__':
    rows = {'control #1042 (k8 word, #1042 map)': glob(K8, IDLE1042),
            'g14 (580/581) what-if: +6 bins of 48': glob(G14, add(IDLE1042, {48: 6})),
            'g14 what-if: +18 bins of 16': glob(G14, add(IDLE1042, {16: 18})),
            'g14, #1042 map kept: +288 rank-1 holes': glob(G14, add(IDLE1042, {1: 288})),
            'g13+g14 what-if: +12 bins of 48': glob(G1314, add(IDLE1042, {48: 12})),
            'g13+g14 what-if: +36 bins of 16': glob(G1314, add(IDLE1042, {16: 36}))}
    for a in sys.argv[1:]:
        lab, loc, idle = a.split('|'); rows[lab] = glob({'K8': K8, 'G14': G14, 'G1314': G1314}[loc], {int(k): v for k, v in json.loads(idle).items()})
    for k, G in rows.items():
        b, kap, bind = price(G); print(f'{k:55s} b={b} kappa={kap} d_kappa_vs_1042={float(kap)-0.001672188119205421:+.4e} {bind}')
