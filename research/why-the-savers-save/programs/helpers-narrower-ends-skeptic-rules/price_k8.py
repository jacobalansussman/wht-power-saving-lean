import subprocess, json, sys
P = __import__('os').path.join(__import__('os').environ.get('PRICER_DIR', '.'), 'pricer.py')
def price(G, W=14336, D=5120):
    out = subprocess.run([sys.executable, '-I', P, '--profile', json.dumps({str(k): v for k, v in sorted(G.items())}), '--m', '73', '--W', str(W), '--deficit', str(D)], capture_output=True, text=True, check=True).stdout
    d = json.loads(out); return d['b_decimal'], d['kappa_decimal'], d['binding']
F1042 = {1:349577,2:60193,3:62608,4:20448,5:11808,7:3,8:1152,16:47,48:4850}
K8 = {1:2412,2:418,3:432,4:142,5:82,8:8}
idle = dict(F1042)
for r,x in K8.items(): idle[r] -= 144*x
idle = {r:x for r,x in idle.items() if x}
print('idle of #1042:', idle)
def glob(local, extra):
    G = {r:144*x for r,x in local.items()}
    for r,x in idle.items(): G[r] = G.get(r,0)+x
    for r,x in extra.items(): G[r] = G.get(r,0)+x
    return G
G13 = {1:2406,2:432,3:424,4:142,5:82,8:8}
G1314 = {1:2400,2:446,3:416,4:142,5:82,8:8}
rows = {
 'control #1042': price(F1042),
 'g13 (=#1090 pair) on #1042 map, freed 288 as 6 bins of 48 (what-if)': price(glob(G13, {48:6})),
 'g13+g14 on #1042 map, freed 576 as 12 bins of 48 (what-if)': price(glob(G1314, {48:12})),
 'g13+g14, #1042 map kept, 576 rank-1 holes': price(glob(G1314, {1:576})),
}
for k,v in rows.items(): print(f'{k:75s}', *v)
