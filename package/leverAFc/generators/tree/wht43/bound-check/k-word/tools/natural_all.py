"""Natural-exit tier-2 price of every distinct hub word (lower_exit + natural exits + ceil/cyclic price)."""
import sys, os, json, gzip, importlib.util
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lower_exit as LE
from exit_price import ceil_unit
from route_price import price, unit
gxp, wdir, out = sys.argv[1], sys.argv[2], sys.argv[3]
spec = importlib.util.spec_from_file_location('gx', gxp); gx = importlib.util.module_from_spec(spec); spec.loader.exec_module(gx)
idx = json.load(open(os.path.join(wdir, 'INDEX.json')))
res = {}
for k, s in idx.items():
    if (s['v'], s['h'], s['cst']) != (112, 9, 64): continue
    c0 = LE.read(os.path.join(wdir, k + '.json.gz'))
    name = k + ' ' + s['origins'][0].split('research/', 1)[-1][:90]
    try:
        cert, ok, total = LE.lower(c0, gx)
    except AssertionError as e:
        print('LOWERING FAILED', name, repr(e)[:120], flush=True); res[name] = dict(error=repr(e)[:300]); continue
    v, R = cert['v'], cert['R']; n = 2 * v + R
    fr = [tuple(f) for f in cert['frames']]; dim = [len(f) for f in fr]
    Bcore = cert['B'][:-n]
    cur = list(cert['start'])
    for g in cert['A'] + Bcore:
        for r in gx.gate_regs(g)[0]: cur[r] = g[1]
    final = list(cert['final'])
    for r in range(2 * v, n): final[r] = cur[r]
    H = Counter(); cur = list(cert['start']); bad = 0
    def climb(r, f):
        global bad
        if cur[r] != f:
            if not (dim[f] > dim[cur[r]] and gx.inside(fr[cur[r]], fr[f])): bad += 1
            H[dim[f] - dim[cur[r]]] += 1; cur[r] = f
    for g in cert['A']:
        for r in gx.gate_regs(g)[0]: climb(r, g[1])
    for _, r, f in cert['ret']: H[dim[f]] += 1
    for g in Bcore:
        for r in gx.gate_regs(g)[0]: climb(r, g[1])
    for r in range(n): climb(r, final[r])
    resid = [dim[final[r]] - dim[cert['start'][r]] for r in range(2 * v, n)]
    N = sum(a * b for a, b in H.items())
    if bad or N - sum(resid) != 1856 or min(resid) < 1:
        print('SKIP', name, 'bad', bad, 'N-res', N - sum(resid), flush=True); res[name] = dict(error='identity'); continue
    u = ceil_unit(dict(H), resid, copies=False); pc = price(u)
    u2 = unit(dict(H), sorted(resid, reverse=True), copies=False); py = price(u2)
    res[name] = dict(hist=dict(sorted(H.items())), residual_sum=sum(resid), residuals=dict(sorted(Counter(resid).items())),
                     ceil=dict(W=u['W'], banks=u['banks'], pad=u['pad'], **pc), cyclic=dict(W=u2['W'], **py),
                     hub_profile_lowering_matches=ok)
    print(f"{pc['dft']} cyc {py['dft']} S={sum(resid)} B={u['banks']} P={u['pad']} {name}", flush=True)
json.dump(res, open(out, 'w'), indent=1)
