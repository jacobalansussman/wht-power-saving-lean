# fc-gate route 1: lever-reprice's variants() (lever_price.py + mixed.py, Decimal prec 60) on design B and controls.
import sys, io, contextlib, json
sys.path.insert(0, '<work>/wht45/frontier/lever-reprice/tools')
with contextlib.redirect_stdout(io.StringIO()):
    from lever_price import variants
from mixed import root, add, sep, D
P = {1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920, 48: 137952}
W = 322560
out = {}
for name, heads, non in (('designB', 80640, 57312), ('pub3945', 57312, 80640)):
    v = variants(name, P, W, heads, non, quiet=True)
    out[name] = {k: ('%.12e' % v[k]['b'], v[k]['separated']) for k in ('uniform', 'half_cut', 'full_cut', 'unsep_cut')}
    out[name]['breakeven_full'] = v['breakeven_full']
    out[name]['best_d'] = (v['best_d']['d'], '%.12e' % v['best_d']['b'])
print(json.dumps(out, indent=1))
