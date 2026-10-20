"""k-combine: the exact packed mirror (m-gen's zmirror.suite, mutable kinds) of the scalar kernel checks of a
mutable-source certificate at Shea's-rule SPar, the smallest E7 digit width, and the sw-1 margin.
usage: spar_search.py <cert.gcert.json.gz> d,sw,ux,us,uy <out.json>"""
import sys, os, json, gzip
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zmirror as zm
c0 = json.loads(gzip.open(sys.argv[1], 'rt').read())
sp = tuple(map(int, sys.argv[2].split(',')))
c = zm.lean_data(c0); c['ext'] = c0['ext']
gl = [r for r, _ in c0['ext']]
batches = [gl[i:i + 80] for i in range(0, len(gl), 80)]
d, sw, ux, us, uy = sp
esw = 8
while not all(zm.evaluate(zm.eCheck, c, (d, esw, ux, us, uy), b, True)['value'] for b in batches):
    esw += 1; assert esw < 64
ep = (d, esw, ux, us, uy)
main = zm.suite(c, sp, 4, ep, 80, True)
marg = zm.suite(c, (d, sw - 1, ux, us, uy), 4, (d, esw - 1, ux, us, uy), 80, True)
res = dict(spar=sp, epar=ep, e_batches=[len(b) for b in batches], main=main, sw_minus_1=marg,
           main_all_true=main['all_true'], margin_scal_false=[x['value'] for x in marg['scal']],
           margin_echk=[x['value'] for x in marg['echk']])
json.dump(res, open(sys.argv[3], 'w'), indent=1)
print(json.dumps(dict(spar=sp, epar=ep, all_true=main['all_true'], scal=[x['value'] for x in main['scal']],
                      ychk=main['ychk']['value'], echk=[x['value'] for x in main['echk']],
                      margin_scal=res['margin_scal_false'], margin_echk=res['margin_echk'])))
