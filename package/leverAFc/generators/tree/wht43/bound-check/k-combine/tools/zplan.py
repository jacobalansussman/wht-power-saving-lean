"""k-lower: the batches of the zero checks (GS.zCheck, Work/KLow/ZCheck.lean) that give the two exit conditions,
and an EXACT mirror of the kernel term (no source-class test: y fillers allowed, as Lean's zCheck) (built on m-gen's zmirror.Run, itself a function-by-function mirror of
Sussman's / Bridge3's gPre, gInitL, gGates, gScat; classic kinds okAK / okBK).

Batches (sources gl, rows):
  bit batches   gl = [e in el | testBit(idxOf e, b) == val], b < nb (el.length <= 2^nb), val in {F, T};
                rows = every slot register.  => ExitCols, and ExitRows at exit columns.
  chunks        gl = consecutive pieces (<= v) of the non-y, non-exit registers; rows = el.  => ExitRows elsewhere.
usage: zplan.py <cert.gcert.json.gz> <out.json> [sw]
"""
import sys, json, gzip, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zmirror as zm
from zmirror import Run, No, guard, okAK, okBK

def bitB(el, b, val):
    return [e for e in el if ((el.index(e) >> b) & 1) == val]

def plan(c):
    v, R = c['v'], c['R']
    n = 2 * v + R
    el = sorted(c['exits'])
    nb = max(1, (len(el) - 1).bit_length())
    slots = list(range(2 * v, n))
    bits = [(b, val, bitB(el, b, val)) for b in range(nb) for val in (0, 1)]
    els = set(el)
    nel = [r for r in range(n) if (r < v or 2 * v <= r) and r not in els]
    k = -(-len(nel) // v)
    size = -(-len(nel) // k)
    chunks = [nel[i:i + size] for i in range(0, len(nel), size)]
    return dict(el=el, nb=nb, slots=slots, bits=bits, nel=nel, chunks=chunks)

def zCheck(c, p, gl, rows):
    Rn = Run(c, p, 0, len(gl))
    lo = 0
    for r in gl:
        guard(lo <= r and r + 1 <= Rn.nr, "incK")
        lo = r + 1
    guard(Rn.us + 1 <= 2 ** (Rn.sw - 1 if Rn.sw >= 1 else 0), "us >= 2^(sw-1)")
    guard(Rn.us + 1 <= 2147483648, "us >= 2^31")
    Rn.gInitL(Rn.us, gl)
    Rn.items(okAK, okBK)
    g = set(gl)
    for r in rows:
        if r in g:
            continue
        guard(Rn.key_ok(r), "zFin: register outside the trie")
        x = Rn.tr.get(r, 0)
        if ((x % Rn.m2) >> 64) != ((x >> Rn.K) >> 64):
            raise No("zFin: register %d holds something" % r)
    return Rn.adds

def run(c0, sw, d=10, ux=1, us=8, uy=8):
    c = zm.lean_data(c0)
    pl = plan(c0)
    p = (d, sw, ux, us, uy)
    res = []
    for b, val, gl in pl['bits']:
        res.append(dict(kind='bit', b=b, val=val, size=len(gl), **zm.evaluate(zCheck, c, p, gl, pl['slots'])))
    for i, gl in enumerate(pl['chunks']):
        res.append(dict(kind='chunk', i=i, size=len(gl), **zm.evaluate(zCheck, c, p, gl, pl['el'])))
    return pl, res

if __name__ == '__main__':
    c0 = json.loads(gzip.open(sys.argv[1], 'rt').read())
    out = sys.argv[2]
    sws = [int(sys.argv[3])] if len(sys.argv) > 3 else [8, 12, 16, 20, 24, 28, 32]
    allres = {}
    for sw in sws:
        pl, res = run(c0, sw)
        ok = all(r['value'] for r in res)
        allres[sw] = dict(all_true=ok, res=res)
        print(sw, ok, [ (r['kind'], r.get('b', r.get('i')), r['value']) for r in res if not r['value']][:6], flush=True)
        if ok and len(sws) > 1:
            # margin: sw - 1
            pl2, res2 = run(c0, sw - 1)
            allres['%d-1' % sw] = dict(all_true=all(r['value'] for r in res2))
            print('margin', sw - 1, allres['%d-1' % sw]['all_true'])
            break
    json.dump(dict(plan=dict(el=pl['el'], nb=pl['nb'], chunks=pl['chunks'], bits=[[b, v, g] for b, v, g in pl['bits']]),
                   results=allres), open(out, 'w'))
