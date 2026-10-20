"""vderived.py -- v-ex2-modp: are ex2-modp's DERIVED systems (sub_*: chart relaxation + t(x-L)-1; ext_*: chart
relaxation + t e - 1) built correctly from the raw exports?  Independent parse; compares polynomials mod each prime."""
import os, re, glob, ast, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import *
M = os.environ.get('MODP_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
def remap(P, src, dst):
    out = {}
    for k, c in P.items():
        e = [0] * len(dst)
        for i, x in enumerate(k):
            if x: e[dst.index(src[i])] += x
        out[tuple(e)] = c
    return out
bad = 0; n = 0
# sub systems
for subdir, decdir in [(f'{M}/data/sub_f27_s39', f'{M}/data/dec/f27_Fp2_m3_s39'),
                       (f'{M}/data/hard/sub_g27_Fany4_s17', f'{M}/data/hard/dec/g27_Fany4_s17'),
                       (f'{M}/data/hard/sub_g27_Fany4_s19', f'{M}/data/hard/dec/g27_Fany4_s19')]:
    for f in sorted(glob.glob(f'{subdir}/sub_*.p1.ms')):
        b = os.path.basename(f)[4:-6]; chart, x = b.rsplit('_', 1)
        log = open(f'{decdir}/full_c_{chart}.log').read()
        for j, p in ((1, P1), (2, P2)):
            sn, sp, SP = read_ms(f'{subdir}/sub_{b}.p{j}.ms')
            cn, cp, CP = read_ms(f'{decdir}/c_{chart}.p{j}.ms')
            L = ast.literal_eval(re.search(rf"p{j}: linear substitutions from the relaxation: (\{{.*\}})", log).group(1))[x]
            ok = sp == p == cp and sn == cn + ['t'] and len(SP) == len(CP) + 1
            ok = ok and all(remap(P, cn, sn) == Q for P, Q in zip(CP, SP[:-1]))
            # expected last row: t*x - t*L - 1
            exp = {}
            ti = sn.index('t'); xi = sn.index(x)
            k = [0] * len(sn); k[ti] = 1; k[xi] = 1; exp[tuple(k)] = 1
            for v, c in L.items():
                k = [0] * len(sn); k[ti] = 1
                if v != '1':
                    for w in v.split('*'): k[sn.index(w)] += 1
                exp[tuple(k)] = (exp.get(tuple(k), 0) - c) % p
            exp[tuple([0] * len(sn))] = p - 1
            exp = {k: v for k, v in exp.items() if v % p}
            ok = ok and SP[-1] == exp
            n += 1
            if not ok: bad += 1; print('BAD sub', subdir, b, j)
print('sub systems checked', n, 'bad', bad)
# ext systems
n = bad = 0
D = f'{M}/data/hard/dec/g27_Fany4_s17'
for f in sorted(glob.glob(f'{M}/data/hard/ext_s17/ext_*.p1.ms')):
    b = os.path.basename(f)[4:-6]; m = re.match(r'(.*)_a(\d)_(\d+)$', b); chart, arr, q = m.group(1), m.group(2), int(m.group(3))
    k = int(re.search(rf"p1 array {arr}:.*extras (\d+)", open(f'{D}/full_c_{chart}.log').read()).group(1))
    for j, p in ((1, P1), (2, P2)):
        en, ep, EP = read_ms(f'{M}/data/hard/ext_s17/ext_{b}.p{j}.ms')
        cn, cp, CP = read_ms(f'{D}/c_{chart}.p{j}.ms')
        fn, fp, FP = read_ms(f'{D}/full_c_{chart}.a{arr}.p{j}.ms')
        e = FP[len(FP) - k + q]
        ok = en == cn + ['t'] and len(EP) == len(CP) + 1 and all(remap(P, cn, en) == Q for P, Q in zip(CP, EP[:-1]))
        ok = ok and all(not any(x for i, x in enumerate(mon) if not fn[i].startswith('a')) for mon in e)
        te = remap(e, fn, en); ti = en.index('t')
        exp = {}
        for mon, c in te.items():
            mm = list(mon); mm[ti] += 1; exp[tuple(mm)] = c
        exp[tuple([0] * len(en))] = p - 1
        ok = ok and EP[-1] == exp
        n += 1
        if not ok: bad += 1; print('BAD ext', b, j)
    # and the full system's extras are exactly rows len-k.. : are there k ext files for this chart?
print('ext systems checked', n, 'bad', bad)
for chart in ['0_0_3_0', '0_0_3_1', '0_0_3_2', '0_0_3_3']:
    k = int(re.search(r"p1 array 0:.*extras (\d+)", open(f'{D}/full_c_{chart}.log').read()).group(1))
    have = len(glob.glob(f'{M}/data/hard/ext_s17/ext_{chart}_a0_*.p1.ms'))
    print(chart, 'extras', k, 'ext systems', have)
