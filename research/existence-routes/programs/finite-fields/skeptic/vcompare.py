"""vcompare.py -- v-ex2-modp: rebuild every full-chart core of the 8 complete skeletons EXACTLY over Q (fullsysQ.py:
nag-rerun's builder with rational arithmetic, no prime, extras off; s39's linsub a6 = 1 applied exactly) and compare with
the system ex2-modp certified (raw mod-p exports reconstructed over Q, appended extras dropped)."""
import os, sys, json, re, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import *
from fractions import Fraction
M = os.environ.get('MODP_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
R = os.path.join(os.environ.get('NAG_RERUN_DIR', 'nag-rerun'), 'recipes')
UV = os.path.expanduser('~/.local/bin/uv')
def canon(names, P):
    return tuple(sorted((tuple(sorted((names[i], e) for i, e in enumerate(k) if e)), Fraction(c)) for k, c in P.items()))
def readQ(path):
    L = open(path).read().split('\n')
    if L[0].startswith('# EMPTY'): return None, []
    names = [v for v in L[0].split(',') if v]
    body = '\n'.join(L[2:]); polys = [parse_poly_str(ch, names) for ch in body.split(',\n') if ch.strip()]
    return names, polys
SK = sys.argv[1:] or ['f27_Fp2_m3_s11', 'f27_Fp2_m3_s13', 'g27_Fp2b_m3_s2', 'g27_mid3_m3_s22', 't27_mid3_m2_s2', 't27_Fmid3_m3_s5',
      't27_mid3_m2_s1', 'f27_Fp2_m3_s39', 'cal_plain4F_s1']
tot = same = 0
for sk in SK:
    V = json.load(open(f'{M}/data/dec/{sk}/verdict.json')) if os.path.exists(f'{M}/data/dec/{sk}/verdict.json') else None
    if V and 'full' in V:
        ents = [(e['chart'], e['keep']) for e in V['full']]
    else:
        ents = []
        for lg in sorted(os.listdir(f'{M}/data/dec/{sk}')):
            if lg.startswith('full_c_') and lg.endswith('.log'):
                m = re.search(r'keep=(\S+)', open(f'{M}/data/dec/{sk}/{lg}').read())
                ents.append((lg[7:-4].split('_'), m.group(1)))
    for chart, keep in ents:
        tag = '_'.join(chart)
        log = open(f'{M}/data/dec/{sk}/full_c_{tag}.log').read()
        subs = ''
        m = re.search(r"p1: linear substitutions from the relaxation: (\{.*\})", log)
        if m and m.group(1) != '{}':
            import ast
            d = ast.literal_eval(m.group(1))
            parts = []
            for v, e in d.items():
                assert list(e) == ['1'], ('non-constant substitution', e)
                parts.append(f"{int(v[1:]) - 1}={ratrec(e['1'], P1)}")
            subs = ','.join(parts)
        out = f'fq/{sk}__{tag}'
        cmd = [UV, 'run', '--quiet', '--with', 'numpy', 'python3', 'fullsysQ.py', f'{R}/{sk}.json', out, '--keep', keep, '--mode', 'n', '--arrays', '0,1'] + (['--subs', subs] if subs else [])
        r = subprocess.run(cmd, capture_output=True, text=True)
        for l in (0, 1):
            pre = f'{M}/data/dec/{sk}/full_c_{tag}.a{l}'
            if not os.path.exists(pre + '.p1.ms'):
                continue
            tot += 1
            qn, qp = readQ(f'{out}.a{l}.Q.txt')
            first = open(pre + '.p1.ms').readline()
            if first.startswith('# EMPTY'):
                ok = qn is None; print(f'{sk} {tag} a{l}: raw EMPTY CORE; exact {"EMPTY" if qn is None else "NOT EMPTY"}'); same += ok; continue
            names, F, _ = to_Q(read_ms(pre + '.p1.ms'), read_ms(pre + '.p2.ms'))
            mm = re.search(rf"p1 array {l}:.*extras (\d+)", log); k = int(mm.group(1)) if mm else 0
            Fc = F[:len(F) - k] if k else F
            A = [canon(names, P) for P in Fc]; B = [canon(qn, P) for P in qp] if qn else []
            ok = (A == B); okset = sorted(A) == sorted(B)
            same += ok or okset
            print(f'{sk} {tag} a{l}: rows raw {len(Fc)} exact {len(B)} | identical (same order) {ok} | same set {okset}' + ('' if ok or okset else '  <-- DIFFERENT'), flush=True)
print(f'TOTAL {same}/{tot} full-chart cores identical to the exact rebuild over Q')
