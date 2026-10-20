"""vcompare_hard.py -- v-ex2-modp: exact rebuild (fullsysQ.py, affine linsub from the logs reconstructed over Q and checked
at p2) of every full chart exported in the s17/s19 snapshot; compare with the raw cores (extras dropped)."""
import os, sys, re, ast, json, subprocess, glob
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
    return names, [parse_poly_str(ch, names) for ch in '\n'.join(L[2:]).split(',\n') if ch.strip()]
tot = same = 0
for sk in ['g27_Fany4_s17', 'g27_Fany4_s19']:
    D = f'{M}/data/hard/dec/{sk}'
    for lg in sorted(glob.glob(f'{D}/full_c_*.log')):
        tag = os.path.basename(lg)[7:-4]; log = open(lg).read()
        keep = re.search(r'keep=(\S+)', log).group(1)
        s1 = ast.literal_eval(re.search(r"p1: linear substitutions from the relaxation: (\{.*\})", log).group(1))
        s2 = ast.literal_eval(re.search(r"p2: linear substitutions from the relaxation: (\{.*\})", log).group(1))
        sub = {}
        for v, e in s1.items():
            d = {}
            for m, c in e.items():
                r = ratrec(c, P1); assert r is not None and (r.numerator - r.denominator * s2[v][m]) % P2 == 0
                key = '' if m == '1' else ','.join(str(int(w[1:]) - 1) for w in m.split('*'))
                d[key] = str(r)
            sub[str(int(v[1:]) - 1)] = d
        sf = f'fq/{sk}__{tag}.subs.json'; json.dump(sub, open(sf, 'w'))
        out = f'fq/{sk}__{tag}'
        subprocess.run([UV, 'run', '--quiet', '--with', 'numpy', 'python3', 'fullsysQ.py', f'{R}/{sk}.json', out, '--keep', keep, '--mode', 'n',
                        '--arrays', '0,1', '--subsfile', sf], capture_output=True, text=True)
        for l in (0, 1):
            pre = f'{D}/full_c_{tag}.a{l}'
            if not os.path.exists(pre + '.p1.ms'): continue
            tot += 1
            qn, qp = readQ(f'{out}.a{l}.Q.txt')
            if open(pre + '.p1.ms').readline().startswith('# EMPTY'):
                ok = qn is None; same += ok; print(sk, tag, l, 'EMPTY', ok); continue
            names, F, _ = to_Q(read_ms(pre + '.p1.ms'), read_ms(pre + '.p2.ms'))
            k = int(re.search(rf"p1 array {l}:.*extras (\d+)", log).group(1))
            Fc = F[:len(F) - k] if k else F
            A = [canon(names, P) for P in Fc]; B = [canon(qn or [], P) for P in qp]
            ok = A == B or sorted(A) == sorted(B); same += ok
            print(f'{sk} {tag} a{l}: raw core {len(Fc)} exact {len(B)} -> {"IDENTICAL" if ok else "DIFFERENT"}', flush=True)
print(f'TOTAL {same}/{tot}')
