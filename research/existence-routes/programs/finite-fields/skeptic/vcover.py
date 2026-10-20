"""vcover.py -- v-ex2-modp: does every exclusion step of decide.py, for each skeleton claimed complete, have a
certificate that MY checker (out_vall.txt) accepted?  Steps are recomputed from verdict.json with decide.py's logic."""
import json, os, re, itertools, glob
M = os.environ.get('MODP_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
ok = set()
for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_vall.txt')):
    m = re.match(r'(\S+?)/(\S+): OK ', line)
    if m: ok.add((m.group(1), m.group(2)))
def has(skel, sys_):
    alias = {'t27_Fmid3_m3_s5': 't27F5_', 't27_mid3_m2_s1': 't27s1_'}
    cands = [(d, f'{skel}__{sys_}') for d in ('cert3', 'certm', 'cert5', 'certh', 'certh2')]
    if skel in alias:
        a = alias[skel]
        cands += [('cert', a + sys_.replace('full_c_', '')), ('cert', a + 'scr_' + sys_)]
    if skel == 'f27_Fp2_m3_s39':
        cands += [('cert5', f'{skel}__{sys_}')]
    return [c for c in cands if c in ok]
SK = ['f27_Fp2_m3_s11', 'f27_Fp2_m3_s13', 'g27_Fp2b_m3_s2', 'g27_mid3_m3_s22', 't27_mid3_m2_s2', 't27_Fmid3_m3_s5',
      't27_mid3_m2_s1', 'f27_Fp2_m3_s39']
tot = good = 0
for sk in SK:
    V = json.load(open(f'{M}/data/dec/{sk}/verdict.json'))
    ns = V['nsumm']; idx = []; b = 0
    for n in ns: idx.append(list(range(b, b + n))); b += n
    Z = V['forced_zero']
    K = [[k for k in range(n) if idx[t][k] not in Z] for t, n in enumerate(ns)]
    steps = [f'fz{i+1}' for i in Z]
    mo = V['marginal_options']
    for t in range(len(ns)):
        for opt in K[t] + ['zero']:
            if opt not in mo[t]:
                steps.append(f'm{t}_{opt}')
    feas = [tuple(c) for c in V['charts_relaxation_feasible']]
    for combo in itertools.product(*mo):
        c = tuple(map(str, combo))
        if c not in feas:
            steps.append('c_' + '_'.join(c))
    fulls = []
    for c in feas:
        fulls.append('full_c_' + '_'.join(c))
    missing = []
    for s in steps:
        tot += 1
        if has(sk, s): good += 1
        else: missing.append(s)
    for f in fulls:
        tot += 1
        if has(sk, f + '.a0') or has(sk, f + '.a1'): good += 1
        else: missing.append(f)
    # linsub substitutions
    subs = []
    for lg in glob.glob(f'{M}/data/dec/{sk}/full_c_*.log'):
        t = open(lg).read()
        if re.search(r"linear substitutions from the relaxation: \{.+\}", t):
            m = re.search(r"p1: linear substitutions from the relaxation: \{'(a\d+)'", t)
            subs.append('sub_' + os.path.basename(lg)[7:-4] + '_' + m.group(1))
    for s in subs:
        tot += 1
        if ('cert5', 'SUB__' + s) in ok: good += 1
        else: missing.append(s)
    res = [e['result'] for e in V.get('full', [])]
    print(f"{sk}: verdict={V.get('result')} | Z={Z} options={mo} | steps {len(steps)} + full {len(fulls)} + subs {len(subs)}"
          f" | missing {missing if missing else 'none'} | zero-move charts {sum('zero' in c for c in feas)}")
print(f'TOTAL {good}/{tot} steps have a certificate accepted by my independent checker')
