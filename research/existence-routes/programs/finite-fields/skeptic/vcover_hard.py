import os, re, glob, ast, itertools, json
M = os.environ.get('MODP_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
ok = set()
for line in open('out_vall.txt'):
    m = re.match(r'(\S+?)/(\S+): OK ', line)
    if m: ok.add((m.group(1), m.group(2)))
for sk in ['g27_Fany4_s17', 'g27_Fany4_s19']:
    D = f'{M}/data/hard/dec/{sk}'
    L = open(f'{D}/log.txt').read()
    Z = ast.literal_eval(re.search(r'forced zeros \(0-based original alpha indices\): (\[.*?\])', L).group(1))
    rem = ast.literal_eval(re.search(r'remaining summands per move: (\[.*\])', L).group(1))
    mo = [ast.literal_eval(re.search(rf"move {t}: surviving first-nonzero options (\[.*?\])", L).group(1)) for t in range(len(rem))]
    jm = re.search(r'joint charts: (\d+) tested, (\d+) feasible in the relaxation: (\[.*\])', L)
    feas = [tuple(str(x) for x in c) for c in ast.literal_eval(jm.group(3))]
    steps = [f'fz{i+1}' for i in Z]
    for t in range(len(rem)):
        for o in rem[t] + ['zero']:
            if o not in mo[t]: steps.append(f'm{t}_{o}')
    for c in itertools.product(*mo):
        c = tuple(map(str, c))
        if c not in feas: steps.append('c_' + '_'.join(c))
    miss = [s for s in steps if (('certh2', f'{sk}__{s}') not in ok)]
    # full charts present in snapshot with GB [1]
    fulls = sorted(os.path.basename(f)[:-9] for f in glob.glob(f'{D}/full_c_*.a0.p1.ms'))
    fdone = [f for f in fulls if ('certh', f'{sk}__{f}.a0') in ok or ('certp', 'withext_' + f[7:] + '.a0') in ok]
    subs = []
    for lg in glob.glob(f'{D}/full_c_*.log'):
        m = re.search(r"p1: linear substitutions from the relaxation: (\{.*\})", open(lg).read())
        if m:
            for x in ast.literal_eval(m.group(1)):
                subs.append(f"sub_{os.path.basename(lg)[7:-4]}_{x}")
    smiss = [s for s in subs if ('certh2', f'{sk}__SUB__{s}') not in ok]
    print(f'{sk}: Z={Z} options={mo}; {len(feas)} relaxation-feasible charts of {jm.group(1)}; screen/relax steps {len(steps)}, '
          f'missing {miss}; full charts exported {len(fulls)}, certified {len(fdone)}: uncertified {sorted(set(fulls)-set(fdone))}; '
          f'linsub relations {len(subs)}, missing {smiss}')
