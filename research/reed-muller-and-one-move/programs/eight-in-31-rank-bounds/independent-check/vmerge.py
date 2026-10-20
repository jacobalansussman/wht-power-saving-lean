# vmerge.py (v-r3-classes): merge my vr3 cascade outputs (jobs/vr3a/remote-out/*/out_[ABC]_*.txt) with the claim's
# verdicts.txt.  Writes ../vverdicts.txt: class, claimed verdict+stage, my verdict (stage A/B/C or PASS-all / not-run).
import glob, os, sys, collections
V = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(os.path.dirname(V), 'r3-classes')
claim = {}
for l in open(os.path.join(R, 'verdicts.txt')):
    t = l.split(); claim[t[0]] = ('P' if t[1] == 'P' else t[1]) + ':' + t[-1] if t[1] != 'P' else 'undecided'
mine = {}   # class -> (stage, line)
seen = collections.defaultdict(set)
for f in sorted(glob.glob(os.path.join(V, 'jobs', 'vr3a', 'remote-out', '*', 'out_[ABC]_*.txt')) +
                glob.glob(os.path.join(V, 'jobs', 'vr3a', 'collected', 'out_[ABC]_*.txt'))):
    st = os.path.basename(f)[4]
    for l in open(f):
        t = l.split()
        if len(t) < 2: continue
        seen[t[0]].add(st)
        if t[1] != 'P':
            if t[0] not in mine or mine[t[0]][0] > st: mine[t[0]] = (st, l.strip())
out = []; tally = collections.Counter()
for c, v in claim.items():
    if c in mine: mv = 'EXCL-' + mine[c][0]
    elif 'C' in seen[c]: mv = 'PASS-ABC'
    elif seen[c]: mv = 'PASS-' + ''.join(sorted(seen[c])) + '(C not run)'
    else: mv = 'not-run'
    tally[(v.split(':')[-1] if v != 'undecided' else 'undecided', mv)] += 1
    out.append('%s claim=%s mine=%s %s' % (c, v, mv, mine[c][1].split(' ', 1)[1] if c in mine else ''))
open(os.path.join(V, 'vverdicts.txt'), 'w').write('\n'.join(sorted(out)) + '\n')
for k in sorted(tally): print(k, tally[k])
