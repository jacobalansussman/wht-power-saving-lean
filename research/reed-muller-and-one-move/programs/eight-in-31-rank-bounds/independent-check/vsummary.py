# vsummary.py: table of claimed stage x my cascade stage from ../vverdicts.txt; writes ../extra_exclusions.txt
import os, collections
V = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = collections.Counter(); extra = []; kinds = collections.Counter()
for l in open(os.path.join(V, 'vverdicts.txt')):
    t = l.split(); c, cl, mv = t[0], t[1][6:], t[2][5:]
    st = cl.split(':')[-1] if cl != 'undecided' else 'undecided'
    T[(st, mv)] += 1
    if mv.startswith('EXCL'): kinds[(st, t[3])] += 1
    if st == 'undecided' and mv.startswith('EXCL'): extra.append(l.strip())
for st in ['stage1', 'stage2', 'stage3', 'undecided']:
    row = {k[1]: v for k, v in T.items() if k[0] == st}
    print(st, sum(row.values()), dict(sorted(row.items())))
print('kinds of my contradiction:', dict(sorted(kinds.items())))
open(os.path.join(V, 'extra_exclusions.txt'), 'w').write('\n'.join(extra) + '\n')
print('extra exclusions (claimed undecided):', len(extra))
