"""merge.py (r3-classes): one verdict per survivor from the stage outputs (first verdict found wins in stage order).
usage: merge.py surv8_final.txt stage1_files... -- stage2_files... -- ...   prints per stage counts, writes verdicts."""
import sys
surv = [l.strip() for l in open(sys.argv[1]) if l.strip()]
stages = []; cur = []
for a in sys.argv[2:]:
    if a == '--': stages.append(cur); cur = []
    else: cur.append(a)
stages.append(cur)
verdict = {}
for si, files in enumerate(stages):
    got = {}
    for f in files:
        for ln in open(f):
            t = ln.split()
            if len(t) >= 3 and t[0] == 'S': got.setdefault(t[1], (t[2], ln.strip()))
    n_ex = 0; n_new = 0
    for k, (v, ln) in got.items():
        if k in verdict and verdict[k][0] != 'P': continue
        if v != 'P': verdict[k] = (v, ln, si); n_new += 1
        elif k not in verdict: verdict[k] = ('P', ln, si)
    print('stage', si, 'files', len(files), 'sets seen', len(got), 'newly excluded', n_new)
ex = [k for k in surv if k in verdict and verdict[k][0] != 'P']
ps = [k for k in surv if k not in verdict or verdict[k][0] == 'P']
print('survivors', len(surv), 'excluded', len(ex), 'not excluded', len(ps))
import collections
print('by kind:', collections.Counter(verdict[k][0] for k in ex), 'by stage:', collections.Counter(verdict[k][2] for k in ex))
def line(k):
    if k in verdict and verdict[k][0] != 'P': return "%s %s stage%d\n" % (k, verdict[k][1].split(" ", 2)[2], verdict[k][2] + 1)
    return "%s P undecided (passed every stage it was given)\n" % k
open('verdicts.txt', 'w').write("".join(line(k) for k in surv))
open('undecided.txt', 'w').write("".join(k + "\n" for k in ps))
