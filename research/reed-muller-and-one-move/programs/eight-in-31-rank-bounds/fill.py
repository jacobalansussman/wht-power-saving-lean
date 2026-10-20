"""fill.py: fill the NOTE placeholders from verdicts.txt (after merge.py) and stage counts."""
import collections, sys, re
v = [l.split() for l in open('verdicts.txt')]
tot = len(v); ex = [x for x in v if x[1] != 'P']; un = [x for x in v if x[1] == 'P']
st = collections.Counter(x[-1] for x in ex)
kinds = collections.Counter(x[1] for x in ex)
s3t = int(sys.argv[1]); s3x = st.get('stage3', 0)
names = {'L': 'lower bound above upper bound', 'Z': 'star sum too small', 'M': 'matching test', 'C': 'clique test'}
kt = "; ".join("%s %d" % (names.get(k, k), n) for k, n in kinds.most_common())
rep = {'@@TOTX@@': "{:,}".format(len(ex)), '@@TOTU@@': "{:,}".format(len(un)), '@@S3T@@': str(s3t), '@@S3X@@': str(s3x),
       '@@S3L@@': str(s3t - s3x), '@@KINDS@@': kt + "."}
t = open('NOTE.md').read()
for a, b in rep.items(): t = t.replace(a, b)
open('NOTE.md', 'w').write(t)
print(rep, 'stages', dict(st))
