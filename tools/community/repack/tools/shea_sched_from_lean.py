"""Read Ryan Shea's published bank table (C588Rows0..17.lean, C588Sched.lean) back into a schedule JSON
(rows, pads, widths) -- the control input.  Usage: shea_sched_from_lean.py <.../Work/GCert/Data/Gen> <out.schedule.json.gz>"""
import sys, re, json, gzip, ast, os
gen, out = sys.argv[1], sys.argv[2]
rows = []
for i in range(18):
    src = open(os.path.join(gen, f'C588Rows{i}.lean')).read()
    m = re.search(r'def rows%d : List \(List Nat\) := (\[.*\])' % i, src, re.S)
    rows += ast.literal_eval(m.group(1))
src = open(os.path.join(gen, 'C588Sched.lean')).read()
pads = ast.literal_eval(re.search(r'schedPad : List Nat := (\[[0-9, ]*\])', src).group(1))
res = ast.literal_eval(re.search(r'schedRes : List Nat := (\[[0-9, ]*\])', src).group(1))
assert len(rows) == 72 and len(res) == 616
open(out, 'wb').write(gzip.compress(json.dumps(dict(rows=rows, pads=pads, widths=res), separators=(',', ':')).encode(), mtime=0))
print('rows', len(rows), 'banks', len(pads), 'res', len(res))
