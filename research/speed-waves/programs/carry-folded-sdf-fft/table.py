# table.py (as-5): one row per (chip, design) from the PER-SEED lines of fpga_measure.py (median [min-max]);
# per-seed lines are used because two runs were split (a stalled heavy-slot run was replaced, see the report).
import json, sys, statistics
rows = {}
import os
for fn in [a for a in sys.argv[1:] if os.path.exists(a)]:
    for ln in open(fn):
        d = json.loads(ln)
        if d.get('summary') or 'seed' not in d or not d.get('harness_ok', True):
            continue
        k = (d['chip'], d['tag'])
        rows.setdefault(k, {})[d['seed']] = d          # one line per seed (a rerun of a seed replaces it)
def f(ds, key, nd=0):
    v = [d[key] for d in ds if key in d]
    if not v: return '-'
    fmt = '%.' + str(nd) + 'f'
    m = statistics.median(v)
    return fmt % m if min(v) == max(v) else (fmt + ' [' + fmt + '-' + fmt + ']') % (m, min(v), max(v))
print('| chip | design | seeds | LUT4 positions (of which carry) | FF | slices, lower bound | slices touched | DSP | EBR | LUT-RAM | clock MHz |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for (chip, tag), sd in rows.items():
    ds = [sd[s] for s in sorted(sd)]
    dsp = f(ds, 'dut_dsp') if 'dut_dsp' in ds[0] else f(ds, 'dut_dsp_prim') + ' prim.'
    print('| %s | %s | %s | %s (%s) | %s | %s | %s | %s | %s | %s | %s |' % (chip, tag, ','.join(str(s) for s in sorted(sd)),
          f(ds, 'dut_comb'), f(ds, 'dut_comb_ccu2'), f(ds, 'dut_ff'), f(ds, 'dut_slices_lower_bound'), f(ds, 'dut_slices_touched'),
          dsp, f(ds, 'dut_bram'), f(ds, 'dut_ramw'), f(ds, 'fmax_mhz', 1)))
