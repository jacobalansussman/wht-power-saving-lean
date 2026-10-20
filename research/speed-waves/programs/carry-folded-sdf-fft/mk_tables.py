# mk_tables.py (as-5): the report's tables with percentages against the generic R2^2SDF, from the per-seed lines.
import json, statistics, os
def load(fns):
    rows = {}
    for fn in fns:
        if not os.path.exists(fn): continue
        for ln in open(fn):
            d = json.loads(ln)
            if d.get('summary') or 'seed' not in d: continue
            rows.setdefault((d['chip'], d['tag']), {})[d['seed']] = d
    return rows
R = load(['fpga_ecp5-25k.jsonl', 'fpga_ecp5-85k.jsonl', 'fpga_ecp5-25k_64.jsonl', 'fpga_nexus-lifcl40.jsonl', 'fpga_chain_pnr.jsonl'])
def med(chip, tag, k):
    v = [d[k] for d in R.get((chip, tag), {}).values() if k in d]
    return statistics.median(v) if v else None
def rng(chip, tag, k, nd=0):
    v = [d[k] for d in R.get((chip, tag), {}).values() if k in d]
    if not v: return 'pending'
    fmt = '%.' + str(nd) + 'f'; m = statistics.median(v)
    if min(v) == max(v): return fmt % m
    return (fmt + ' [' + fmt + '-' + fmt + ']') % (m, min(v), max(v))
def pct(a, b):
    return '%+.1f %%' % (100.0 * (a - b) / b) if a is not None and b is not None else '-'
names = {'sdf': 'R2^2SDF, generic RTL (baseline)', 'sc': 'serial-butterfly SC (the idea)', 'sdfx': 'R2^2SDF, muxes in carry cells (variant)',
         'sdfn': 'R2^2SDF, muxes in carry cells (variant)', 'r22sdf': 'nanamake/r22sdf FFT64 (published, MIT)', 'chain_sdf': 'chain with R2^2SDF FFTs',
         'chain_sdfx': 'chain with fused FFTs'}
chipname = {'ecp5-25k': 'ECP5-25F', 'ecp5-85k': 'ECP5-85F', 'nexus-lifcl40': 'CrossLink-NX LIFCL-40'}
def table(spec):
    out = ['| chip, size | design | seeds | LUT4 positions | FF | slices touched | DSP | EBR | clock MHz | vs baseline: LUT4 / slices touched / clock |',
           '|---|---|---|---|---|---|---|---|---|---|']
    for chip, size, tags, base in spec:
        for t in tags:
            tag = t if t.startswith('chain') else '%s_%s' % (t, size)
            if (chip, tag) not in R:
                out.append('| %s, %s | %s | - | pending | | | | | | |' % (chipname[chip], size, names[t])); continue
            dsp = rng(chip, tag, 'dut_dsp') if 'dut_dsp' in next(iter(R[(chip, tag)].values())) else rng(chip, tag, 'dut_dsp_prim') + ' prim.'
            b = base if base.startswith('chain') else '%s_%s' % (base, size)
            cmp = '' if tag == b else '%s / %s / %s' % (pct(med(chip, tag, 'dut_comb'), med(chip, b, 'dut_comb')),
                                                         pct(med(chip, tag, 'dut_slices_touched'), med(chip, b, 'dut_slices_touched')),
                                                         pct(med(chip, tag, 'fmax_mhz'), med(chip, b, 'fmax_mhz')))
            out.append('| %s, %s | %s | %d | %s | %s | %s | %s | %s | %s | %s |' % (chipname[chip], size, names[t], len(R[(chip, tag)]),
                       rng(chip, tag, 'dut_comb'), rng(chip, tag, 'dut_ff'), rng(chip, tag, 'dut_slices_touched'), dsp,
                       rng(chip, tag, 'dut_bram'), rng(chip, tag, 'fmax_mhz', 1), cmp))
    return '\n'.join(out)
print(table([('ecp5-25k', '256', ['sdf', 'sdfx', 'sc'], 'sdf'), ('nexus-lifcl40', '256', ['sdf', 'sdfn', 'sc'], 'sdf'),
             ('ecp5-85k', '1024', ['sdf', 'sdfx', 'sc'], 'sdf'), ('ecp5-25k', '64', ['sdf', 'sdfx', 'sc', 'r22sdf'], 'sdf')]))
print()
print(table([('ecp5-85k', 'chain', ['chain_sdf', 'chain_sdfx'], 'chain_sdf')]))
