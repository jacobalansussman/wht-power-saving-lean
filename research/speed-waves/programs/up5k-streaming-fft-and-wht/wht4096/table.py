#!/usr/bin/env python3
"""table.py [fpga.jsonl] -- one row per tag from the summary lines of common/fpga_measure.py."""
import json, sys
f = sys.argv[1] if len(sys.argv) > 1 else 'fpga.jsonl'
rows = {}
for ln in open(f):
    d = json.loads(ln)
    if d.get('summary'):
        rows[d['tag']] = d
def fm(st, k, fmt='%g'):
    s = st.get(k)
    if not s: return '-'
    if s['min'] == s['max']: return fmt % s['med']
    return (fmt + ' [' + fmt + ' - ' + fmt + ']') % (s['med'], s['min'], s['max'])
print('| design | logic cells | of which carry | FF | SB_MAC16 | EBR | SPRAM | clock MHz | seeds ok | harness ok |')
print('|---|---|---|---|---|---|---|---|---|---|')
for t, d in rows.items():
    st = d['stats']
    print('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (t, fm(st, 'dut_lc'), fm(st, 'dut_carry'), fm(st, 'dut_ff'),
          fm(st, 'dut_dsp'), fm(st, 'dut_bram'), fm(st, 'dut_spram'), fm(st, 'fmax_mhz', '%.1f'),
          '%d/%d' % (len(d['seeds']) - len(d['failed_seeds']), len(d['seeds'])) if isinstance(d['seeds'], list) else d['seeds'], d.get('harness_ok_all')))
