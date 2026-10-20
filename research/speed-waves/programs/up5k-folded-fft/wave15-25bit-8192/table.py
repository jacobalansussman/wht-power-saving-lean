#!/usr/bin/env python3
"""table.py FILE: one line per design from fpga_measure summary lines (median [min-max] over seeds)."""
import json, sys
for l in open(sys.argv[1]):
    d = json.loads(l)
    if not d.get('summary'): continue
    st = d['stats']; f = lambda k: ('%g' % st[k]['med'] + ('' if st[k]['min'] == st[k]['max'] else ' [%g-%g]' % (st[k]['min'], st[k]['max']))) if k in st else '-'
    print('%-18s seeds %s lc %s ff %s carry %s dsp %s ebr %s spram %s fmax %s harness_ok %s' % (d.get('tag'), d.get('n_ok', len(d.get('seeds', [])) if 'seeds' in d else '?'),
          f('dut_lc'), f('dut_ff'), f('dut_carry'), f('dut_dsp'), f('dut_bram'), f('dut_spram'), f('fmax_mhz'), d.get('harness_ok_all')))
