#!/usr/bin/env python3
"""table.py FILE.jsonl : one row per design (the LAST summary line of each tag): cells, FF, carry, DSP, EBR, SPRAM, nextpnr MHz"""
import json, sys
rows = {}
for l in open(sys.argv[1]):
    d = json.loads(l)
    if d.get('summary'):
        rows[d['tag']] = d
for t, d in rows.items():
    s = d['stats']; m = lambda k: s[k]['med']
    f = s['fmax_mhz']
    print('%-12s cells %5d  FF %4d  carry %4d  DSP %d  EBR %2d  SPRAM %d  nextpnr %.1f [%.1f - %.1f] MHz  seeds %s failed %s harness_ok %s' % (
        t, m('dut_lc'), m('dut_ff'), m('dut_carry'), m('dut_dsp'), m('dut_bram'), m('dut_spram'), f['med'], f['min'], f['max'],
        d['seeds'], d['failed_seeds'], d['harness_ok_all']))
