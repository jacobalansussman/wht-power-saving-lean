#!/usr/bin/env python3
"""tmg_paths.py REPORT.json [...] -- hand timing of the SB_MAC16 / SPRAM / EBR paths of fold_dsp6 (as-6), which
nextpnr-ice40 times with placeholder delays.  Routed net delays come from nextpnr's detailed report (the max of each
endpoint's 'delay' pair, as ao-1's dsp_paths.py and ar-4 read it); cell delays are the worst corner of icestorm's
timings_up5k.txt:
  SPRAM CLOCK->DATAOUT 1.82, DATAIN setup 0.14;  EBR RCLK->RDATA 1.18;  LC clk->Q 1.39, LUT 1.24 (worst input)
  MAC16 CLK->O 1.82 (registered outputs);  MAC_U_16X16_BYPASS setup A 6.34 / B 6.50 (multiply + 32-bit accumulate)
  ADS_U_16P16_BYPASS B->O 4.13, D->O 4.02 (combinational adder; CI taken as B);  ACC_U_16P16 setup A/B 2.76
  C/D -> 32-bit adder -> register (um2/um4): 5.26 (ap-6 / ar-4's figure);  ADDSUB setup 2.64, OHOLD 0.43, CHOLD 0.26
Each path = sum over its segments of the WORST routed delay of that segment (an upper bound over bits)."""
import json, sys, re, collections


def load(fn):
    r = json.load(open(fn))
    seg = collections.defaultdict(float)
    cnt = collections.Counter()
    for n in r['detailed_net_timings']:
        d, p = n['driver'], n['port']
        for e in n['endpoints']:
            c, q, t = e['cell'], e['port'], max(e['delay'])
            key = None
            if re.search(r'sp[ri]_RAM$', d) and p.startswith('DATAOUT'):
                if re.search(r'ub[ri]_DSP$', c):
                    key = 'spram->ubr ' + ('C/D(reg)' if q[0] in 'CD' else 'A/B/CI')
                else:
                    key = 'spram->other'
            elif re.search(r'ub[ri]_DSP$', d):
                bit = int(p.split('_')[1])
                if bit < 16:
                    key = 'ubr.sum->' + ('DATAIN LUT' if 'DATAIN' in c else c.split('.')[-1][:12])
                else:
                    key = 'ubr.d->um.A' if re.search(r'um\d_DSP$', c) else 'ubr.d->other'
            elif re.search(r'um[13]_DSP$', d):
                key = 'um1/3->um2/4 C/D' if re.search(r'um[24]_DSP$', c) else 'um1/3->other'
            elif re.search(r'um[24]_DSP$', d):
                key = 'um2/4.Y->DATAIN LUT' if 'DATAIN' in c else 'um2/4->other'
            elif 'rom' in d and p.startswith('RDATA'):
                key = 'rom->um.B' if re.search(r'um\d_DSP$', c) else 'rom->other'
            elif 'DATAIN' in d and re.search(r'sp[ri]_RAM$', c):
                key = 'DATAIN LUT->spram'
            elif re.search(r'(um\d|ub[ri])_DSP$', c):
                key = 'fabric->MAC ' + re.sub(r'_\d+$', '', q)
            elif re.search(r'sp[ri]_RAM$', c):
                key = 'fabric->spram ' + re.sub(r'_\d+$', '', q)
            if key:
                seg[key] = max(seg[key], t); cnt[key] += 1
    return r, seg, cnt


def paths(seg, comb_sum, inreg=False):
    """NOTE (checked on seed 1's critical path): the 'delay' pair of nextpnr's detailed_net_timings is the [min, max]
    ARRIVAL time at the endpoint (cumulative from the launching register), not the route alone.  For nets launched by
    a MAC / RAM output nextpnr uses a placeholder clock-to-out (about 0.1 ns), so the arrival is about the route; we add
    the real clock-to-out on top (an upper bound).  For endpoints fed by fabric logic, the arrival already holds the
    fabric clock-to-out and LUTs; where that logic is fed by a MAC output, +1.72 ns corrects the placeholder (ar-4)."""
    g = lambda k: seg.get(k)
    P = []
    def add(name, parts):
        tot = 0.0; miss = []
        for x in parts:
            if isinstance(x, str):
                v = g(x)
                if v is None:
                    miss.append(x)
                else:
                    tot += v
            else:
                tot += x
        P.append((name, tot, miss))
    if comb_sum:
        add('SPRAM -> ubr adder (comb) -> DATAIN mux LUT -> SPRAM (bound)', [1.82, 'spram->ubr A/B/CI', 4.13, -0.1, 'DATAIN LUT->spram', 0.14])
    add('MAC output register -> DATAIN mux LUT -> SPRAM (bound)', ['DATAIN LUT->spram', 1.72, 0.14])
    add('SPRAM -> ubr top adder -> d register (in the MAC)', [1.82, 'spram->ubr A/B/CI', 2.76])
    add('SPRAM -> ubr C/D input register', [1.82, 'spram->ubr C/D(reg)', 0.29])
    if inreg:   # dsp6rr: A/B input registers in the twiddle MACs (setup taken as 0.5 ns; the IM entries give 0.01)
        add('ubr d register -> um A input register', [1.82, 'ubr.d->um.A', 0.5])
        add('ROM (EBR) -> um B input register', [1.18, 'rom->um.B', 0.5])
        add('inside um: input register -> multiply (MUL_S 9.05) -> 32-bit add (ADS_32P32 5.14) -> register (pessimistic)', [9.05, 5.14])
    else:
        add('ubr d register -> um A -> multiply-accumulate -> register', [1.82, 'ubr.d->um.A', 6.50])
        add('ROM (EBR) -> um B -> multiply-accumulate -> register', [1.18, 'rom->um.B', 6.50])
    add('um1/um3 -> um2/um4 C/D -> 32-bit adder -> register', [1.82, 'um1/3->um2/4 C/D', 5.26])
    add('fabric -> MAC ADDSUB', ['fabric->MAC ADDSUBTOP', 2.64])
    add('fabric -> MAC OHOLD', ['fabric->MAC OHOLDTOP', 0.43])
    add('fabric -> MAC CHOLD', ['fabric->MAC CHOLD', 0.26])
    add('fabric -> SPRAM ADDRESS (real setup 0.27)', ['fabric->spram ADDRESS', 0.27])
    return P


if __name__ == '__main__':
    allw = []
    for fn in sys.argv[1:]:
        r, seg, cnt = load(fn)
        fm = list(r['fmax'].values())[0]['achieved']
        print('== %s  (nextpnr fmax %.2f MHz, DSP/SPRAM paths not timed by it)' % (fn.split('/')[-1], fm))
        for k in sorted(seg):
            print('   route %-28s worst %.2f ns over %d endpoints' % (k, seg[k], cnt[k]))
        w = 0
        for name, t, miss in paths(seg, comb_sum=('dsp6r' not in fn and 'dsp6' in fn), inreg=('dsp6rr' in fn)):
            print('   %-62s %6.2f ns = %5.1f MHz%s' % (name, t, 1000 / t, ('   (segment missing from report: %s)' % ', '.join(miss)) if miss else ''))
            if not miss:
                w = max(w, t)
        print('   WORST hand-timed DSP/SPRAM path %.2f ns = %.1f MHz; design bound min(nextpnr, hand) = %.1f MHz' % (w, 1000 / w, min(fm, 1000 / w)))
        allw.append(min(fm, 1000 / w))
    if allw:
        s = sorted(allw)
        print('OVER SEEDS: clock bound median %.1f [%.1f - %.1f] MHz' % (s[len(s) // 2], s[0], s[-1]))
