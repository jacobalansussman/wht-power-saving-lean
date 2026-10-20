#!/usr/bin/env python3
"""scale16.py -- av-5: as-6's fold_dsp6rr (16-bit halving, 4,096 points) rewritten for 8,192 points (13-bit addresses,
12-bit butterfly and twiddle counters, a 4,096-entry Q14 half-wave cosine table): the 16-bit AREA FLOOR baseline.
Pure text substitution on as-6's generated file; every substitution must hit exactly the expected number of times."""
import math, re, sys
src = open('../fold_dsp6rr.v').read()
subs = [("input [11:0] addr_in", "input [12:0] addr_in", 1), ("reg [10:0] j; reg [11:0] th;", "reg [11:0] j; reg [12:0] th;", 1),
        ("reg [11:0] bprev, aprev;", "reg [12:0] bprev, aprev;", 1), ("reg [10:0] k, step;", "reg [11:0] k, step;", 1),
        ("wire [11:0] a, oh;", "wire [12:0] a, oh;", 1), ("x < 12", "x < 13", 1), ("x == 11", "x == 12", 1),
        ("assign a[11] = th[10] ? 1'b0 : j[10]; assign oh[11] = th[10];", "assign a[12] = th[11] ? 1'b0 : j[11]; assign oh[12] = th[11];", 1),
        ("wire [11:0] b = a | oh;", "wire [12:0] b = a | oh;", 1), ("wire [11:0] caddr", "wire [12:0] caddr", 1),
        ("{2'b00, busy ? caddr : addr_in}", "{1'b0, busy ? caddr : addr_in}", 1), ("reg [10:0] radr;", "reg [11:0] radr;", 1),
        ("radr = {k[10] ^ (ph == 2'd3), k[9:0]}", "radr = {k[11] ^ (ph == 2'd3), k[10:0]}", 1),
        ("sub_r = ~k[10], sub_i = k[10]", "sub_r = ~k[11], sub_i = k[11]", 1), ("reg [15:0] rom [0:2047];", "reg [15:0] rom [0:4095];", 1),
        ("th <= 12'h7ff", "th <= 13'h0fff", 1), ("step <= 11'd1", "step <= 12'd1", 1), ("(j == 11'd0) ? 11'd0", "(j == 12'd0) ? 12'd0", 1),
        ("j <= j + 11'd1", "j <= j + 12'd1", 1), ("th <= {1'b0, th[11:1]}", "th <= {1'b0, th[12:1]}", 1),
        ("step <= {step[9:0], 1'b0}", "step <= {step[10:0], 1'b0}", 1), ("th == 12'd0", "th == 13'd0", 1),
        ("4,096-point", "8,192-point (scaled by av-5 scale16.py)", 1)]
for a, b, n in subs:
    c = src.count(a)
    assert c == n, (a, c)
    src = src.replace(a, b)
rom = ''.join("    rom[%d] = 16'h%04x;\n" % (i, int(round(16384 * math.cos(math.pi * i / 4096))) & 0xFFFF) for i in range(4096))
src, n = re.subn(r"  initial begin\n(    rom\[\d+\] = 16'h[0-9a-f]{4};\n)+  end\n", "  initial begin\n" + rom.replace('\\', '\\\\') + "  end\n", src)
assert n == 1, n
left = [l for l in src.split('\n') if re.search(r"\b(11|10|2047|1024|2048)\b|\[11:|\[10:|\[9:", l) and 'rom[' not in l and 'SB_MAC16' not in l and '.C(' not in l]
print('lines still mentioning old widths (check by eye):'); print('\n'.join(left))
open('f16_dsp6rr_8k.v', 'w').write(src)
