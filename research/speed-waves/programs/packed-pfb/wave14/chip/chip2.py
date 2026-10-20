# as-4 chip-level figures from the 4-signal x 4-lane slice syntheses (Yosys synth_xilinx, SYNTHESIS ONLY).
# SKA-Low station FPGA (XCKU040): [read via ap-4: Comoretto et al. arXiv 2002.10307 Table 5] 16 signals x 4 samples/clock,
#   DSPs used 1,748 of 1,920, channeliser 1,343, PFB 896 (56 per signal); LUT 150k of 242,400; FF 238k of 484,800.
#   chip = 4 slices (32 pair-lanes).
# 7-series: per signal-lane figures (the slice is 8 pair-lanes = 16 signal-lanes).
import re, os
S = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'synth')
A4 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'wave13', 'synth')
def stat(path):
    t = open(path).read().split('cells')[-1]
    d = {m.group(2): int(m.group(1)) for m in re.finditer(r"^\s+(\d+)\s+(\S+)\s*$", t, re.M)}
    lut = sum(n for k, n in d.items() if k.startswith('LUT') or k == 'INV')
    return dict(dsp=d.get('DSP48E2', 0) + d.get('DSP48E1', 0), lut=lut, inv=d.get('INV', 0),
                carry=d.get('CARRY4', 0) + d.get('CARRY8', 0), ff=d.get('FDRE', 0))
rows = [('baseline 1 product/DSP, 18-bit (as built)', f"{S}/stat_gbslice_e2_w18_xcu.txt", 'xcu'),
        ('ap-4 packed whole branch (1 extraction)', f"{A4}/stat_packed_w12_xcu.txt", 'xcu'),
        ('ap-4 packed groups of 7 (2 extractions)', f"{A4}/stat_packedg7_w12_xcu.txt", 'xcu'),
        ('WP486-style per-tap extraction (14 extractions)', f"{A4}/stat_pertap_w12_xcu.txt", 'xcu'),
        ('as-4 graded c2s: centre pair x2 (2 extractions)', f"{S}/stat_gslice_e2_c2s_w12_xcu.txt", 'xcu'),
        ('as-4 graded g3 (3 extractions)', f"{S}/stat_gslice_e2_g3_w12_xcu.txt", 'xcu'),
        ('baseline 1 product/DSP48E1, 18-bit', f"{S}/stat_gbslice_e1_w18_xc7.txt", 'xc7'),
        ('as-4 graded g4 on DSP48E1 (4 extractions)', f"{S}/stat_gslice_e1_g4_w10_xc7.txt", 'xc7'),
        ('as-4 graded g5 on DSP48E1 (5 extractions)', f"{S}/stat_gslice_e1_g5_w10_xc7.txt", 'xc7'),
        ('as-4 u4 on DSP48E1: centre pair unpacked + 3 graded groups', f"{S}/stat_gslice_e1_u4_w13_xc7.txt", 'xc7')]
print("| variant | family | slice DSP | slice LUT (INV) | CARRY | FF | DSP per signal-lane | LUT per signal-lane | KU040 channeliser DSPs | KU040 chip DSPs | KU040 LUT added (use) |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for name, p, fam in rows:
    if not os.path.exists(p): print(f"| {name} | {fam} | missing {p} |"); continue
    d = stat(p); pfb = 4 * d['dsp']; dd = pfb - 896
    chip = f"{1343 + dd:,} ({dd / 1343:+.1%}) | {1748 + dd:,} ({(1748 + dd) / 1920:.0%} full) | +{4 * d['lut']:,} ({(150000 + 4 * d['lut']) / 242400:.1%})" if fam == 'xcu' else "- | - | -"
    print(f"| {name} | {fam} | {d['dsp']} | {d['lut']} ({d['inv']}) | {d['carry']} | {d['ff']} | {d['dsp'] / 16:.1f} | {d['lut'] / 16:.1f} | {chip} |")
