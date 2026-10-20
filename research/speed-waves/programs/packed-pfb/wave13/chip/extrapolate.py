# Chip-level extrapolation for the SKA-Low station FPGA (Xilinx KU040) from the 4-signal x 4-lane slice synthesis.
# [read] Comoretto et al., arXiv 2002.10307, Table 5 and sec. 8.1: per FPGA 16 signals, 4 samples/clock; DSP48 used 1,748
#        (91 %), channeliser 1,343, PFB 56 multipliers per signal (= 896); LUT 150k (62 %), FF 238k (49 %), BRAM 535 (89 %).
# [recalled] KU040 totals: 242,400 LUT, 484,800 FF, 1,920 DSP48E2 (consistent with the paper's percentages).
# [ran] slice figures: ../synth/stat_<variant>_w<W>_xcu.txt (Yosys 0.69 synth_xilinx -family xcu; SYNTHESIS ONLY).
# The slice is 8 pair-lanes (4 signals x 4 lanes); the chip has 32 pair-lanes (16 signals x 4 lanes): factor 4.
import re, os
S = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'synth')
dsp_used, dsp_chan, dsp_avail, pfb_base = 1748, 1343, 1920, 896
lut_used, lut_avail, ff_used, ff_avail = 150_000, 242_400, 238_000, 484_800
def stat(v):
    t = open(f"{S}/stat_{v}_xcu.txt").read().split('cells')[-1]
    d = {m.group(2): int(m.group(1)) for m in re.finditer(r"^\s+(\d+)\s+(\S+)\s*$", t, re.M)}
    lut = sum(n for k, n in d.items() if k.startswith('LUT') or k == 'INV')
    return d.get('DSP48E2', 0), lut, d.get('INV', 0), d.get('CARRY4', 0) + d.get('CARRY8', 0), d.get('FDRE', 0)
rows = [('base_w18', 'baseline: 1 product/DSP, 18-bit coefficients (SKA-Low as built)'),
        ('base_w12', 'precision-planned baseline: 1 product/DSP, 12-bit'),
        ('pertap_w12', '(i) packed, per-tap extraction + fabric adder trees'),
        ('packedg4_w12', '(ii) packed, cascade groups of 4 taps (exact L1 guard), 4 extractions'),
        ('packedg7_w12', '(ii) packed, cascade groups of 7 taps (exact L1 guard), 2 extractions'),
        ('packed_w12', 'packed, whole 14-tap cascade (exact L1 guard), 1 extraction')]
print("| variant | slice DSP48E2 | slice LUT (of which INV) | slice CARRY | slice FF | chip PFB DSP | channeliser DSP | chip DSP (of 1,920) | chip LUT added | chip LUT use | chip FF added |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for v, name in rows:
    d, l, inv, c, f = stat(v)
    pfb = 4 * d; dd = pfb - pfb_base
    print(f"| {name} | {d} | {l} ({inv}) | {c} | {f} | {pfb} ({dd/pfb_base:+.0%}) | {dsp_chan+dd} ({dd/dsp_chan:+.1%}) | "
          f"{dsp_used+dd} ({dd/dsp_used:+.1%}; {(dsp_used+dd)/dsp_avail:.0%} full) | +{4*l:,} | {(lut_used+4*l)/lut_avail:.1%} | +{4*f:,} |")
