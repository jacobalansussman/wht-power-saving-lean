| variant | slice DSP48E2 | slice LUT (of which INV) | slice CARRY | slice FF | chip PFB DSP | channeliser DSP | chip DSP (of 1,920) | chip LUT added | chip LUT use | chip FF added |
|---|---|---|---|---|---|---|---|---|---|---|
| baseline: 1 product/DSP, 18-bit coefficients (SKA-Low as built) | 224 | 0 (0) | 0 | 0 | 896 (+0%) | 1343 (+0.0%) | 1748 (+0.0%; 91% full) | +0 | 61.9% | +0 |
| precision-planned baseline: 1 product/DSP, 12-bit | 224 | 0 (0) | 0 | 0 | 896 (+0%) | 1343 (+0.0%) | 1748 (+0.0%; 91% full) | +0 | 61.9% | +0 |
| (i) packed, per-tap extraction + fabric adder trees | 112 | 4336 (112) | 1696 | 9264 | 448 (-50%) | 895 (-33.4%) | 1300 (-25.6%; 68% full) | +17,344 | 69.0% | +37,056 |
| (ii) packed, cascade groups of 4 taps (exact L1 guard), 4 extractions | 112 | 1139 (112) | 336 | 1552 | 448 (-50%) | 895 (-33.4%) | 1300 (-25.6%; 68% full) | +4,556 | 63.8% | +6,208 |
| (ii) packed, cascade groups of 7 taps (exact L1 guard), 2 extractions | 112 | 432 (112) | 160 | 928 | 448 (-50%) | 895 (-33.4%) | 1300 (-25.6%; 68% full) | +1,728 | 62.6% | +3,712 |
| packed, whole 14-tap cascade (exact L1 guard), 1 extraction | 112 | 120 (112) | 40 | 304 | 448 (-50%) | 895 (-33.4%) | 1300 (-25.6%; 68% full) | +480 | 62.1% | +1,216 |
