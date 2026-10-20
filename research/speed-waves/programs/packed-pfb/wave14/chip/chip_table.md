| variant | family | slice DSP | slice LUT (INV) | CARRY | FF | DSP per signal-lane | LUT per signal-lane | KU040 channeliser DSPs | KU040 chip DSPs | KU040 LUT added (use) |
|---|---|---|---|---|---|---|---|---|---|---|
| baseline 1 product/DSP, 18-bit (as built) | xcu | 224 | 0 (0) | 0 | 0 | 14.0 | 0.0 | 1,343 (+0.0%) | 1,748 (91% full) | +0 (61.9%) |
| ap-4 packed whole branch (1 extraction) | xcu | 112 | 120 (112) | 40 | 304 | 7.0 | 7.5 | 895 (-33.4%) | 1,300 (68% full) | +480 (62.1%) |
| ap-4 packed groups of 7 (2 extractions) | xcu | 112 | 432 (112) | 160 | 928 | 7.0 | 27.0 | 895 (-33.4%) | 1,300 (68% full) | +1,728 (62.6%) |
| WP486-style per-tap extraction (14 extractions) | xcu | 112 | 4336 (112) | 1696 | 9264 | 7.0 | 271.0 | 895 (-33.4%) | 1,300 (68% full) | +17,344 (69.0%) |
| as-4 graded c2s: centre pair x2 (2 extractions) | xcu | 112 | 432 (112) | 160 | 944 | 7.0 | 27.0 | 895 (-33.4%) | 1,300 (68% full) | +1,728 (62.6%) |
| as-4 graded g3 (3 extractions) | xcu | 112 | 776 (112) | 296 | 1872 | 7.0 | 48.5 | 895 (-33.4%) | 1,300 (68% full) | +3,104 (63.2%) |
| baseline 1 product/DSP48E1, 18-bit | xc7 | 224 | 0 (0) | 0 | 0 | 14.0 | 0.0 | - | - | - |
| as-4 graded g4 on DSP48E1 (4 extractions) | xc7 | 112 | 1008 (112) | 416 | 2016 | 7.0 | 63.0 | - | - | - |
| as-4 graded g5 on DSP48E1 (5 extractions) | xc7 | 112 | 1336 (112) | 536 | 2368 | 7.0 | 83.5 | - | - | - |
| as-4 u4 on DSP48E1: centre pair unpacked + 3 graded groups | xc7 | 128 | 1144 (96) | 392 | 2272 | 8.0 | 71.5 | - | - | - |
