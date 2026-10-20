| design (prototype Bs, linear phase) | DSP per signal-lane | extractions per pair-lane | ripple ±dB (need <= 0.2) | alias edge [0.685, 1.5) ch dB (need <= -60) | worst beyond 1.5 ch dB (target <= -80) | median beyond 1.5 ch | share of stop-band <= -80 dB | legal | file |
|---|---|---|---|---|---|---|---|---|---|
| baseline: 1 product/DSP, 18-bit (as built) | 14 | 0 | 0.190 | -62.29 | -85.31 | -112.6 | 99.9 % | - | `rounded 18-bit` |
| baseline: 1 product/DSP, 12-bit | 14 | 0 | 0.190 | -61.66 | -77.45 | -98.2 | 99.8 % | - | `rounded 12-bit` |
| WP486-style per-tap packing, s=19 (12-bit) | 7 | 14 | 0.190 | -61.77 | -81.38 | -96.9 | 99.9 % | True | `c3_Bs_pt_s19_c2.npy` |
| ap-4 whole branch, s=19 | 7 | 1 | 0.191 | -61.25 | -77.56 | -92.2 | 98.8 % | True | `c3_Bs_whole_s19_c2.npy` |
| ap-4 groups of 7, s=19 | 7 | 2 | 0.190 | -61.46 | -79.29 | -94.4 | 99.6 % | True | `c3_Bs_g7_s19_c2.npy` |
| as-4 graded c2s, s=19 (DSP48E2) | 7 | 2 | 0.190 | -62.08 | -83.18 | -98.9 | 99.9 % | True | `c3_Bs_c2s_s19_c2.npy` |
| as-4 graded g3, s=19 (DSP48E2) | 7 | 3 | 0.190 | -62.21 | -82.54 | -100.2 | 99.9 % | True | `c3_Bs_g3_s19.npy` |
| WP486-style per-tap packing, s=17 (DSP48E1) | 7 | 14 | 0.193 | -60.56 | -67.28 | -83.9 | 70.7 % | True | `c3_Bs_pt_s17_c2.npy` |
| ap-4 whole branch, s=17 (DSP48E1) | 7 | 1 | 0.199 | -57.97 | -58.66 | -80.6 | 53.1 % | True | `c3_Bs_whole_s17.npy` |
| ap-4 groups of 7, s=17 (DSP48E1) | 7 | 2 | 0.197 | -60.50 | -65.56 | -83.0 | 65.1 % | True | `c3_Bs_g7_s17_c2.npy` |
| as-4 graded g4, s=17 (DSP48E1) | 7 | 4 | 0.190 | -61.72 | -78.98 | -91.8 | 98.8 % | True | `final_Bs_g4_s17.npy` |
| as-4 graded g5, s=17 (DSP48E1) | 7 | 5 | 0.190 | -61.75 | -78.76 | -91.8 | 99.1 % | True | `final_Bs_g5_s17.npy` |
| as-4 u4: centre pair unpacked + graded, s=17 (DSP48E1) | 8 | 3 | 0.190 | -62.01 | -82.69 | -97.5 | 99.9 % | True | `c3_Bs_u4_s17_c2.npy` |
