# Filter-bank coefficients (section 6a)

NumPy `.npy` files (float64 arrays) used by the test benches and tables of [../../programs/packed-pfb](../../programs/packed-pfb/).

| file | what it is |
|---|---|
| `wave13/proto_h_w20_400.npy` | prototype A, the 14,336-tap low-pass prototype written by `proto.py 20 400` |
| `wave13/proto_h_B.npy` | prototype B (more alias-edge margin) |
| `wave13/c2_A_packed_s19_L1.npy` | wave 13's packed coefficient set on prototype A (-60.20 dB alias edge, -75.48 dB far stop-band) |
| `wave14/proto_h_Bs.npy` | prototype B re-centred to exact linear phase (`recentre.py`) |
| `wave14/c3_Bs_c2s_s19_c2.npy` | **the headline set**: graded tap groups on DSP48E2, ripple ±0.19 dB, alias edge -62.08 dB, far stop-band -83.18 dB |
| `wave14/final_Bs_g4_s17.npy` | the DSP48E1 (7-series) set, groups of 4: -61.72 / -78.98 dB |
| `wave14/c3_Bs_u4_s17_c2.npy` | the DSP48E1 set with the centre pair unpacked: meets every requirement at -43 % DSPs |

Each file is 114,816 bytes. The searches that wrote them are deterministic (a re-run of the headline search gave a
bit-identical file in the original work).
