# A telescope filter bank at half the multipliers (section 6a)

Two 8-bit signals share one 27 × 18 DSP48E2 multiply when they share a coefficient (the packing of Xilinx WP486).
The new part is the precision: an exact guard from each branch's or tap group's coefficient sum, an offset through
the DSP's C port, and tap groups graded by coefficient size. Synthesis only (Yosys); no Vivado place-and-route.

| folder | what it holds |
|---|---|
| `wave13/coef/` | `proto.py` (the 14,336-tap prototype filter), `quant.py` and `quant2.py` (coefficient quantisation and the packed-guard search), `table.py` (the stop-band table), `ffa_stack.py` (the fast-FIR probe of 6d), `stopband_table.md` |
| `wave13/rtl/` | `gen.py` (lane and slice RTL), `tb.py` (bit-exact test bench against numpy), `dsp48e2_sim.v` (a simulation model of the DSP48E2 subset used, written from the user guide) |
| `wave13/proof/` | the unpacking logic (`unpack_new.v`, the reference `unpack_ref.v`, a broken copy, a helper) and `exh/` (the exhaustive check of all 2^28 inputs with Verilator) |
| `wave13/chip/`, `wave13/synth/` | the whole-chip extrapolation and the synthesis statistics it reads |
| `wave14/coef/` | the graded-group search: `recentre.py` (prototype B made exactly linear-phase), `gq.py` (graded quantisation with guards), `base18.py`, `summarize.py` (the stop-band table), and helpers |
| `wave14/rtl/` | `gen2.py`, `tb2.py`, simulation models of the DSP48E1 and DSP48E2 subsets |
| `wave14/proof/exh/` | the exhaustive unpacking check for the DSP48E1 variant |
| `wave14/chip/`, `wave14/synth/` | the chip table (224 → 112 DSPs per slice; -33.4 % of the channeliser's DSPs) and the synthesis outputs |

Coefficient files (prototypes and the final quantised sets that the test benches check) are in
[../../data/packed-pfb](../../data/packed-pfb/). Copy them into the `coef/` folders before running the test benches.

```sh
export PATH=$OSS_CAD_SUITE/bin:$PATH
# wave 13: the exhaustive unpacking check (all 2^28 inputs of one DSP)
cd wave13/proof/exh
verilator --cc --exe --build -O3 -Wno-fatal --top-module top top.v ../unpack_new.v ../unpack_ref.v ../unpack_broken.v \
  ../../rtl/dsp48e2_sim.v main.cpp -o vexh && ./obj_dir/vexh
#   vectors 268435456  new-vs-ref mismatches 0  broken-vs-ref mismatches 41180672
# wave 13: prototype, quantisation, RTL test bench
cd ../../coef && python3 proto.py 20 400 && cp proto_h.npy proto_h_w20_400.npy
PROTO=A python3 quant.py w18_round w12 packed_s19_L1   # rows of the stop-band table
PROTO=A python3 quant2.py packed_s19_L1 90          # -> c2_A_packed_s19_L1.npy: -60.20 / -75.48 dB
python3 table.py > stopband_table.md
cd ../rtl && python3 gen.py && python3 tb.py packed 12 100000 ../coef/c2_A_packed_s19_L1.npy        # PASS
# wave 14: the graded groups (the headline set c3_Bs_c2s_s19_c2.npy)
cd ../../wave14/coef && python3 recentre.py proto_h_B.npy proto_h_Bs.npy
PROTO=Bs NDITH=100 DLEV=0.15 python3 gq.py c2s 19 0
PROTO=Bs INIT=c3_Bs_c2s_s19.npy G0=3070.1554022229852 KA=96 MULTI=4 WIDE=1 python3 gq.py c2s 19 60 -60.5 _c2   # -62.08 / -83.18 dB
C2=1 G4=final_Bs_g4_s17.npy python3 summarize.py > stopband_table.md
cd ../rtl && python3 gen2.py e2 c2s 12 && python3 tb2.py e2 c2s 12 100000 ../coef/c3_Bs_c2s_s19_c2.npy                  # PASS
cd ../synth && ./synth.sh glane_e2_c2s_w12.v gslice_e2_c2s_w12.v gslice_e2_c2s_w12 xcu
cd ../chip && python3 chip2.py
```

Status: wave 13 checked by computer, one critic; wave 14 checked by computer, one program, not repeated by anyone.
Needs: Python 3 with numpy and scipy; Verilator and Yosys. Time: the exhaustive check about 40 seconds; the 60-step
coefficient searches tens of minutes each; synthesis minutes. Memory under 2 GB.

Re-run for this revision: the wave-13 exhaustive check gave output identical to `wave13/proof/exh/exh.out`.
