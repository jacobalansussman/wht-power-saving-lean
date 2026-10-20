# FPGA measuring tools

The three tools every FPGA number of the speed-first waves came from. Standard library Python 3.9 or later; the
FPGA tools themselves come from `$OSS_CAD_SUITE/bin` when that variable is set, else from `PATH`.

## fpga_measure.py: synthesis and place-and-route on 8 chips

One design through synthesis and place-and-route on one or more chips, several placement seeds, one JSON line per
seed and a summary line (median, lowest, highest).

```sh
python3 -B fpga_measure.py --list-chips
python3 -B fpga_measure.py design.v --chip ecp5-85k,nexus-lifcl40 --seeds 5 --work ./fpga_work --out results.jsonl
```

- Chips: iCE40-HX8K, iCE40-UP5K (with and without `-dsp`; **use `ice40-up5k-nodsp` for designs that instantiate
  SB_MAC16 blocks by hand**, because `synth_ice40 -dsp` silently rewrites them; the UP5K results of the speed waves
  used that setting), ECP5-25F and -85F, CrossLink-NX LIFCL-40, Gowin GW1N-9C (place-and-route), and Xilinx 7-series
  and Intel Cyclone V (synthesis only).
- **The harness.** A design has more inputs and outputs than a package has pins, so it sits in a harness: a serial
  shift register feeds all its inputs, and a signature register takes all its outputs. The harness cells are counted
  apart. A design without registers of its own gets an output register, so the timed path is register to register.
- Clocks are nextpnr's estimates. For UP5K designs with DSP blocks use `up5k_sta.py` (below).
- Optional resource wrappers: if executable scripts `heavy.sh` and `bg.sh` sit next to the tool, every tool call goes
  through them. Without them the tools run directly.
- Time and memory: from seconds (small iCE40 designs) to 30 minutes per seed (large ECP5 designs); nextpnr takes
  0.35 GB on the ECP5-85F and 0.7 GB on CrossLink-NX even for a tiny design, more for large ones.

## fpga_equiv.py: is a new design equal to a reference?

A formal check (SymbiYosys, bounded check plus k-induction) or a simulation (Icarus Verilog, corner cases first, then random
vectors), and a deliberately broken copy (`mutate`) that must fail.

```sh
python3 -B fpga_equiv.py sim   --new new.v --ref ref.v --vectors 100000
python3 -B fpga_equiv.py prove --new new.v --ref ref.v --lat-new 6
python3 -B fpga_equiv.py mutate --new new.v -o broken.v
```

A warning learned the hard way: two different adder networks compared bit by bit time out in every solver engine,
even at 8 points. Check through a helper that has the same lines at full width. `selftest` runs that recipe on the
8-point designs of the earlier circuit study (in the repository wht-fewer-operations, `tools/circuits/fpga/rtl`):
set `WHT41_RTL` to that folder.

## up5k_sta.py: timing of UP5K designs with DSP and large-RAM blocks

nextpnr-ice40 treats every non-clock pin of the UP5K's SB_MAC16 DSP block and SPRAM as a register pin with a 0.1 ns
placeholder, so a DSP used as a combinational adder or multiplier is timed far too fast. This script re-times the
routed design from nextpnr's SDF file and packed netlist, with Lattice's characterised delays from icestorm's
`timings_up5k.txt`:

```sh
nextpnr-ice40 --up5k ... --sdf D.sdf --write D.packed.json
python3 -I up5k_sta.py D.sdf D.packed.json --paths 3
```

It prints nextpnr's own figure (recomputed; it must equal nextpnr's log), the corrected maximum clock and the
critical paths. Remember also the datasheet's 50 MHz cap for an SB_MAC16 multiply without its pipeline register,
which no timing model shows. Limits: one clock domain; the DSP model is exact for the characterised configurations
and derived for mixed ones (marked in the output).

None of the three was re-run for this revision.
