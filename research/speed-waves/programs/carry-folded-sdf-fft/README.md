# A streaming FFT with its multiplexers folded into carry cells (section 6b)

A radix-2² single-path delay-feedback (R2²SDF) FFT, one sample per clock, 16-bit, every butterfly halving with the
carry-in rounding fix. In the folded variant (`sdfx`), each stage's multiplexers are absorbed into the adders' Lattice
CCU2 carry cells. Surviving result: **about -20 % slices** against the critic's gated generic RTL (`sdfg`) on
CrossLink-NX at 256 points, open tools only. Status: checked by computer, one critic.

| file | what it does |
|---|---|
| `gen.py` | the generator: `--arch sdf` (generic R2²SDF), `sc` (serial commutator), `sdfx` (folded carries); `--n 6`/`8`/`10`; `--broken` perturbs one twiddle; also the bit-true model |
| `sim.py` | bit-exact simulation against the model (corner frames, then random frames); `--lib ccu2c_model.v` for `sdfx` |
| `ccu2c_model.v`, `ccu2_nexus_model.v` | simulation models of the ECP5 and Nexus carry cells, written from Yosys's cell semantics |
| `check_cells.py` | 20,000 random checks of the four fused carry-cell functions against the arithmetic |
| `count.py`, `count.txt` | the resource count by hand |
| `gen_chain.py` | a minimal range-Doppler chain: window, range FFT, corner turn, Doppler FFT, |·|² |
| `sqnr64.py` | signal-to-noise of the 64-point model |
| `table.py`, `mk_tables.py`, `table_all.md` | the tables from the measurement lines |
| `critic/gen.py` | the critic's generator with the `sdfg` arch (the same gated arithmetic written as plain RTL) |
| `critic/frame_check.py`, `critic/sim.py`, `critic/check_cells.py`, `critic/table.py` | the critic's checks and tables |
| `results/`, `critic/results/` | the original measurement lines (`fpga_nexus-lifcl40.jsonl`, `fpga_ecp5-85k.jsonl`, `critic/results/fpga_nx_256.jsonl`) and simulation results |

```sh
export PATH=$OSS_CAD_SUITE/bin:$PATH
python3 check_cells.py                                 # bad 0 / nexus bad 0
python3 gen.py --arch sdfx --n 8 -o rtl_sdfx_256.v
python3 sim.py rtl_sdfx_256.v --frames 400 --lib ccu2c_model.v
python3 gen.py --arch sdfx --n 8 --broken -o rtl_sdfx_256_broken.v && python3 sim.py rtl_sdfx_256_broken.v --frames 20 --lib ccu2c_model.v   # must FAIL
python3 critic/gen.py --arch sdfg --n 8 -o rtl_sdfg_256.v
python3 -B ../fpga-measure/fpga_measure.py rtl_sdfx_256.v --chip nexus-lifcl40,ecp5-85k --seeds 5 --out f.jsonl
python3 -B ../fpga-measure/fpga_measure.py rtl_sdfg_256.v --chip nexus-lifcl40,ecp5-85k --seeds 5 --out f.jsonl
python3 table.py f.jsonl
```

Needs: Python 3 with numpy; Icarus Verilog (or Verilator), Yosys, nextpnr. Time: simulation minutes; each seed of
place-and-route a few minutes. Memory: under 1 GB per nextpnr run at 256 points.

Re-run for this revision: `check_cells.py` printed `bad 0` and `nexus bad 0`, as in the original run. The rest was
not re-run.
