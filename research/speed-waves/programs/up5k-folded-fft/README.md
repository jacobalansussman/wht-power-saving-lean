# Folded FFTs on the iCE40-UP5K (section 6c)

**4,096 points, 16-bit (`ffold.py`).** An in-place radix-2 decimation-in-frequency FFT, every stage halving with the
carry-in rounding fix, all arithmetic in 6 SB_MAC16 DSP blocks, data in the large single-port RAMs. Surviving result:
235 logic cells, -40 % against the best DSP-multiplier rival, 50 MHz (the datasheet's cap). Status: checked by
computer, one critic.

**8,192 points, 25-bit (`wave15-25bit-8192/gen24.py`, wave 15).** The same plan widened so that nothing wraps: 274
cells, the noise floor read within 0.73 dB. Status: checked by computer, one program, not repeated by anyone.

| file | what it does |
|---|---|
| `ffold.py` | bit-true model, Verilog generator (variants `dsp6rr` the headline, `dsp6r`, `dsp6`, `dspbf`, `dsp`, `lc1`, `lcs`), simulation check, precision report (`prec`); `BRK_TW=1` / `BRK_CARRY=1` write broken copies |
| `lcs_core.v.tmpl` | the datapath of the logic-only serial rival |
| `tmg_paths.py` | hand timing of the DSP and RAM paths |
| `table.py`, `precision.txt`, `results/fpga.jsonl` | the tables and the original measurements |
| `critic/` | the critic's tables (`fpga.jsonl`, `fpga_b.jsonl`, `fpga_c.jsonl`), timing reports (`results/tmg/`) and `partial_check.py` (checks netlist frames as they finish; run it from this folder's parent so it finds `ffold.py`) |
| `wave15-25bit-8192/` | `gen24.py` (model, generator, simulation), `dr.py` to `dr6.py` (dynamic-range tests and their outputs), `scale16.py` and `sim16.py` (the 16-bit design scaled to 8,192 points, for comparison), `table.py` |

```sh
export PATH=$OSS_CAD_SUITE/bin:$PATH
python3 ffold.py prec                                         # identical to precision.txt
python3 ffold.py gen --var dsp6rr -o fold_dsp6rr.v
python3 ffold.py sim fold_dsp6rr.v --frames 9 --work sim_a    # RESULT PASS
BRK_TW=1 python3 ffold.py gen --var dsp6rr -o brk_tw.v && python3 ffold.py sim brk_tw.v --frames 8 --work sim_b   # FAIL
python3 -B ../fpga-measure/fpga_measure.py fold_dsp6rr.v --chip ice40-up5k-nodsp --seeds 5 --out fpga.jsonl
python3 table.py fpga.jsonl                                   # 235 cells
# re-time the DSP paths: run nextpnr-ice40 with --sdf and --write, then ../fpga-measure/up5k_sta.py
cd wave15-25bit-8192
python3 gen24.py gen --var dsp --log 13 -o f24_dsp.v
python3 gen24.py sim f24_dsp.v --log 13 --frames 4 --work sim24
python3 dr.py                                                 # +0.73 dB (25-bit) against +33.3 dB (16-bit) at -90 dBFS
python3 -B ../../fpga-measure/fpga_measure.py f24_dsp.v --top f24 --chip ice40-up5k-nodsp --seeds 5 --out fpga.jsonl
```

Needs: Python 3 with numpy; Icarus Verilog, Yosys, nextpnr-ice40. Time: an RTL simulation of 9 frames takes minutes;
each seed of place-and-route 1 to 15 minutes. Memory under 1 GB.

Re-run for this revision: `ffold.py prec` gave output identical to `precision.txt`. The rest was not re-run.
