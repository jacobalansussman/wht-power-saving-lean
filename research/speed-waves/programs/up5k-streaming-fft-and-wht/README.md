# Streaming 4,096-point FFT and WHT on the iCE40-UP5K (section 6c)

**FFT (`fft4096/fft4k.py`).** A 4,096-point streaming FFT that fits the chip: 4,849 logic cells (92 %), 30.9 MHz over
12 seeds, 47.7 dB signal-to-noise. A capability; its only rivals are ablations of itself (the critic's correction of the first claim that no
baseline fits). Status: checked by computer, one critic.

**WHT (`wht4096/gen.py`).** A 4,096-point streaming Walsh-Hadamard transform whose adders and phase multiplexers sit
in the 8 idle SB_MAC16 DSP blocks, one stage per block: 253 against 633 cells (-60 %), about 58 MHz in the registered
form (preset `T8r`). Status: checked by computer, one critic.

| file | what it does |
|---|---|
| `fft4096/fft4k.py` | model (precision report), generator (`--round V6` the carry-in rounding, `--mem ebr|spram`, `--rot mac|lc`, `--break twiddle|round`) and simulation check |
| `fft4096/w16_prec.py`, `lam_search.py` and their outputs | word-length and scaling studies |
| `wht4096/gen.py` | the WHT generator; presets `B0` (logic baseline), `T8b`, `T8r` (every DSP output registered), `T8Sb`, `B0S` and more; `--break-stage` writes a broken copy |
| `wht4096/simcheck.py` | RTL simulation against numpy |
| `wht4096/table.py`, `seeds15.py` | the tables |
| `wht4096/rounding/` | the rounding-mode precision study (`ha_precision.py`) |
| `*/results/` | the original measurement lines |

```sh
export PATH=$OSS_CAD_SUITE/bin:$PATH
cd fft4096
python3 fft4k.py model --frames 16                                             # -> model_precision_final.txt
W16PIPE=0 python3 fft4k.py gen -o fft_v6s.v --round V6
python3 fft4k.py sim fft_v6s.v --round V6 --samples 110000 --work fftsim       # RESULT PASS
python3 -B ../../fpga-measure/fpga_measure.py fft_v6s.v --chip ice40-up5k-nodsp --seed-list 1,2,3,4,5 --out fpga.jsonl
cd ../wht4096
python3 gen.py 12 T8r -o w12_T8r.v; python3 gen.py 12 B0 -o w12_B0.v
python3 simcheck.py w12_T8r.v 12 --samples 100000 --work sim                  # RESULT PASS
python3 -B ../../fpga-measure/fpga_measure.py w12_T8r.v --chip ice40-up5k-nodsp --seeds 5 --out fpga.jsonl
python3 -B ../../fpga-measure/fpga_measure.py w12_B0.v --chip ice40-up5k-nodsp --seeds 5 --out fpga.jsonl
python3 table.py fpga.jsonl
```

UP5K clocks with DSP blocks must be re-timed with [../fpga-measure/up5k_sta.py](../fpga-measure/up5k_sta.py); the
nextpnr figure alone is too optimistic. Needs: Python 3 with numpy; Icarus Verilog, Yosys, nextpnr-ice40. Time: the
FFT placement takes 10 to 20 minutes per seed. Memory under 1 GB.

Re-run for this revision: `fft4k.py model --frames 16` gave output identical to
`fft4096/model_precision_final.txt`. The rest was not re-run.
