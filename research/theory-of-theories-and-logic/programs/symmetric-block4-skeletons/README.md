# The 652 symmetric block-4 skeletons with every coefficient free: the skeptic's own rebuild

Section 13. Unpack the recipes first: `tar -xzf ../../data/symmetric-block4-skeletons/recipes_all678.tar.gz` (678 JSON files: the 652 skeletons and 26 calibrations, each the list of letters of one program). `vnr.py` rebuilds a skeleton's relaxation and full closure system in exact arithmetic with every coefficient free; `vnr_run.py relax NAME` decides the relaxation with Singular over the rationals (HIT when the Gröbner basis is not {1}); `chart.py NAME SECONDS` decides the full system chart by chart (a chart is excluded when Singular returns {1}; the skeleton has no unit when the last line reads `{'open': 0, 'FEASIBLE': 0}`). Before `chart.py`, run `python3 vnr_run.py fz NAME`, which writes the forced zeros to `out/NAME.fz`; `python3 vnr_run.py full NAME SECONDS` decides the full system without charts. Singular must be on the path. **msolve warning:** if you use msolve instead, write powers as `x^k`; msolve 0.6.5 reads `x*x` as `x`. The verdict tables of the original re-run are in `data/symmetric-block4-skeletons/`.

Outputs of the original runs: [data/symmetric-block4-skeletons/](../../data/symmetric-block4-skeletons/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `vnr_run.py` | the skeptic's rebuild of a symmetric skeleton (vnr.py) and Singular decision of its relaxation or full system | `tar -xzf ../../data/symmetric-block4-skeletons/recipes_all678.tar.gz && python3 vnr_run.py relax NAME` | seconds to minutes | under 1 GB | numpy, Singular | not re-run |
| `chart.py` | chart-by-chart exact decision over Q | `python3 chart.py NAME 1800` | minutes | under 1 GB | numpy, Singular | not re-run |

Other files here: `ptcheck.py`, `vnr.py`.

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
