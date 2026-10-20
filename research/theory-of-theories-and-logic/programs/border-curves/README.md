# Logic, part 2: the 80-digit refinement of an exact border curve on a 2-in-7 skeleton

Section 15. `refine.py` reads the relaxation `g27_Fp2b_m3_s5.u.p1.ms` and the curve start `curve_starts.jsonl` from `../../data/border-curves/`. `search.py` is the skeptic's random-start search on the first-order systems; `ms2sing.py` converts the msolve-format systems to Singular input (the converted files are in `data/border-curves/sing/`). The `.ms` files are input systems, not msolve outputs; if you feed them to msolve, check first that no variable is written as a repeated product.

Outputs of the original runs: [data/border-curves/](../../data/border-curves/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `refine.py` | 80-digit refinement of the g27_Fp2b_m3_s5 border curve | `python3 refine.py` | 12 s | under 100 MB | numpy, mpmath | identical apart from timing |

Other files here: `ms2sing.py`, `myfo.py`, `search.py`.

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
