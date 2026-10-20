# frame-game-completion: every frame-game program is an engine unit

`complete.py` is the completion rule Λ. It turns a frame-game element into an engine element that is right in all
four column-count blocks. `lift.py` writes out a whole frame-game program as an explicit engine word: the same arrays
and paid moves in the same order, one fixed six-term free step before each paid move, corrections on the free
letters, and one scratch array for the multi-term free steps. It also has a plain evaluator, which applies each
group element to every digit column, in floating point and exactly modulo two primes. `fpa.py` is an unchanged copy
of `../frame-game/fpa.py`.

To certify a find: build it as a `FrameProg`, check `P.solve() == 'yes'` (exact readouts over the rationals), then
`word = lift(E, P)` and replay it with `run_float` and `run_modp` at small f. The usage example is at the top of
`lift.py`.

Run every command from this folder with Python 3.9 or later and numpy. The outputs of the original runs are in `../../data/frame-game-completion/` with the same names. Some programs write a JSON file next to themselves or into `out/`; run them from a copy if you want to keep the shipped files.

| program | what it does | command | time |
|---|---|---|---|
| `t_units.py` | lifts the test programs (plain schedules, exchange, two-input mixes, the order-3 controls, random programs with free letters) and replays them on raw arrays at f = 1 to fmax, and exactly modulo two primes up to fexact | `python3 t_units.py 2 8 6`; `python3 t_units.py 3 7 6`; `python3 t_units.py 4 4 3` | 1 s; seconds; minutes |
| `t_formal.py` | at block 2, the lifted words are exact identities in the engine's group algebra (all 46,080 group elements enumerated) | `python3 t_formal.py` | 17 s |
| `t_hom.py` | calibration of `complete.py`: conventions, orbits, cocycle identities, the homomorphism at f = 1 to fmax, and the failure of the f = 1 lift alone | `python3 t_hom.py 2 8`; `python3 t_hom.py 3 6` | 1 s; seconds |
| `t_dsize.py` | how many terms the free corrections have, for every letter type at blocks 2 to 5 | `python3 t_dsize.py` | 4 s |
| `t_fastpath.py` | the evaluators' fast Pauli path against explicit matrices | `python3 t_fastpath.py` | seconds |
| `t_units_cases.py` | the test programs, as a module | imported | |

The outputs named `_H` come from a first version of the rule, written in another basis. They give the same verdicts.
