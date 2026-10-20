# The separation toll of packed two-row arrays: why Theorem S8 was refuted as written

The wave-1 packing track claimed that separating packed two-row arrays costs at least 8 per array of the block-3
hook's room (Theorem S8), in a model whose free steps are the *Levi letters* of the two-row target R (free steps h
with R·h·R⁻¹ again free). Its skeptic refuted the proof as written and repaired a narrower statement by hand. The
skeptic's computation of the group is here.

| File | What it shows |
|---|---|
| `vlib.py` | The skeptic's own hook module (V ⊗ Λ⁴V at f = 5, built by iterated wedges, arithmetic modulo 786433). |
| `v1_struct.py` | (a) the dimensions (N 280, E_R 136, E_K 56) and the failure ranks of the plain moves; (b) the full Levi group of R in the widened game, enumerated directly: 73,728 elements up to a global phase, nine times the 8,192 the proof used; (c) the orbit of the failure kernel under the full group: 9 kernels, not the 2 × 2 grid of 4 the proof needed. |

Command (run from this folder; numpy needed):

```
python3 v1_struct.py        # about 7 minutes on a Colab machine; not re-run for this revision
```

Expected output: `../../data/separation-toll-strict-model/out_v1_struct.txt`. The skeptic's other saved outputs there:
`out_v3_twoletter.txt` (two-letter words through the full Levi group: 73,090 of the first 128,000 fail the tightness
condition) and `out_v2_W2_g16.txt`, `out_v2_W3_g16.txt` (random separable programs from the full group: the minimum
toll is exactly 8.000 per array for two and three arrays).

`v1_struct.py` writes two pickle files (`v_levi.pkl`, `v_classes.pkl`) into the folder it runs in. Comments inside the
programs use the agents' track labels ("v-ms-round2" is the skeptic); the section's Sources maps them.
