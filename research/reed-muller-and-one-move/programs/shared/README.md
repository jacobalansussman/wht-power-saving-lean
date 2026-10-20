# Shared modules

`refl.py` and `refl2.py` come from an earlier round of the research (the commutative-slice families). They build
genuine numerical U_3 units (2 arrays on all 7 directions at block 3, the one known program without phase letters in
which arrays really share moves) and evaluate their matrices, so that the zero-based bounds of RM7 and RM8 can be
tested on a real program. They are imported by `../eight-in-31-rank-bounds/u3bounds.py`, `u3trace.py`, its
`independent-check/u3check.py`, and `../natural-form-star-tool/independent-check/u3true.py`.

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/shared/`, where the data files of `../../data/shared/` sit next to the programs. Python programs need Python 3.9 or later and numpy.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `refl.py` | the commutative slice as a reflection problem. Commutative slice at block m: moves t_z (z in F_2^m \ 0) generate an abelian group A; shifts t_z^2 are free; no address maps, no phase letters. ... |
| `refl2.py` | batched version of refl.py (same system, forms grouped by support size; same residual ordering). Optional subset of forms (screening: a subset is still a NECESSARY condition). |
