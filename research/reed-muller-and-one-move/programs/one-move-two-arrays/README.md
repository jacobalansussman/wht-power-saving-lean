# One move, two arrays (tier 2)

What these programs show: section 5 of the section README (the best two-array programs one move short of plain miss
by 21 to 29 %) and the coverage correction of section 4, item 3. `vp.py` to `vp6.py` are the successive versions of
the solver that solves for the readout exactly at every step; `norms.py` and `sfield.py` are helpers (`sfield.py` is
the Singer-ansatz track's field). The skeptic's programs are in `independent-check/`: `vt2.py` (its own solver),
`placements.py` (the 2, 10 and 92 placements at blocks 3, 4, 5), `fibcount.py`, `gauge.py`. Lemmas Z, S and B(k) are
written out in the section README (RM12); their full proofs are in the track's working note, which is not included.

Nothing here was re-run for this revision (all runs are numerical, on remote machines).

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/one-move-two-arrays/`, where the data files of `../../data/one-move-two-arrays/` sit next to the programs. Python programs need Python 3.9 or later, numpy and torch.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `independent-check/fibcount.py` | number of H-orbits of nonzero binary quadratic forms (the fibres) in the lifted game (Singer-15 on the first 4 bits, identity on bit 5) and at blocks 4, 5. Burnside over the cyclic group; pure Python. |
| `independent-check/gauge.py` | gauge-fixed size of a block-4 2-in-7 point. The system is invariant under v_t -> lam_t v_t (B_t, Gamma_t. scaled by lam_t; Gamma_.t, A_.t by 1/lam_t). Minimise the squared parameter norm over the ... |
| `independent-check/placements.py` | how many inequivalent 'Singer games' are there at block m? The no-phase game with address group H = <sigma> (Singer) and target [I] is the group algebra C[Sym_m x/ H] with moves [zz^T] and target ... |
| `independent-check/vt2.py` | an INDEPENDENT rebuild of the no-phase Singer system of an earlier construction track Lemma C , written from the description, not from vp.py/sfield.py. Game: C[Sym_m x/ H], H = <sigma> (sigma = ... |
| `lbf.py` | L-BFGS exploration for the canonical Singer system (same system and fibres as lma.py). All fibres are evaluated together in one batched dense solve: X_f = Atil diag(pi_f) (I + Gamma diag(pi_f))^-1 ... |
| `norms.py` | total target norm sum_f W/J_f/ per fibre system (pure Python, small): Singer-(2^m-1) at m = 3, 4, 5 and the lifted group at block 5. Used only to read costs as fractions of the target. |
| `sfield.py` | F_{2^m}, a self-dual basis, the Singer cycle in self-dual coordinates, binary quadratic forms, the linear part a_Q, and the canonical Singer data (K_Q, J_Q) of every form. Pure Python, exact over ... |
| `vp.py` | VARIABLE-PROJECTION L-BFGS for the canonical Singer no-phase system . The readout Atil enters linearly, and with circulant blocks its least-squares optimum splits by Fourier frequency of the cha |
| `vp2.py` | vp.py with a fast numpy fibre builder (blocks 5 and 6). vp.py: VARIABLE-PROJECTION L-BFGS for the canonical Singer no-phase system . The readout Atil enters linearly, and with circulant blocks |
| `vp3.py` | vp2.py plus exact ALS steps on Btil (per-frequency least squares, like the A step) between L-BFGS chunks (--als k rounds of B-step + A-step per chunk). vp2.py: vp.py with a fast numpy fibre ... |
| `vp4.py` | vp3.py plus structured starts (--init plain / merge:p / del:p / merge2h:p,q,hpos). vp3.py : vp2.py plus exact ALS steps on Btil (per-frequency least squares, like the A step) between L-BFGS chunks ... |
| `vp5.py` | vp4.py on the LIFTED group (Singer of F_16 on bits 0-3, identity on bit 4) at block 5. vp4.py : vp3.py plus structured starts (--init plain / merge:p / del:p / merge2h:p,q,hpos). vp3.py : vp2.py ... |
| `vp6.py` | vp4.py plus fibre-count continuation (--sched 8,16,40: nested subsets, warm starts). vp4.py : vp3.py plus structured starts (--init plain / merge:p / del:p / merge2h:p,q,hpos). vp3.py : vp2.py ... |
