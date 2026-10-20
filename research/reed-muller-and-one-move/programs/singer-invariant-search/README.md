# The Singer-invariant ansatz

What these programs show: RM11 (the canonical Singer placement; `sfield.py` builds the field, a self-dual basis and
the canonical data of every form; `struct_check.py` checks the formulation exactly at blocks 2 to 5; `exact_plain.py`
checks in exact rational arithmetic that the system accepts the plain schedule and refuses broken variants;
`delcheck.py` and `uscan_sd.py` are the deletion and rank-cover tests in the two placements) and the numerical search
of section 5 (`lbf.py`, `lbfH.py`, `lma.py`, `lmc.py`, `lmd.py`: evidence only). The proofs are in
`../../proofs/singer-invariant-ansatz.md`.

Commands:

```
python3 struct_check.py > out_struct_check.txt   # 76 s; re-run for this revision, identical
python3 exact_plain.py  > out_exact_plain.txt    # about 7 s; not re-run
```

The solver runs (hours, remote 8-processor machines) were not re-run; their results are summarised in the section
README, section 5.

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/singer-invariant-search/`, where the data files of `../../data/singer-invariant-search/` sit next to the programs. Python programs need Python 3.9 or later, numpy, and torch for the solvers.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `delcheck.py` | fibre data (/K/, /J/) per H-orbit for the two placements of the Singer group (sd: self-dual basis, sigma symmetric; poly: polynomial basis), and the exact deletion test M/K/ >= W/J/ for given shapes. |
| `exact_plain.py` | EXACT check (rational arithmetic) that the canonical Singer system accepts the plain schedule and refuses broken variants. Plain schedule for one array at block m, all moves in direction 1 with ... |
| `independent-check/vfield.py` | independent rebuild of the field side of the Singer-ansatz track's Theorem S2, Props S3, S4, S5 and Remark S9. Written from the statements only . Field GF(2^m) with my own primitive polynomials; ... |
| `independent-check/vfieldvec.py` | vectorised second code path for S3/S4/S5 over ALL forms at m = 3..7. Uses F_2-linearity in the form coefficients: the move vector K(Q) (bits k: Q(a^-k)) and the target vector J(Q) (bits k: <Q, ... |
| `independent-check/vplain.py` | exact (integer) check of the canonical Singer system on the plain schedule, own code. System (Lemma C of an earlier construction track, all paid moves in direction z0 = 1): X_Q = Atil Pi_Q (I + ... |
| `independent-check/vpolyalt.py` | deletion bound max/J///K/ in the polynomial-basis placement for EVERY primitive polynomial, m=3..5 (robustness of S5(b), which is claimed for the polynomials of an earlier construction track/w8). |
| `independent-check/vsolve.py` | independent L-BFGS search of the canonical Singer system (own code; the Singer-ansatz track's lbf.py, lma.py, sfield.py NOT imported). Fibres: one representative per Z_n-orbit of nonzero forms (my ... |
| `lbf.py` | L-BFGS exploration for the canonical Singer system (same system and fibres as lma.py). All fibres are evaluated together in one batched dense solve: X_f = Atil diag(pi_f) (I + Gamma diag(pi_f))^-1 ... |
| `lbfH.py` | the L-BFGS search of lbf.py for the no-phase game with an ARBITRARY finite address group H, given by matrices over F_2 in self-dual coordinates : channels a in H; P_g/a> = /g a>; a free ... |
| `lma.py` | complex Levenberg-Marquardt with an ANALYTIC Jacobian for the canonical Singer system (same system and fibres as lmc.py / lmd.py). The residual X_f = Atil_D (I + Gamma_DD)^-1 Btil_D - I_W (x) P_J ... |
| `lmc.py` | Levenberg-Marquardt with conjugate gradients (matrix-free: torch jvp/vjp) on the canonical Singer system. H = <sigma^p> (order d), moves in field directions z_t, W arrays. For every H-orbit of ... |
| `lmd.py` | DENSE Levenberg-Marquardt (Jacobian by torch jacfwd, Marquardt scaling) on the canonical Singer system. H = <sigma^p> (order d), moves in field directions z_t, W arrays. For every H-orbit of ... |
| `sfield.py` | F_{2^m}, a self-dual basis, the Singer cycle in self-dual coordinates, binary quadratic forms, the linear part a_Q, and the canonical Singer data (K_Q, J_Q) of every form. Pure Python, exact over ... |
| `struct_check.py` | exact checks of the canonical Singer formulation at m = 2..5. (1) sigma is symmetric in self-dual coordinates; so the virtual targets sigma^-k sigma^-kT = sigma^-2k are the trace forms Tr(c x y), ... |
| `uscan_sd.py` | the decoupled multi-target rank-cover test (cell partitions) when the virtual direction counts are UNIFORM (forced when H is transitive on directions: Singer cycle |
