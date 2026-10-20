# 4 in 15 at block 4 without phase letters, group by group

What these programs show: RM2 and RM6. `forms.py` (forms and supports), `hgroups.py` and `placements.py` (the cyclic
address groups and their placement classes relative to the target), `mtlp.py` and `mtlp2.py` (the multi-target
rank-cover linear program) with `certify.py` (rebuilds every row from its label and sums the rationalised
multipliers exactly) and `showcert.py`, `epat.py` (forced zeros plus matching), `rcover.py` (rectangle partitions),
`chain_classes.py` (all exact tests on every split for one group per class), `allsweep.py` (every direction used
once, all groups), and `solver/clm.py` (an exact-model Levenberg-Marquardt solver, evidence only). The proofs
(Theorems MT, O, E, Corollary L3, Propositions U, N, D) are in `../../proofs/four-in-15-channels-and-rank-depth.md`.

The skeptic's own programs are in `independent-check/`: `vchain.py` (classes and the four exact tests, no shared
code), `thmE.py` (Theorem E's numbers), `t10check.py` with its negative control `t10check_neg.py` (the whole chain
end to end on a real unit with address maps), `l3n.py`, `misc1.py` to `misc4.py` (class membership, symmetric
generators, 1 in 3, labels), and in `independent-check/lp/` its own linear program `mylp.py` and exact certificates
`mycert.py`.

Commands re-run for this revision (identical):

```
cd independent-check
python3 thmE.py                                    # under 1 s: rank cover 6 < slack 8 at x1x2 + x1x3
python3 vchain.py classes > out_classes.txt        # 40 s: 616, 336, 336 groups in 28, 10, 10 classes
```

Not re-run: `python3 chain_classes.py` over all classes and splits, and `vchain.py chain 4:15 <generator> <split>`
for the symmetric classes (minutes to hours each; outputs `out_chain_415.txt`, `out_chain_311.txt`, `out_chain_27.txt`
and the skeptic's `out_vchain_415_*.txt`); the solver runs of `solver/`.

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/four-in-15-by-address-group/`, where the data files of `../../data/four-in-15-by-address-group/` sit next to the programs. Python programs need Python 3.9 or later, numpy and scipy.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `allsweep.py` | for EVERY cyclic subgroup H of GL(4,2) of the given orders, 4 in 15 with the moves spread evenly (all 15 directions once: every direction gets n = /H/ virtual moves), run the three exact tests in ... |
| `certify.py` | exact, independent check of a cell-partition multi-target LP bound (mtlp.py --dump). Rebuilds every row from its label with its own form/target code (no import of forms.py); takes only the float ... |
| `chain_classes.py` | for one generator per fixed-point-free placement class of order 3 and 5 and a shape W in M, run the three exact tests on EVERY split of the M moves over the direction orbits: fixc.py (deletion + ... |
| `epat.py` | an earlier families round Corollary 10' (E-pattern matching test) on the virtual multi-target system. E[i,j] = 0 is forced when virtual moves i, j (i = j allowed) lie in a common tight support ... |
| `epat_check.py` | independent re-derivation (no shared code) of the E-pattern exclusion of 4 in 15 with the order-3 address group sigma = alpha^5 (alpha a root of x^4 + x + 1, polynomial basis), an earlier ... |
| `epat_sweep.py` | deletion + E-pattern test (epat.py) for 4 in 15 with the moves spread evenly (all 15 directions once, so every direction gets exactly n = /H/ virtual moves, for any H), over EVERY cyclic subgroup ... |
| `fixc.py` | the multi-target rank-cover test with the virtual direction counts FIXED by the physical directions. With c fixed the LP of mtlp.py decouples: it is feasible iff (i) slack(G) = c(G) - a_G >= 0 for ... |
| `fixc_classes.py` | fixc.py (deletion + rank cover, exact, all orbit splits) on one generator per placement class of the fixed-point-free cyclic groups of order 3 and 5 (out_placements.txt), shapes given on the ... |
| `forms.py` | block-m binary quadratic forms, GL(m,2) action, the virtual-target distribution of Lemma C (targets h h^T, h in H), and per-orbit support size s(F) and odd fraction a(F). Conventions: direction x ... |
| `gridcover.py` | rank cover with GRID partitions (independent row and column partitions of supp F) and forced zeros, for fixed direction counts; supports of size <= --smax (default 6). For a grid {R_i} x {C_j}: ... |
| `hgroups.py` | address groups H for mtlp*.py (exec'd; needs args, D, m) |
| `independent-check/l3n.py` | Corollary L3 (block 3, Singer Z_7) and Proposition N orbit ratios (blocks 3, 4), own code. --- L3 |
| `independent-check/lp/mycert.py` | exact LP lower bound from my own dual. Build the cell LP (as mylp.py), solve with HiGHS, rationalise the duals, and check weak duality EXACTLY with Fractions: y <= 0, A^T y <= c componentwise, ... |
| `independent-check/lp/mylp.py` | own build of the multi-target cell LP (Proposition L) with scipy HiGHS; float values. |
| `independent-check/lp/vchain.py` | independent re-implementation, written from the statements in PROOF.md only (no code shared with the wave-8 4-in-15 track). Tuple matrices over F_2, forms as dicts of coefficients, own set ... |
| `independent-check/misc1.py` | symmetric generators per order-15 class 1 in 3 at orders 3 and 5 (one-orbit-free splits: all splits of 3) |
| `independent-check/misc2.py` | helper (no description in the file) |
| `independent-check/misc3.py` | helper (no description in the file) |
| `independent-check/misc4.py` | helper (no description in the file) |
| `independent-check/sound.py` | helper (no description in the file) |
| `independent-check/t10check.py` | end-to-end check of Theorem MT + Theorem 10' on a genuine unit that USES address maps. Physical unit (1 in 4, H = Singer Z_15): r = [1]; for j = 1..4: r <- P_{g_j} [e1 e1^T] P_{g_j}^-1 r with g_j ... |
| `independent-check/t10check_neg.py` | t10check.py : end-to-end check of Theorem MT + Theorem 10' on a genuine unit that USES address maps. Physical unit (1 in 4, H = Singer Z_15): r = [1]; for j = 1..4: r <- P_{g_j} [e1 e1^T] ... |
| `independent-check/thmE.py` | check every number in Theorem E with own GF(16) (x^4 + x + 1, polynomial basis, bit i = coeff of alpha^i) |
| `independent-check/vchain.py` | independent re-implementation, written from the statements in PROOF.md only (no code shared with the wave-8 4-in-15 track). Tuple matrices over F_2, forms as dicts of coefficients, own set ... |
| `mrc.py` | minimum-rank-completion test on the tight forms (brief step 3), Z_n Singer power. Virtual moves (t,k), channel k sees direction sigma^-k z_t; virtual arrays (l,k), target sigma^-k sigma^-kT. A ... |
| `mtlp.py` | the multi-target rank-cover LP. By Lemma C of an earlier construction track, a unit of C[Sym_m x/ H] (moves z_t) gives, in the channel representation Ind chi_Q = C^H, the order-2 system (*) of an ... |
| `mtlp2.py` | multi-target rank-cover LP with general partitions (pieces R x C of supp F). Same relaxation as mtlp.py (virtual slice of the H-game: virtual targets h h^T). For each form F and each chosen ... |
| `placements.py` | classes of cyclic address groups H of GL(4,2) RELATIVE TO THE TARGET I, i.e. orbits of the orthogonal group O_4 (stabiliser of I under S -> g S g^T) acting on subgroups by conjugation. Two ... |
| `rcover.py` | rectangle-cover rank bound with forced zeros, exact, for fixed direction counts. For a form F with support P (/P/ <= --smax), rank E[D_F, D_F] = slack(F) (Theorem 10'). E vanishes on cells (x,y) ... |
| `showcert.py` | print the nonzero dual rows of an mtlp.py dump in words (rationalised multipliers). |
| `solver/clm.py` | complex Levenberg-Marquardt with an ANALYTIC Jacobian for the exact no-phase game C[Sym_m x/ H], H = <sigma> cyclic of order n (Theorem A of an earlier construction track; transfer form (**) of ... |
| `soundness_classes.py` | the three exact tests must PASS every known unit. For one generator per fixed-point-free class of order 3 and 5, run fixc.py, epat.py and rcover.py on the split of the plain schedules 1 in 4 ... |
| `sp.py` | Lemma N's Pauli-trivial block C[Sp(2m,2)] - deletion and orbit-union Theorem R on permutation modules. V = F_2^{2m}, u = (x/b) (x bits 0..m-1, b bits m..2m-1), omega(u,v) = x.b' + b.x'. Move: the ... |
| `sweep.py` | the cell multi-target LP (mtlp.py) for EVERY cyclic subgroup of GL(4,2) of the given orders (deduplicated; the LP depends only on the multiset of virtual targets h h^T). The constraint matrix is ... |
| `sweep3.py` | the cell multi-target LP (mtlp.py) at block 3 for EVERY subgroup of GL(3,2) (all subgroups are generated by at most two elements; deduplicated). Prints, per subgroup, its order, the number of ... |
| `symclass.py` | which placement classes of cyclic H (order 3, 5, 15) contain a SYMMETRIC generator (g = g^T; then every virtual target h h^T = h^2 lies in H itself, e.g. the Singer cycle written in a self-dual ... |
| `t12.py` | an earlier families round Theorem 12 (slack-s covering, via Lemma 11) on the virtual multi-target system, fixed direction counts. If a form Q_e with n_odd(Q_e) = 0 (E[D,D] nonsingular, Theorem ... |
| `trsub.py` | Theorem R over GL(m,2)-submodules of the permutation blocks of C[Sym_m x/ GL(m,2)]. Block O (a GL-orbit of binary quadratic forms): V_O = C[O], move [z z^T] = diag((-1)^{Q(z)}), target [I] = ... |
| `uscan.py` | the decoupled multi-target rank-cover test (cell partitions) when the virtual direction counts are UNIFORM (forced when H is transitive on directions: Singer cycle Z_{2^m-1}, GL). For c_x = ... |
| `whichclass.py` | which placement class (out_placements.txt) contains an earlier construction track's Singer powers (x^4 + x + 1, polynomial basis) and the six fixed-point-free order-3 survivors. |
