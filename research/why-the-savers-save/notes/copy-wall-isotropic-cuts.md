> Working note written by an AI agent during the round of 13 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# d-isocut (agent d-g-break-1, 13 Oct 2026): confined totals at isotropic cuts on complex P112

d-g-break-1 is a symlink to this folder. Everything here is our own code; no hub script was run.

## Step 0: the realisability rule adopted (gate)

Sources read: hub format-limits-corrected (rule table rows 1-4, 6, 8; A.1 (Bat), (Sch); C.1 sandwich; C.4 Lemma Q),
s2-decoder/REPORT.md item 2.

(a) Degenerate frames. Rule-table row 2: the spec allows any invertible D_v; nondegenerate subspaces are only the
    paper's instance; #130's general frame lemma covers degenerate binary subspaces with d(U,V) = dim U + dim V -
    2 dim(U cap V). Adopted: h-perp and <h> (and any t.i. g, g-perp) are legal frames, priced by d.
(b) A copy that ends at g != 0. Row 4: exact restoration is required of helpers and dirty registers; scratch is
    arbitrary by spec ("rho fixes scratch" is convention). Adopted, two-sided: OPTIMISTIC = terminal frame g costs
    nothing; CONSERVATIVE = the scratch must be back at its call-start frame, +d(g, 0) = dim g per copy per call.
    Under the conservative charge a line cut saves 7 + 1 - 8 = 0, so confined totals gain nothing at all.
(c) Reader chains. Rows 1 and 3 (common frame at every gate; climbs along a register's life, drops only as copies)
    give: a target reads copy j only while its own frame equals the copy's frame; its frame climbs 0 -> u_S-perp.
    So the cuts of the copies one target reads form a chain. A copy may branch and climb at paid cost d(g, g').
    Consequence: with 1-dim cuts a target reads cut-0 (star) copies plus copies of AT MOST ONE isotropic line <h>;
    several cuts only along a totally isotropic flag <h> < <h,h'> < ... inside u_S-perp (at most one cut per dimension).

Cost model (#446 Theorem C as used in the hub): total j costs dim P_j - dim g_j, g_j subset of P_j cap (served span)-perp.
Equivalently cost_j = rank_F2 B0[served_j, read_j]. Cut-0 totals are star totals.

## Step 1: counting (exact unless marked)

`s1_classes.py` (8 s) -> out/s1_classes.txt: every t.i. cut g (dims 1..4: 255 / 5355 / 11475 / 2295 subspaces):
|C_g|, omega, octads inside, octads met, a partition of C_g into orthogonal cliques (certificate).
Line cuts (|C_h|, omega, clique-partition parts): 48/8/6, 52/8/7, 56/7/8 (the 8 special h), 56/8/7, 60/8/9.
So |T| <= 9 omega(T) for every T inside any C_h.

`s1_k8.py` -> out/s1_k8.txt: the 8 omega-7 classes are exactly the isotropic h orthogonal to all 8 removed-octad
labels; their complement has omega 4. Complement of C_h: omega 8 for the 48/52/56-omega-8 classes, omega 4 for the
56-omega-7 and the 60 classes.

Lemma (M1). Model M1: chain rule (c); every total is either a star total with read span of dim >= 8 (cost >= 8) or a
line-cut total at an isotropic <h> whose read span is h-perp (cost 7). Then every fitting decoder with K >= 9 totals
costs >= 64. Proof: cost = 8 K0 + 7 K1 >= 7K + K0 >= 63 + K0, so cost <= 63 forces K0 = 0 and K = 9. With no star
totals each target reads exactly one class h(S); the class-h totals give rank B[T_h, T_h] <= K_h, and an orthogonal
clique in T_h forces K_h >= omega(T_h) >= |T_h| / 9 (clique-partition certificate). Summing, K >= ceil(112/9) = 13 > 9.
Pure line-class designs therefore need K >= 13 totals (any read spans), and cost >= 91 = 7 x 13 with full spans.
K = 8 (cost 56 + K0 <= 63 would beat the wall): killed EXACTLY for the 48/52/56-omega-8 classes with K0 = 7 (the
complement holds an octad, rank <= 7 there); for the 120 classes whose complement has omega 4 counting does not
kill it. Numerically: `s2_als.py` (masked ALS, P112, pattern cx) -- control `coord` fits at residual 6e-12;
`free8` (8 unmasked totals) converges to a natural B = D1 Gram D2 on 6/6 restarts (rank-8 uniqueness, now seen on
P112, not only E8); `k8special`, `k8c60`, `k8c48`, `k8two` (7 or 6 free totals + 1 or 2 confined) stay at residual
6.6 / 7.2 / 7.0 / 9.9 on 6/6 restarts each. Natural totals cannot be confined (u-ours (ii), exact). So K = 8 is dead
up to the numerical uniqueness.

Not closed: (i) flags (2-, 3-, 4-dim t.i. cuts, cost 5 / 3 / 1 with full spans): the chain rule allows one cut per
dimension per target, and the per-bottom-class clique count does not bound cheap high-dim cuts (a 4-dim cut total
reads a single orthogonal clique C_g and decodes it at x = 1). Plain octad counting is far too weak here (met/cost up
to 501, out/s1_classes.txt). (ii) totals with partial read spans (cheap stars, partial class totals): this contains
the open general star floor (OPEN-72) and d-sparse's search.

## Pricing (pricer control reproduced: control_b_matches_799 true, b 1.6644003e-3; out/pricer_control.txt)
- What-if, nothing realised: 9 x 7 in place of {8:576} -> b 1.7106e-3 (+2.77 %), out/pricer_9x7.txt.
  By the Lemma no M1 decoder realises this profile.

## Is a pure line-class decoder realisable at all? (unfinished, numerical)
`s2_pureclass.py extra restarts iters`: classes h = 5, 9, 12 (|C_h| 60 each, covering P112), T_h sizes 60/26/26,
omega 8/4/4, K_h = omega + extra, reads C_h, serves T_h. Best residuals: K 16: 3.7; K 22: 0.17; K 28: 9.3e-4;
K 40: 6.7e-3 (800 iterations, 2 restarts). No candidate below 1e-10 yet. Even if one converges it costs >= 7 x 28 = 196
with full spans: a realisability check of the mechanism, not a gain.
Resume: `cd checks/wht43/copywall/d-isocut && checks/wht43/bound-check/g-p112-decoder/venv/bin/python -I -W ignore s2_pureclass.py 4 4 5000`
