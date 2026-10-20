<!-- header: existence-routes -->
# The τ-theorem, the Φ barrier and the truncation lemma

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** Theorem 1 (exact τ-transfer) and Lemma 2 (lockstep) are argued by hand, not checked; two AI skeptics of the recursion-lever round found no error in their extension to grouped paid blocks. Theorem 3 is argued by hand with no dedicated skeptic; its corollary "any saving must come from one finite unit" (section 3 and section 8 of this file) is refuted: it holds only when no paid block is bigger than a target (see [duality.md](duality.md)). Lemma 6 is argued by hand with no skeptic; its what-ifs on community circuits failed (record-track section). Propositions 4, 5 and 7 are argued by hand, not checked.

---

# ex-laser: asymptotic existence (laser method / tau-theorem) in our game

meta/ex-laser, 14 October 2026, about 07:55 to 08:50 EDT. One agent; no skeptic has read this; nothing posted; no Lean.
Labels: [proved] (proof below), [computed], [argued], [guess], [read], [what-if]. Files in this folder (scripts and
their `out_*.txt`); every job ran under a memory-capped runner except two arithmetic one-liners (RUNLOG.txt).

## 0. The question and the answer in five lines

Question (brief): does the engine's recursion admit a laser-method / tau-theorem transfer, what does a non-uniform
family of units cost, and does any existing small object give a positive asymptotic saving through it?

1. **The tau-theorem holds exactly** (Theorem 1) [proved in the job model of section 1; the engine realisation is
   argued, section 1.3]. A mixed unit (W_k arrays with rank-k targets, M rank-1 moves) gives time O(n^tau_c) per entry,
   n = log2 N, where sum_k W_k k^tau_c = M. No epsilon is lost: the only splitting happens once, at the top.
   So the literature round's targets are worth what was claimed: 1 + 4 in 10 at block 3 is worth 0.1189,
   7 + 38 in 88 is worth 0.1250.
2. **Nothing beyond the tau-theorem exists in this model** (Theorem 3, the Phi barrier) [proved]. Combining rules,
   powering them by substitution, Kronecker products, splitting: none beats the best single rule. Kronecker powers
   with any zeroing out or hashing tend to exponent 1 (Proposition 4) [proved]. So, unlike matrix multiplication,
   **powering a non-saving object can never create a saving here.** An existence proof must produce ONE finite
   (possibly mixed) unit with tau_c < 1, and it may do so non-constructively.
3. **Truncation lemma** (Lemma 6) [proved]: a saving unit that ends (or starts) some array with a block done on that
   array alone is worth strictly more as a mixed unit with the block removed. What-if on the community's #1042
   profile: 1.674989e-3 becomes 1.821393e-3 (+8.7 %), IF its 7,550 idle children are such per-array tails and IF
   Theorem 1 extends to grouped paid blocks (open, section 5) [what-if].
4. **Non-uniform families** (Proposition 7): a fixed member suffices for the Lean goal (an existence statement;
   Classical.choice is a permitted axiom). Packaging a mixed unit inside today's engine pricing (moment root, macro-word
   at block B) loses about 0.09/ln B in the exponent [computed]: B = 96 gives 0.0992 for 1 + 4 in 10 and 0.1053 for
   7 + 38 in 88. Only a new recursion (Theorem 1's schedule) gets the full value.
5. **Existing objects give nothing through the transfer** (Proposition 5): plain schedules have tau_c = 1 exactly; the
   E8 / published words are already priced at their own tau_c and no power or mixture beats it (only truncation, item 3);
   the Fano-gadget blueprint is a dual (certificate-side) object, not a word, so there is nothing to transfer.

## 1. The job model

### 1.1 Definitions

- A **job of size s** is: apply H to s chosen address bits of a register holding 2^nu entries (nu >= s), i.e. a batch of
  2^(nu-s) independent transforms of size 2^s. Its cost is counted per entry: T(s) = number of full passes over the
  register (one pass per single-bit Hadamard leaf, one per free pass). A pass over a fraction w of the entries costs w.
- A **rule** (mixed unit) at block m: a word, correct for every column count f, with W_k registers carrying a
  **rank-k target** (k of the m rows transformed across the f columns, the other m - k rows untouched), M rank-1 paid
  moves (one row across the f columns) and F register-passes of free steps. Applied at column c it turns
  {W_k jobs of size k c} into {M jobs of size c}, on registers of at least 2^(m c) entries (spare rows come from the
  register's other bits). W = sum_k W_k. Target ranks k >= 2 (K = set of ranks used, m = max K).
- **tau_c** of a rule: the root of sum_k W_k k^tau = M in (0, 1] (it exists and is unique when W < M <= sum_k W_k k).
- **Split**: a job of size a + b is the job of size a followed by the job of size b (H on disjoint bits commute), cost
  T(a) + T(b). **Plain leaf**: a job of size s costs s passes.
- **Phi_tau** of a multiset of jobs (with weights = fractions of entries) = sum weight * size^tau.

### 1.2 Why this is the engine's recursion

The engine's recursion is exactly this model for pure targets: a paid move is a recursive call on a whole register;
the W arrays of the next level are W slices of that one call's batch. The published pricing (the "moment root",
`wave16/g-stages/stages.py` root(): sum_r n_r r^tau = W m^tau for paid blocks of rank r) is Theorem 1's statement on
the PAID side. What a mixed-TARGET rule needs, and the engine does not have, is jobs of several sizes at the same time.
They cannot come from one call's batch, so the schedule must coordinate different calls. That is section 2.

### 1.3 What the realisation in the engine needs [argued]

(i) Units applied at any column c (true: units are correct for every f). (ii) A register's entries split into
fractional portions of fibres; selecting and relabelling a portion is one O(N) pass in OpenAI's RAM model. Portions are
dyadic; rounding the weights below costs a relative error of at most (number of lattice points) x 2^-(batch bits), which
is negligible once every job has at least n/2 batch bits (arranged by first splitting the top job into two halves,
a factor 2^(1-tau) in the constant). (iii) Paid moves inside a unit are sequential; Lemma 2 orders everything without
deadlock. (iv) In Lean this is a new recursion (fluid portions + lockstep phases); today's engine proves only the pure
case. None of (i)-(iv) is checked in Lean.

## 2. Theorem 1 (exact tau-transfer)

**Theorem 1** [proved in the job model]. Let a rule have target ranks K subset {2..m}, counts W_k, M rank-1 moves,
F free register-passes, and W < M. Then for every n,
  T(n) <= C n^tau_c,  C = 2^(1-tau_c) (1 + F/(M - W) + O(1)) / (1 - (1 - 1/m)^tau_c).

*Proof.*

**Step 1 (a balanced flow on the index lattice).** Index vectors a in N^K; column c(a) = prod_k k^(a_k). Different
indices may share a size; jobs are labelled by index, so that does not matter. Put p_k = W_k/M and
  u(a) = multinomial(|a|; a) prod_k p_k^(a_k)     (|a| = sum_k a_k).
A unit "at index a" works at column c(a): it consumes W_k jobs at index a + e_k (size k c(a)) for each k and produces
M jobs at index a (size c(a)). The multinomial Pascal rule gives, for a != 0,
  M u(a) = sum_k W_k u(a - e_k)     (terms with a_k = 0 omitted).
Fix a cutoff L and let R = {a : m c(a) <= L} (down-closed). Units sit at the indices of R. Then:
- a in R, a != 0: produced M u(a) = consumed sum_k W_k u(a - e_k) (all a - e_k lie in R). Balanced.
- a = 0: produced M u(0) = M jobs of size 1 = the **leaves**; nothing consumes them.
- a not in R with some a - e_k in R: consumed, never produced: a **source** of weight
  Source(a) = sum_{k: a - e_k in R} W_k u(a - e_k), and size c(a) in (L/m, L].
Checked in exact rationals for 1+4 in 10, 1+1 in 4, 7+38 in 88, 1+6 in 14 at L = 30, 200, 3000: 0 unbalanced sizes
(`frontier_flow.py`, `out_frontier_flow_480.txt`) [computed].

**Step 2 (Phi conservation).** At tau = tau_c each unit has Phi_out = M c^tau = sum_k W_k (k c)^tau = Phi_in. Summing
over all units, every consumed job is produced except the sources and every produced job is consumed except the leaves,
so sum_a Source(a) c(a)^tau_c = M (the leaves, size 1). Checked to 10 digits in the same runs [computed].

**Step 3 (free steps).** sum_a u(a) <= sum_j (sum_k p_k)^j = M/(M - W), so the free passes are at most F/(M - W) per
leaf.

**Step 4 (order: Lemma 2 below).** The flow is a valid sequential program.

**Step 5 (the top).** Take a job of size s on a register with at least n/2 batch bits (section 1.3). Give each
source index a a portion w(a) = Source(a)/S of its entries (S = sum Source, L = s), split that portion as c(a) + (s - c(a)), send the
c(a)-part into the flow and do the (s - c(a))-part by the same procedure, recursively. By Steps 2-3 the flow costs
(1 + F/(M-W)) sum_a w(a) c(a)^tau_c <= A s^tau_c with A = 1 + F/(M-W) (since c(a) <= s), and every rest has size
< (1 - 1/m) s. By induction T(s) <= C' s^tau_c as soon as A + O(1) + C'(1 - 1/m)^tau_c <= C', i.e.
C' = (A + O(1))/(1 - (1 - 1/m)^tau_c). The first halving of the batch-free top job gives the factor 2^(1-tau_c). ∎

**Computed constants** (`gamma_bound.py`, `out_gamma_bound.txt`): the induction ratio
gamma(n) = sum w c^tau / (n^tau - sum w (n-c)^tau) is at most 1.21 (1+4 in 10), 1.49 (1+1 in 4), 1.23 (7+38 in 88),
1.16 (1+6 in 14) for all n <= 10^7 (sampled), and falls slowly for large n; the analytic bound is
1/(1 - (2/3)^tau_c) = 3.33 for 1+4 in 10. The schedule's own cost (`frontier_flow.py`, without the halving, free
passes not counted): T(n)/n^tau_c = 1.160, 1.162, 1.165, 1.168 at n = 60, 120, 240, 480 for 1+4 in 10, against the
literature LP's optimum 1.096 ... 1.090 (meta/literature) [computed]. The LP is the optimum over all fluid schedules;
ours is one explicit schedule with a proof. The rising digits are pre-asymptotic: the induction caps them.

**Lemma 2 (lockstep: no deadlock)** [proved; simulated]. Generation g(a) = |a|. Units at generation g consume jobs of
generation g + 1 and produce jobs of generation g. Let J(a) = M u(a)[a in R] + Source(a) (total supply; by Step 1
J(a) = sum_k W_k u(a - e_k) = consumption). Procedure Process(g, lambda, fluid), called with
fluid = lambda J restricted to generation g + 1:
1. form lambda u(a) units at every generation-g index a; they consume the fluid exactly (balance);
2. for k = 1, ..., M: issue move k of all these units; this gives lambda u restricted to generation g; add
   lambda Source/M at the generation-g source indices; the result is (lambda/M) J restricted to generation g;
   call Process(g - 1, lambda/M, that fluid), which completes every job of move k (jobs at index 0 are leaves);
3. return.
Root: Process(g_max, 1, Source at generation g_max + 1). A move k is issued only after every job of moves 1..k-1 has
been consumed by units that have returned, and every unit's inputs exist when it starts; so the order is a valid
straight-line program. The sources used total exactly sum Source, and the leaves total M.
*Simulated literally in exact rationals* (`lockstep_sim.py`, `out_lockstep_sim.txt`): 1+4 in 10 at L = 12, 30, 54
(111, 1,111 and 11,111 classes), 1+1 in 4 at L = 30 (85 classes): every class fed exactly, no surplus or deficit,
sources used = Source, leaves = M [computed]. (A fifth case, 7+38 in 88 at L = 30, needs 88^3 classes; I stopped it
at 1.2 GB by its own pid before the 3 GB guard; not re-run.)

Why the order matters: a careless schedule deadlocks. If a unit took one input from move 1 of unit A and another from
move 2 of A, while A's move 2 waits for that very unit, nothing can start. Lockstep classes (each child class takes
only k-th moves) rule this out.

## 3. Theorem 3 (the Phi barrier: no asymptotic gain beyond the best rule)

**Theorem 3** [proved, job model]. Let the available operations be: any set of rules, each turning
{W_k jobs of size k c} into {P_r jobs of size r c} (any ranks, grouped paid blocks allowed); splits; plain leaves.
For a rule put D(tau) = sum_r P_r r^tau - sum_k W_k k^tau. If every rule used has D(tau) >= 0 at some tau in [0, 1],
then every schedule for a job of size n costs at least n^tau per entry.

*Proof.* Phi_tau never decreases: a rule changes it by c^tau D(tau) >= 0; a split by a^tau + b^tau - (a+b)^tau >= 0
(0 <= tau <= 1); and a leaf of size s costs s >= s^tau. So cost >= Phi_tau(end) >= Phi_tau(start) = n^tau. ∎

**Consequences.**
- For rules with rank-1 paid moves, D(tau) = M - sum W_k k^tau >= 0 exactly for tau <= tau_c. So the exponent of any
  schedule is >= min over the rules used of tau_c, and Theorem 1 attains it with the best rule alone. Mixing rules,
  substitution powers (they are compositions of rule applications), Kronecker products (a split followed by rule
  applications on the parts) and restrictions all live inside this model.
- **Plain rules** (M = sum W_k k, or grouped blocks that just add up the ranks) have D(tau) = sum W_k (k - k^tau) >= 0 for
  all tau <= 1: no combination of plain rules saves anything.
- This is the precise sense in which **the laser method collapses to the tau-theorem here.** In matrix multiplication
  the helper's value is supermultiplicative under the tensor product while its cost (rank) is submultiplicative, and
  restriction acts on the factors' variables. Here the multiplicative product (substitution) acts inside jobs, where
  nothing can be zeroed out (a job is a dense transform), and the product that does act on arrays (Kronecker) is
  additive in cost (Proposition 4).

## 4. Proposition 4 (Kronecker powers with hashing tend to exponent 1)

The N-fold Kronecker product of a word U (block m, W registers, M rank-1 moves) is run factor by factor: for each factor
i and each slice sigma (an assignment of the other N - 1 coordinates) run U on the W registers of the slice. After any
restriction (zeroing inputs, dropping outputs, hashing) keep a set S of tuples whose targets are rank-k targets
(k <= N m). If every kept tuple is acted on in every factor, the restricted product word still runs at least one move
per (factor, slice that contains a kept tuple), and a slice holds at most W tuples, so
  M_V >= N |S| / W,   value at tau = sum_{t in S} k_t^tau <= |S| (N m)^tau.
Hence tau_c(V) >= ln(N/W) / ln(N m) -> 1 as N -> infinity [proved]. (Kept tuples untouched in some factor only lower N.)
The exponential counting that drives the laser method (|S| ~ W''^N against a cost R^N) cancels here: cost and tuple
count share the base W. Unrestricted powers (`kron_check.py`, `out_kron_check.txt`) [computed]:
1+4 in 10: tau_N = 0.8811, 0.9360, 0.9496, 0.9603, 0.9668, 0.9709 for N = 1, 2, 3, 5, 8, 12; 7+38 in 88 and 5 in 14
the same shape.

## 5. Lemma 6 (truncation) and the open extension to grouped paid blocks

**Lemma 6** [proved]. Let a rule (paid blocks of any ranks) give some array w a target of rank k_w, where the array's
last (or first) operation is a block of rank r acting on that array alone. Removing the block and giving w the rank
k_w - r target changes D by -r^tau + k_w^tau - (k_w - r)^tau < 0 for 0 < tau < 1 (strict subadditivity of x^tau).
Since D > 0 below tau_c and D < 0 above (for a saving rule), tau_c strictly decreases: the truncated mixed rule is worth
more. Example: 5 in 14 at block 3 written as "1 + 4 rank-2 in 10, then 4 single moves" is worth 0.0628; the truncated
1 + 4 in 10 is worth 0.1189 (if the transfer of Theorem 1 is used).

**What-if for the community's best profile** (`reprice_hub.py`, `out_reprice_hub.txt`; control reproduces
stages.py's 1.674989e-3 for #1042): if each of the 7,550 idle children (ranks 1, 2, 3, 7, 16, 48) were a per-array
tail, mixed pricing would give **1.821393e-3 (+8.7 %)** [what-if]. Two things are unverified:
1. that idle children are per-array tails or heads (symmetry/NOTE.md says "each remainder child sits on one register";
   whether that register holds a single array at that point I did not check) [guess: plausible];
2. Theorem 1 for rules with GROUPED paid blocks. The flow of Step 1 then solves (W - P) * u = Source on a multiplicative
   lattice whose inverse has mixed signs, so positivity of u is not automatic. A weaker sufficient condition [argued]:
   any deficit-free flow whose excess (jobs left over and done recursively) carries a Phi fraction
   epsilon < 1 - (1 - 1/m)^tau_c of the top closes the induction, because leaves cost 1 per unit of Phi against a budget
   C. Fluid LP (`grouped_lp.py`, `out_grouped_lp.txt`): targets {3:1, 2:4}, paid {1:6, 2:2} (tau_c 0.8214): optimum
   cost/n^tau_c = 1.128, 1.128, 1.140, 1.154 at n = 60 ... 480: bounded so far but rising, so not settled [computed].
   (Two other test rules there had tau_c = 1 and test nothing.)

## 6. Proposition 5 (existing objects through the transfer)

- **Plain schedules**: tau_c = 1; Theorem 3: no power, mixture or restriction of them saves [proved].
- **E8 unit / published words** (pure targets, grouped paid blocks): their exponent is their own moment root. Theorem 3
  (with D in grouped form): no power or mixture beats the best of them [proved]. The only gain the transfer offers is
  re-pricing by truncation (Lemma 6; section 5 what-if).
- **Fano-gadget blueprint** (proofs/duality section 7): an optimal DUAL solution, a "fractional unit" bounding
  certificates from above. It is not a word and has no moves to count, so Theorems 1 and 3 do not apply to it. It says
  where a real unit's halving steps must sit, and that Theorem R can never exclude savings <= 0.0107 (block 3) or
  <= 0.0104 (block 4) [read]. It yields no saving through the transfer [argued].

## 7. Proposition 7 (non-uniform families)

(a) [read + argued] OpenAI's goal (`UniformFourier.lean`: DFTGoal := ∃ order solve W, DFTProgram ∧ TimeBounds;
permitted axioms propext, Quot.sound, Classical.choice) is an existence statement. A proof that some finite unit with
tau_c < 1 - s exists, however non-constructive, yields the program (eliminate the ∃ inside the Prop goal; no
computability is required). The user's point is right.
(b) [proved from Theorem 1] A family U_B with tau_c(U_B) -> tau*: each fixed member gives exponent tau_c(U_B) with the
constant of Theorem 1, which grows like m_B/tau (the factor 1/(1 - (1-1/m)^tau)) and like F_B/(M_B - W_B). Every s < 1 - tau*
is achieved by a fixed member; s* = 1 - tau* itself is reached only as n^(tau* + o(1)) with B growing, and then the
program must build U_B at run time. As with omega, the infimum is a limit of fixed constructions.
(c) [computed] Packaging a mixed rule into today's engine (a pure macro-word at block B, priced by the moment root,
`macro_lp.py`, `out_macro_lp.txt`): saving at B = 3, 6, 12, 24, 48, 96:
1+4 in 10: 0.0628, 0.0815, 0.0908, 0.0948, 0.0972, 0.0992 (tau_c value 0.1189);
7+38 in 88: 0.0628, 0.0847, 0.0951, 0.1001, 0.1031, 0.1053 (tau_c value 0.1250).
(tau_B - tau_c) ln B = 0.06 ... 0.09 and still drifting up: the loss is the top split (about 8 % of Phi per macro level),
paid once per macro level, so it decays only like 1/ln B. Theorem 1's schedule pays it once in total.

## 8. What this predicts, and the test

- Prediction (abstract -> concrete): every asymptotic "power" construction from a non-saving word fails (Theorem 3,
  Proposition 4); savings must come from one finite rule; mixed rules are priced exactly by tau_c; truncating
  per-array tails always helps. Tests run at the smallest sizes: exact balance and Phi (L = 30), literal lockstep
  simulation (L = 12 to 54), Kronecker squares at block 3, macro-words from B = 3.
- What would change the picture: an object outside the job model, i.e. something passed between levels that is not an
  exact transform on a set of bits (approximate / border calls in "formally polynomial" form, wave16 d-border), or a
  product under which cost is multiplicative AND arrays can be zeroed. I know of neither [argued].
