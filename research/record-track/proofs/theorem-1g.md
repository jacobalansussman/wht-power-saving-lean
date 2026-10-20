# Theorem 1G: mixed jobs with grouped paid blocks

This note holds the proof behind the recursion lever of the record track: an extension of the τ-theorem (in the
form argued in [../../existence-routes/](../../existence-routes/), "Theorem 1") from paid moves of rank 1 to paid
blocks of any rank. It was written on 14 October 2026 by an AI agent working for me, in a simple model of jobs, and
it is reproduced here with its notation, cleaned of internal file names.

**Status.** *Argued by hand, not checked.* Two independent AI skeptics tried to refute it on 14 October: one
re-derived the logic and simulated the construction with its own code, the other built the construction explicitly
(units, random matching inside the prescribed classes, a generic deadlock test) and checked the hypotheses of the
induction on the actual flows. Neither found an error. No person and no proof checker has read this text. A Lean
form of the theorem inside OpenAI's model was built on 19 October as part of the lever record (section 6); its
status is the record's: *proved in Lean*, as part of each of the five lever packages (in `package/leverF` and `package/leverB`, the modules under `Work/Mixed/` and `Work/MixedCost/`), with the official comparator accepting each result and three skeptics upholding each.

Labels used inside the proof: [argued] (an argument given here in full, by hand; no proof checker has read it),
[computed] (a program, named in section 4), [what-if] (a price that rests on an unchecked step), [open].

## 0. The job model

This is the model of the τ-theorem's proof, recalled so that this note can be read alone.

- A **job of size s** is: apply the Walsh-Hadamard transform H on s chosen address bits of a register holding
  2^ν entries (ν ≥ s), that is, a batch of 2^(ν-s) independent transforms of size 2^s. Its cost is counted per entry:
  T(s) is the number of full passes over the register (one per single-bit Hadamard leaf, one per free pass). A pass
  over a fraction w of the entries costs w.
- A **rule** (a mixed unit) at block m: a word, correct for every number f of digit columns, with W_k registers
  carrying a **rank-k target** (k of the m rows transformed, the other m - k untouched), paid blocks, and F free
  register-passes. Applied at column c it turns jobs of the target sizes kc into jobs of the paid sizes r_j c.
- **Split**: a job of size a + b is the job of size a followed by the job of size b (Hadamards on disjoint bits
  commute), at cost T(a) + T(b). A **plain leaf** of size s costs s passes.
- **Φ_τ** of a set of jobs (weighted by fractions of entries) is the sum of weight times size^τ.
- **τ** of a rule is the root in (0, 1) of Σ_j r_j^τ = Σ_k W_k k^τ (the *balance*). With every r_j = 1 this is the
  τ-theorem's M = Σ_k W_k k^τ.

The community's recursion is this model with one target size: a paid block is a recursive call on a whole register,
and the W arrays of the next level are W slices of that one call's batch. A rule with targets of two sizes needs jobs
of two sizes at once, which no single call's batch provides; that is what the construction below supplies.

## 1. Re-check of the rank-1 case (the τ-theorem)

- **Step 1 (multinomial flow).** With p_k = W_k/M and u(a) = multinomial(|a|; a) Π p_k^(a_k), Pascal's rule gives
  u(b) = Σ over k with b_k ≥ 1 of p_k u(b - e_k) for b ≠ 0, so M u(b) = Σ_k W_k u(b - e_k). The region
  R = {a : m c(a) ≤ L} is down-closed. A frontier point b (not in R, with some b - e_k in R) has size
  c(b) = k c(b - e_k) ≤ kL/m ≤ L and c(b) > L/m. Confirmed [argued]. An exact re-run found 0 unbalanced sizes,
  Φ(sources) = M to 10 digits, and sources in (L/3, L] [computed].
- **Step 2 (Φ conservation).** Each unit is Φ-neutral at τ, so Φ(sources) = Φ(leaves) = M. Confirmed.
- **Step 3 (free passes).** The sum of u over generation j is (W/M)^j, so the sum of u is at most M/(M - W), and the
  free passes come to F/(M - W) per unit of leaf Φ. The passes spent re-grouping portions (one per class level) cost
  O(1) per register used by a unit, so they fold into F.
- **Lemma 2 (lockstep).** Process(g, λ) issues move k of all units in the class. The outputs plus (λ/M) Source at
  generation g equal (λ/M) J at generation g, which is exactly what the child consumes. At level g there are
  M^(g_max - g) calls with λ = M^-(g_max - g), so the sources used add up to Source exactly [argued]. The literal
  simulation was re-run on 1+4 in 10 at L = 12, 30, 54, on 1+1 in 4 at L = 30, and on 7+38 in 88 at L = 12
  (7,833 classes): every class was fed exactly, the sources used equal Source, and the leaves equal M u(0)
  [computed].
- **Step 5 (top split and constant).** Every rest has size s - c(a) < (1 - 1/m)s. The induction closes with
  C' = (A + c0)/(1 - (1 - 1/m)^τ). A sharper constant is C' = (A + c0)/η_m with
  η_m = min over t in [1/m, 1] of (1 - (1-t)^τ)/t^τ, which is at least 1 - (1-1/m)^τ (Lemma A below).
- **Rounding.** The number of fluid portions is the number of lockstep classes, at most M^(g_max+1) = n^O(ln M/ln k_min).
  That is polynomial in n, not the number of lattice points, and still far below 2^(n/2) batch entries, so dyadic
  rounding costs a relative 2^-Ω(n). The claim stands; only the count in the earlier text was wrong.
- **A hypothesis that both theorems need and the earlier text did not state.** Each unit must be a correct word at
  every column count f that occurs: products of target ranks in the rank-1 case, and products of target and paid ranks
  in Theorem 1G. The community's engine also asserts m^2 > 2 max(r)^2; what that halving constraint asks of a new
  recursion was not examined here.

## 2. Lemma A: the top-split induction with leftovers (the exact condition)

**Setting.** T(n) is the per-entry cost of the best schedule for a job of size n (with at least n/2 batch bits). A
*flow for n* consists of:
- source weights w(s) ≥ 0 on sizes 1 ≤ s ≤ n, with Σ w = 1;
- exit weights x(e) ≥ 0 on sizes e < n: jobs the flow leaves over, each completed by a separate call of the whole
  procedure;
- an internal cost I ≤ A Φ̂ (free passes, re-grouping, plain size-1 leaves), where Φ̂ = Σ_s w(s) s^τ;
- a deadlock-free order, given that each exit call returns a completed job.

The schedule splits the top job into portions w(s), each portion being s + (n - s). The s-parts feed the flow, and
the (n - s)-parts and the exits are done recursively. So

    T(n) ≤ c0 + I + Σ_e x(e) T(e) + Σ_s w(s) T(n - s).

**Lemma A** [argued]. Let ψ ≥ 1 be nondecreasing, and let β > 0 and n0 be such that for every n ≥ n0 there is a flow
with

    (E)   Σ_e x(e) e^τ ψ(e) + ψ(n) Σ_s w(s) (n - s)^τ  ≤  ψ(n) (n^τ - β Φ̂).

Then T(n) ≤ C n^τ ψ(n) for all n, with C = max{ max_{n<n0} T(n)/(n^τ ψ(n)), (A + c0)/β }.

*Proof.* Induction on n; for n < n0 the claim holds by the choice of C. For n ≥ n0, the hypothesis applied to every
smaller size and (E) give T(n) ≤ c0 + A Φ̂ + C ψ(n)(n^τ - β Φ̂). Since Φ̂ ≥ 1 (sizes ≥ 1, Σ w = 1) and ψ ≥ 1,
c0 + A Φ̂ ≤ (A + c0) Φ̂ ≤ C β Φ̂ ψ(n). So T(n) ≤ C n^τ ψ(n). ∎

**Corollaries** (K is the ratio such that all sources lie in (n/K, n]; η_K as above).
- A1 (no leftovers, ψ = 1): n^τ - (n - s)^τ ≥ η_K s^τ for s in (n/K, n], so (E) holds with β = η_K. This is the
  rank-1 theorem's Step 5.
- A2 (leftovers measured by Φ, ψ = 1): (E) holds if Φ(exits) ≤ (η_K - β) Φ̂. It also holds under the earlier, weaker
  form Φ(exits) ≤ (1 - (1 - 1/K)^τ - β) n^τ, because the rests carry at most (1 - 1/K)^τ n^τ and Φ̂ ≤ n^τ. So that
  earlier argued condition is a valid sufficient condition [argued].
- A3 (small leftovers, ψ = max(1, ln)^κ): suppose every exit has ln e ≤ θ ln n with θ < 1, and Φ(exits) ≤ Φ̂. Then
  for n ≥ e^(1/θ), ψ(e) ≤ θ^κ ψ(n), so (E) holds with β = η_K/2 once θ^κ ≤ η_K/2.
- A4 (no-go for ψ = 1): if Φ(exits) = Φ̂, then (E) with ψ = 1 needs Σ_s w(s)[s^τ + (n - s)^τ - n^τ] ≤ -β Φ̂. That is
  impossible, since s^τ + (n-s)^τ ≥ n^τ. A flow whose leftovers carry all the Φ can never give exact O(n^τ) through
  this induction.

## 3. Theorem 1G: grouped paid blocks

**Rule.** Targets W_k for k in K (W = Σ W_k, m = max K, k_min = min K). Paid rounds j = 1..P, in word order, with
ranks r_j ≥ 1 (r_max = max r_j). F free passes. τ is the root in (0, 1) of Σ_j r_j^τ = Σ_k W_k k^τ.
**Separation hypothesis (S):** r_max < k_min, that is, every paid block is smaller than every job.

**The dilated lockstep tree.** Let u(a) = multinomial(|a|; a) Π_k W_k^(a_k) on N^K. Then u(0) = 1 and, for b ≠ 0,
u(b) = Σ over k with b_k ≥ 1 of W_k u(b - e_k). Let c(a) = Π k^(a_k), R(ρ) = {a : m ρ c(a) ≤ L}, and
G = max{|a| : a in R(1)}. For every round sequence s = (j_1..j_d) with d ≤ G, put ρ_s = Π_i r_(j_i). **Class s**
consists of u(a) units at column ρ_s c(a), for every a in R(ρ_s) with |a| = G - d. Round j of class s feeds class s.j
and nothing else. All sizes are integers (products of ranks), so no rounding of sizes is needed.

**Lemma B (exact balance)** [argued]. Let s' be a class and s = s'.j its child through round j, so ρ_s = ρ_(s') r_j.
At each point b with |b| = G - d(s') and b in R(ρ_(s')), the child's demand Σ_k W_k u(b - e_k) equals the parent's
round-j supply u(b). Both are jobs of size ρ_s c(b).
*Proof.* The parent's unit at b issues one job of size r_j ρ_(s') c(b) = ρ_s c(b), u(b) times. For every k with
b_k ≥ 1, m ρ_s c(b - e_k) = m ρ_(s') c(b) r_j/k ≤ L r_j/k < L, so b - e_k is in R(ρ_s) by (S). Also b ≠ 0, because
|b| ≥ 1 for every class with a child. So the demand at b is Σ over k with b_k ≥ 1 of W_k u(b - e_k) = u(b). Every
parent output lies at such a point, so every output is consumed exactly once. ∎

**Lemma C (sources only at the top)** [argued]. A demand point b of class s that is not in R(ρ_(s')) is a source of
size ρ_s c(b) in (r_j L/m, L]. For the root (G chosen maximal) the sources have size in (L/m, L].
*Proof.* m ρ_(s') c(b) > L gives ρ_s c(b) > r_j L/m. Some b - e_k lies in R(ρ_s), so ρ_s c(b) ≤ kL/m ≤ L. ∎

**Lemma D (exits)** [argued]. Outputs leave the tree only at generation 0 (d = G). There each class has at most one
unit (at column ρ_s, present when m ρ_s ≤ L), and its P outputs have sizes ρ_s r_j. Φ(exits) = Φ(sources) exactly.
Every exit has size e ≤ r_max · r_max^G ≤ r_max (L/m)^θ with θ = ln r_max / ln k_min < 1, so ln e ≤ θ ln L. Exits of
size 1 are plain leaves.
*Proof.* By Lemma B, outputs at generation ≥ 1 are all consumed. Every unit is Φ-neutral at τ, and every consumed job
was produced or sourced, so Φ(exits) = Φ(sources). G ≤ log_(k_min)(L/m) because c(a) ≥ k_min^|a|, which gives
ρ_s ≤ r_max^G ≤ (L/m)^θ. Finally r_max ≤ m^θ (since r_max = k_min^θ ≤ m^θ), which gives ln e ≤ θ ln L. ∎

**Lemma E (number of units)** [argued]. units ≤ Φ(sources) / (W k_min^τ (1 - k_min^-τ)), with no log factor.
*Proof.* A unit at depth d sits at column ≥ c(a) ≥ k_min^(G-d), and its input Φ is
Σ_k W_k (k · col)^τ ≥ W k_min^τ k_min^(τ(G-d)). The total input Φ at each depth is at most Φ(sources), since Φ only
moves down through the tree or leaves it. Summing the geometric series over G - d gives the bound. ∎

**Lemma F (order: no deadlock)** [argued]. Procedure Process(s), called with the parent's round-j outputs pending:
1. the class's inputs are exactly those outputs plus top sources (Lemmas B and C);
2. for j = 1..P: issue round j of every unit of s, running the word's free steps between moves as the word says; if
   d(s) = G, complete each exit job by a top call Do(e); otherwise call Process(s.j);
3. return.
A round j+1 is issued only after Process(s.j) has returned, which means every unit that consumed a round-j job has
finished. A unit starts only after all its inputs exist. So this is a straight-line program. Its recursion depth is
G + 1, and every exit call is a smaller instance (e < n). ∎

A careless order does deadlock: if a unit took one input from round 1 of unit A and another from round 2 of A, while
A's round 2 waits for that very unit, nothing can start. Lockstep classes, in which each child class takes only the
outputs of one round, rule this out.

**Theorem 1G** [argued, job model]. Under (S): T(n) ≤ C n^τ (ln n)^κ for all n ≥ 3, where
κ = ⌈ln(2/η_m) / ln(1/θ)⌉ (κ = 0 when r_max = 1) and C depends only on the rule. So T(n) = O(n^(τ+δ)) for every
δ > 0.
*Proof.* Take L = n. The flow is the tree normalised to source weight 1. The sources lie in (n/m, n] (Lemma C), so
K = m. The internal cost I ≤ A Φ̂ with A = (F + O(W))/(W k_min^τ (1 - k_min^-τ)) (Lemma E, plus one size-1 leaf per
unit of leaf Φ). The exits satisfy Φ(exits) = Φ̂ and ln e ≤ θ ln n (Lemma D). The order is valid (Lemma F). Lemma A3
applies with ψ = max(1, ln)^κ and β = η_m/2. The batch-bit and dyadic-portion arrangements are as in section 1: the
number of classes is at most Σ_d P^d = n^O(ln P / ln k_min). ∎

*Remarks.*
- With every r_j = 1 the tree is the rank-1 construction (ρ = 1, exits = leaves), and κ = 0 recovers O(n^τ).
- Exact O(n^τ) for grouped rules, with no log factor, would need a flow with Φ(exits) ≤ (η - β) Φ̂ (Lemma A2). This
  tree has Φ(exits) = Φ̂, so by A4 it cannot give that [open].
- The log factor costs nothing in the exponent: every saving below 1 - τ is reached with a fixed constant. OpenAI's
  goal is a big-O statement, so an existential constant is allowed.
- Without (S): Lemma B needs r_j ≤ k for the pairs that meet, and Lemma D needs r_max < k_min. The case r_max ≥ k_min
  is open. An example (targets {3:1, 2:4}, paid {1:6, 2:2}) has a fluid linear programme whose optimum was still rising
  at n = 480.

## 4. Numerical checks [computed]

Every program named here is in [../programs/theorem-1g/](../programs/theorem-1g/) and its output in
[../data/theorem-1g/](../data/theorem-1g/). All were re-run on 19 October 2026 for this revision, on an Apple
M4 Mac mini with Python 3.9, in under a minute each and under 300 MB; every output was identical to the one made on
14 October.

Three test rules, plus a control:
- A: targets {4:1, 3:2}, rounds 1,1,1,1,2,2 (τ 0.746375; θ 0.631, κ 3);
- B: targets {7:1, 4:1}, rounds 1,1,1,1,2,2,2 (τ 0.928234). This is the truncation of the pure rule {7:2} with paid
  {1:4, 2:3} plus a rank-3 idle block;
- C: targets {5:1, 4:3}, rounds 1×6, 2×2, 3×2 (τ 0.930659; θ 0.792, κ 4);
- D: the rank-1 control 1+4 in 10.

1. **Literal lockstep** (`tree_flow.py`): explicit round sequences, each round in word order. A, B and C were run at
   L = 40 and 150 (up to 259 classes). Every class was exactly fed (no surplus, no deficit), the sources used equal the
   merged count, and the exits match the merged count by size.
2. **Merged exact checks** at L = 10^3, 10^5, 10^7 (up to 1.6e10 classes, counted with multiplicity): 0 unbalanced
   points; all sources in (L/m, L]; Φ(exits)/Φ(sources) = 1.000000000000; the largest exit at or below the bound
   r_max (L/m)^θ (16,384 ≤ 21,761 for A at L = 10^7, 2,048 ≤ 2,391 for B, 177,147 ≤ 295,508 for C); the unit count
   within the bound of Lemma E.
3. **Cost of the whole schedule** (`recur.py`; F = 0, L = n, exits by recursive calls), T(n)/n^τ:

   | rule | n = 50 | n = 400 | n = 3,200 | n = 6,400 |
   |---|---|---|---|---|
   | A | 1.612 | 1.896 | 2.034 | 2.052 |
   | B | 1.144 | 1.197 | 1.218 | 1.229 |
   | C | 1.177 | 1.192 | 1.249 | 1.246 |
   | D (rank-1 control) | 1.162 | 1.180 | 1.183 | 1.183 |

   The grouped rules grow slowly, far below (ln n)^κ; A's local slope over the last doubling is 0.759 against
   τ = 0.746 and falling. The control is flat, as the rank-1 theorem says.
4. **Re-check of the rank-1 theorem** (`recheck_t1.py`, which imports the τ-theorem's own programs from
   [../../existence-routes/programs/tau-theorem-and-laser/](../../existence-routes/programs/tau-theorem-and-laser/)).
5. **The community's circuit T74, as a what-if** (`hub_whatif.py`). Control: the balance with one job size reproduces
   the community's figure for T74, 1.6797550e-3 (1.6797549561e-3). The mixed price, assuming every idle block of T74
   is a head or a tail on its own array, is 1.8169468e-3 (+8.17 %), separated (8 < 25). **This assumption turned out
   to be false**: only 11,010 of the 22,748 idle blocks are heads and none is a tail (the README of this section,
   section 4.5), so this figure is **refuted** as a price of T74. It stays here as the test it was.

The skeptics' own programs, written from this text without reusing its code:
- `indep_tree.py` rebuilds the merged tree independently (toy rule A and the T74 mixed rule up to L = 10^8):
  0 unbalanced points, Φ ratio 1, exits within the bound. `partial.py` prices partial truncation, which showed that
  the result depends on how many idle blocks really are heads.
- `explicit_sim.py` builds explicit unit objects, matches outputs to inputs at random inside the prescribed classes,
  and runs a generic topological sort of start, issue and finish events as a deadlock test independent of Lemma F;
  Φ is tracked to 50 digits. Every lockstep run is deadlock-free with Φ(exits) = Φ(sources); a deliberately careless
  matching across classes deadlocks, as it should.
- `cost_fit.py` and `check_E.py` compute the actual flows up to n = 12,800 and check hypothesis (E) of Lemma A3 at
  every n there, the unit bound of Lemma E and the exit bound; for the T74 rule `check_E.py` checks every n up to
  20,000 and every 997th n up to 3,000,000. No violation anywhere. They also
  confirm Lemma A4: (E) with ψ = 1 fails at every n, as predicted.

## 5. What the new recursion needs that the community's uniform recursion does not have

The community's recursion is one call per job: a call of size s slices its own batch into W rank-m arrays at column
s/m, and each paid block is an independent recursive call; no call ever needs a job of another size. Theorem 1G needs:
1. **Several sizes at once.** A unit at column c takes rank-73 arrays of size 73c and smaller ones (rank 49 for the
   half cut, rank 25 for the full cut) in a fixed ratio. These come from different producers' outputs, never from one
   call's batch.
2. **Fluid portions.** Classes are fractional portions of batches, which needs a portion and relabel step in the cost
   model and dyadic rounding with at least n/2 batch bits (the class count is polynomial in n).
3. **A class tree in lockstep**, ordered as in Lemma F.
4. **A top split** of the top job into frontier portions of sizes in (n/m, n] plus rests, done recursively.
5. **Exit calls.** Jobs of size ρ_s r_j are completed by fresh top calls.
6. **Mixed-target words.** The word with the head blocks removed or shortened must be proved correct for the reduced
   targets at every column count f that occurs.
7. **Hypothesis (S).** On the circuits of the record track it holds with room: for the half cut every paid rank is at
   most 48 and the smaller job is 49; for the full cut every paid rank is at most 24 and the smaller job is 25.
8. **The structure of the circuit.** Each block that is cut must be the first operation on its array, acting on that
   array alone (the truncation lemma of [../../existence-routes/](../../existence-routes/)). This is the head test
   of the README, section 4.5.

## 6. The Lean form

On 19 October the recursion was written in Lean inside OpenAI's RAM model as part of the lever record. Two things
changed against the plan in section 5:

- **No log factor.** The Lean cost bound does not use Lemma A3. It applies Lemma D at a rate point z1 and argues that
  for every exponent z above z1 the exits' share of the sources' Φ shrinks by a factor at most 2 (E m/n)^(z - z1),
  where E = r_max^(G+1) bounds the exits; since E m/n tends to 0, Lemma A closes with ψ = 1 for all large n, and a
  finite maximum covers the small ones. The statement is: if the balance holds at z1, then for every z in (z1, 1] there
  is a constant C with cost ≤ C · 2^n · (n + 1)^z. The constants exist but are astronomically large; only their
  existence is used, which OpenAI's big-O goal allows. Lemma A4 is not contradicted: it concerns the exponent τ itself.
- **The program.** Units are coded as sequences of roles (the Pascal bijection of Lemma B as arithmetic on codes);
  each paid block of a class becomes one child call whose answers are read back ("Hadamards on disjoint bits
  commute"); the leaf class hands its exits to one fresh top call; the order is the lockstep order of Lemma F, by
  induction over the word.

Status of this Lean form: *proved in Lean*, as part of each of the five lever packages (in `package/leverF` and `package/leverB`, the modules under `Work/Mixed/` and `Work/MixedCost/`), with the official comparator accepting each result and three skeptics upholding each. Its module and theorem names are listed in the lever package
(package/leverB/).
