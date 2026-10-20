<!-- header: existence-routes -->
# OpenAI's price language

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** Proposition 3 and the 0.036 ceiling are refuted (its hypothesis is inconsistent at every exponent below 1; [skeptic-price-reformulation.md](skeptic-price-reformulation.md)). Proposition 1 is open: argued, with its key step not written out. Lemma 2 and the sparsity check are right (the skeptic).

---

# ex2-reform: proofs and arguments

14 October 2026, evening. Track R of `../BRIEF.md`. No skeptic has read this. Nothing posted. Labels: [proved] = a
complete argument is written here; [argued] = an argument that leaves a step to an unchecked in-house theorem or
to a check not done; [computed] = `calib.py` / `sparse/` with their outputs; [hypothesis] = assumed, not proved.

## 0. The rule world (the setting of Proposition 1)

- A **type** is a kernel K in GL_q(C) that is not monomial (not a permutation times a diagonal). A **task**
  [K, n] asks for K^{(x)n} applied to one array. Its normalised cost c_K(n) is the cost per coordinate.
  The plain row-column method gives c_K(n) <= n c_K(1).
- A **rule** is one exact word, right for every number f of digit columns. Its targets are tasks
  [A_j, m_j f] on live arrays, and its paid calls are tasks [B_t, r_t f]. Scalar mixing of arrays and monomial
  letters (applied f-fold) are free. This is our engine's notion of a unit, with several kernels allowed.
  - Example: the engine's W in M at block m has W targets [H, m f] and M calls [H, f], where H is the 2x2 Hadamard.
  - A block of rank r is one call [H, r f].
  - Mixed targets are targets [H, k f] with k < m.
  - A conversion such as "U from 3 calls of H" is a rule with target [U, f] and three calls [H, f].
- A **program** applies rules recursively. A call [B, n'] is either served by a rule or is a base call of
  bounded size. The **exponent** z_K is the infimum of the z for which some program has c_K(n) = O(n^z).
  The saving is 1 - z.
- A **theta-price** is a function v from types to (0, infinity) with v(U) = 1, U = [[1,0],[1,1]], such that every
  rule satisfies

      sum over targets j of  v(A_j) m_j^theta   <=   sum over calls t of  v(B_t) r_t^theta.        (R_theta)

## 1. Proposition 1 (theta-price duality)

**Statement.** Fix theta in (0, 1]. Some type has a program of exponent below theta if and only if no
theta-price exists. Below theta, read "with exponent at most theta' for some theta' < theta"; boundary cases are
handled by continuity in theta.

**(a) A theta-price rules out exponents below theta** [proved]. Give every pending task [K, n] the potential
v(K) n^theta. By (R_theta), applying a rule never lowers the total potential. A base call of bounded size has
potential at most a constant P_0 and costs at least one operation per coordinate. A program for [K, n] therefore
makes at least v(K) n^theta / P_0 base calls per coordinate, so c_K(n) >= (v(K)/P_0) n^theta.

**(b) Finitely many types and rules, no theta-price: a combination of rules that saves** [proved, standard].
Write rule r as a vector a_r over the types, with a_r(K) = (target weight of type K) - (call weight of type K),
each weight being sum of size^theta. The condition (R_theta) reads a_r . v <= 0. Motzkin's transposition
theorem (B = I, C = -A) says that exactly one of the following holds:
- some v > 0 has a_r . v <= 0 for every r;
- some lambda >= 0 has sum_r lambda_r a_r >= 0 componentwise and not zero.

The second alternative is a mixture of rules. Run in the proportions lambda, it produces at least as much weighted
target as it consumes in weighted calls for every type, and strictly more for at least one.

**(c) Realising the mixture as a program** [argued; this is the step that needs a skeptic]. For one type, this is
ex-laser's Theorem 1 (the tau-transfer, `../../meta/ex-laser/PROOF.md`; no skeptic yet). That proof has three
parts: split the top job once, feed units at many sizes in a balanced flow, and run their moves in lockstep. It
gives a program of exponent tau with sum_k W_k k^tau = M, with no loss. With several types the flow is a
multi-commodity version of the same construction. I have not written it out.
- OpenAI's Lemma 4.2 (packing) is their realisation step, in their tensor-power model. It does not preserve the
  exponent: the win it builds is for a block-diagonal J whose width w must be taken large (for the random schedule),
  with T <= (b - eps) w^(b-1) calls; amplification gives alpha = ln T / ln(b w^(b-1)), so the saving is about
  (eps/b) / ln(b w^(b-1)) and shrinks like 1/ln w [argued from their section 2.2 formula].

**(d) Infinitely many types** [argued]. Fix conversion rules for each type K: K from C_K calls of U, and
U (x) I_q from e_K calls of K (positive generation, OpenAI Lemma 2.1). Under (R_theta) they force
v(K) in [c_K, C_K] with c_K > 0. Every finite subfamily of rules has a theta-price by (b), and the product of these
intervals is compact (Tychonoff). This is OpenAI's section 4.2 argument word for word.

**Corollaries.**
- **C1 (one type)** [proved, given (a) to (c)]. With only H, v(H) cancels from (R_theta). The condition is
  W m^theta <= sum_t r_t^theta for every unit: Theorem C's z formula (meta-spectrum). With mixed targets it is
  ex-laser's sum_k W_k k^tau = M. `calib.py` Part A reproduces 14 recorded values to 4 digits [computed].
  Consequence: in the engine's own game the dual object is a single number, so duality says nothing new there
  (meta-spectrum's finding, re-derived).
- **C2 (theta = 1)** [read]. (R_1) with OpenAI's normalisation is their comparison inequality (4.2):
  sum_i Delta_i p(A_i) <= 0, where Delta = plain calls minus actual calls. In their model their Theorem 4.1 is
  "no finite win implies a 1-price, with exact laws", and their main theorem is "no 1-price". So s > 0 is the
  statement "the theta = 1 point does not extend from the Kronecker family to all kernels".
- **C3 (universality, iso section 2b)** [computed on a two-type toy, argued in general]. Conversion rules of
  shrink 1 in both directions (H from 2 calls of U, U from 3 calls of H) give bounded ratios v(U)/v(H). The price
  operator T(v)_K = min over the rules for K of (call weight)/(target weight) then has one cone spectral radius,
  and its threshold in theta is the same for both types (`calib.py` Part B: 0.585 for both). One kernel's saving
  is every kernel's saving.

## 2. Lemma 2 (the group-ring dictionary) [proved from standard theorems]

Let G be the finite group generated by the engine's letters on m address bits (directions, Paulis, affine maps,
Clifford phases). Let rho_f be its f-fold action g -> g^{(x)f}.

**Statement.** For a fixed word, the following are equivalent:
1. it is right for every f;
2. it is right for f = 0, 1, ..., r - 1, where r is the number of distinct values of the trace character of G;
3. its identity holds in Mat(C[G]);
4. its identity holds in every irreducible representation of G, with the same scalar mixing matrices.

**Proof.**
- (1) => (2) is trivial.
- (2) => (4): the trace character chi of the defining representation is faithful. By the Burnside-Brauer theorem,
  every irreducible representation of G is a constituent of some chi^f with f <= r - 1. The character of rho_f is
  chi^f.
- (4) <=> (3): by Wedderburn, C[G] is the product over the irreducible representations rho of M_{d_rho}(C).
- (3) => (1): rho_f is a representation of C[G].

**Consistency checks.**
- This is w-engine's P1 ("right for every f" = "right as a formal sum"). The bound r explains why a few
  consecutive values of f sufficed there.
- It is ex-dimension's remark that the everything-free game collapses at every fixed f: one f sees only the
  constituents of chi^f, never all of them.

**Reading.**
- A unit is one choice of scalar mixing that factorises the kernel in every irreducible representation at once.
  Each representation alone is a stand-in, and the stand-ins are easy.
- The content of a unit is that the stand-ins share their scalars.
- The certificates of Theorem R live on the Pauli-blind representations (reps Theorem 1), which factor through the
  frame geometry. A construction must also satisfy the faithful ones, f = 1 among them.
- The existence question becomes "a nonempty intersection over all rho of the varieties V_rho, with shared
  coordinates". The general existence theorems for such intersections are of two kinds:
  - dimension counts, which fail here (ex-dimension, Lemmas B to D);
  - the projective dimension theorem. It gives points at infinity, that is, border (approximate) units, which
    count only in formally polynomial form (meta-border Theorem X).

  So no algebraic-geometry existence theorem applies naturally. For one fixed shape the question is decidable:
  1 is not in the ideal (Nullstellensatz). That is the exact-algebra and mod-p ground of ex2-modp.

## 3. Proposition 3 (OpenAI's reflection pair at theta < 1, under a transfer hypothesis)

**Hypothesis H_theta** [hypothesis, unproved]:
- Every law OpenAI proves for the 1-price (Theorem 4.1 laws, corners, feedback, specialization, scalar pivots,
  one-port chains) holds for some theta-price.
- Wherever their proofs use the tensor law, use homogeneity instead: price(X^{(x)k}) / |X^{(x)k}| =
  k^theta price(X) / |X|.
- A matrix that their price-preserving exchanges turn into a tensor power has that power's price.

**Statement** [proved from H_theta; computed]. Put pi = price / width.
- The lower bound becomes pi(R) >= (2d)^theta / 2 - 4. At theta = 1 it is d - 4. It uses:
  - T_d ~ U^{(x)d} by exchanges, so pi(T_d) = d^theta / 2;
  - pi(Gamma) = pi(T_d);
  - M = Gamma (x) Gamma, so pi(M) = 2^theta pi(Gamma) = (2d)^theta / 2;
  - then their chain dw <= p(R) + 4w, which uses only direct sums, symmetry, the pivot inequality and pair maps.
- The upper bound becomes pi(R) <= beta(theta)(d/2 + 1) + 2, with beta(theta) = 3 2^(-theta).
  - Their dense-pair bound 3/2 comes from p(H (x) H) = 4 p(H) and a 4-coordinate map of price at most 6.
  - Under homogeneity, p(H (x) H) = 2^(1+theta) p(H), so beta = 6 / 2^(1+theta).
- A contradiction at size d holds iff (2d)^theta / 2 - 4 > beta(theta)(d/2 + 1) + 2.

**Calibration.** At theta = 1 the first such d is 31. The paper says the bounds are incompatible when d > 30.

**Results** (`out_calib.txt`, Part C) [computed]:
- No d gives a contradiction for theta <= theta_c = 0.9639. So this pair could refute theta-prices, and hence
  prove a saving exists, only up to a saving of 1 - theta_c = 0.036, even granting H_theta.
- At theta = 0.9 the best d (about 120) falls short by 7.47 (normalised units).
- Refuting theta = 0.9 needs a dense-pair price below 1.0388 instead of 1.6077. Equivalently, the cheap side must
  be about 35% cheaper.
- For theta = 0.95 the pair price must be below 1.408.

**Where H_theta is used** [read]. OpenAI's own proofs use the tensor law in six places:
- the one-coordinate intersection lemma (`07-pivots.tex` line 106), on which the scalar pivot law, and so every
  exchange in the price of T_m, is built;
- the corner theorem (`05-corners.tex` line 448);
- feedback and specialization (`06-feedback.tex` lines 384 and 463);
- the dense-pair bound (`07-pivots.tex` line 285);
- the price of M (`08-reflection.tex`);
- the packing lemma's tensor sectors (`04-packing-prices.tex`).

Each of these divides an additive error by a growing count of tensor axes. Under homogeneity of degree theta < 1,
that count enters as k^theta, not k. Whether each proof survives is the open question. So Proposition 3 is a
ceiling for this mechanism, not a bound on s*.

**Lemma 3a (the first link already breaks)** [proved by substitution]. Run OpenAI's proof of the
one-coordinate lemma at theta < 1, with two calls to the same Y of width q that share one coordinate. Their
inequality is p(Y (x) Y) <= 2(q - 1) p(Y) + p(P). Homogeneity gives p(Y (x) Y) = q 2^theta p(Y). So the proof
yields only

    p(P) >= (q (2^theta - 2) + 2) p(Y),

instead of p(P) = 2 p(Y). At theta = 0.9 this is 1.73 p(Y) for q = 2, and it says nothing once
q >= 2 / (2 - 2^theta), which is 15 at theta = 0.9.

Consequence. Scalar pivot invariance is derived from this lemma through chains of copies. Its proof therefore
gives only an inequality with a constant-factor loss at theta < 1. The price of T_d needs 2^d - 1 such pivots
(2^m at the step from T_m to T_(m+1)). Unless a new proof makes each pivot exact at theta < 1, the lower bound (2d)^theta / 2 is not available at all.
H_theta is therefore not a plausible routine extension: at least its pivot part needs a new idea.

**A search for a cheaper side at small d** [computed, `sparse/`]. The upper bound counts the supported off-diagonal
pairs h of a real symmetric A with ker A = E, where E is the positive eigenspace of M. OpenAI's
A0 = (J (x) J)(Gamma (x) I - I (x) Gamma) has h = (d/2 + 1) w.
- Reweighted L1 steered towards A0 returns exactly the support of A0 at d = 1 and 2 (w = 16, 64).
- With unsteered normalisations, L1 collapses to rank-1 matrices. The rank condition is not convex, so L1 cannot
  decide the question.
- Greedy removal from A0's support while keeping the generic rank at w/2 (`greedy.py`): no pair can be removed.
  The support is removal-minimal at d = 1 (24 pairs, 5 random orders) and d = 2 (128 pairs, 2 orders)
  [computed, `sparse/out_greedy.txt`]. A sparser A with a different support is not excluded.

## 4. What a duality proof of the bar would have to supply [argued]

Proposition 1 makes "a saving above 0.1 exists" equivalent to "no 0.9-price exists". A proof in this style would
need two things, the two obstructions:

- **(O1) Laws at theta = 0.9.** Either prove the needed part of H_0.9 (corner, feedback, pivot, dense pair), or
  find replacements. Each of OpenAI's laws is first order: it is exact because the 1-price is additive under
  tensor products, which only theta = 1 prices are.
- **(O2) A gap pair.** This is a family X_d with a law-derived lower bound of order (kappa d)^0.9 w, the price of
  a tensor power, and a law-derived upper bound beta d w + O(w), with beta below about 0.52. OpenAI's pair has
  beta = 0.75 at theta = 1 and 0.80 at theta = 0.9.
  - A constant beta is enough at one d. But the pair must beat the tensor price by the factor (2d)^0.1 plus the
    boundary cost, about 1.94 at d of about 120.
  - The matrices involved have width about 4^120 and are never written down. That is exactly an existence proof
    in a place search cannot reach. The unit it implies is finite but unknown, as in OpenAI's paper.
