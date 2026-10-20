> **Working note, kept as written.** An AI agent wrote this during the theory-of-theories round of 14 October 2026 (the spectrum track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts. One skeptic later upheld Theorem C and Lemmas D and E; it corrected Corollary D1's screen to "a growth exponent of at
> least 0.9 at arbitrarily large blocks" and enlarged the block-6 sample (growth exponent about 0.437). The line "No
> skeptic has read this" is out of date.

# meta-spectrum: proofs

14 October 2026, about 06:20 to 08:30 EDT. One agent (meta-spectrum). No skeptic has read this. Nothing posted, no Lean.
Labels: [proved] (full proof here), [computed] (program and output in this folder), [argued], [guess],
[read] (read today at the place cited). Everything is "given Theorem R" (`wht35/wider/w-bounds.md` section 2).

## 1. The two structures

**Tasks and oracles.** An *array of width w* is a register of 2^(w f) entries, read as f digit columns of w bits.
A *task* [r] on such an array (r ≤ w) means: the array is owed H on r chosen bits of every column, fibrewise (the
kernel of block size r on those bits; the other w − r bits are left as they are). An *oracle call* [r] applies that
same operator once to one array of width ≥ r. *Free steps*: rational (Q(i)) mixing of arrays of equal width, and the
free letters of the game acting on one array, the same map in every column (engine: Paulis, affine address maps, the
engine's phases; widened: every address permutation, every Z4 phase pattern). A *program* must be correct for every
f (equivalently: a formal identity in the group algebra, w-bounds section 1).

**Per-block monoid A_m** [definition]. Elements: finite multisets of group elements g of Γ_m (an array of width m owed
g), modulo Q(i)^x scalars. Sum: disjoint union. Preorder: x ≥ y if a program using each element of x once as a paid
letter serves the tasks of y. This is a preordered commutative monoid (Fritz 2017, the setting of "Resource
convertibility and ordered commutative monoids").

**Composition semiring S** [definition]. Elements: finite formal sums Σ W_i [m_i] (W_i arrays owed the full kernel at
block m_i, each array of width m_i). Sum: disjoint union. Product: [m][k] := [mk], extended bilinearly. Preorder:
x ≥ y if a program using each oracle of x once serves the tasks of y. 1 = [1] (one direction), 0 = nothing.

**Lemma A (axioms)** [proved, except where marked].
- (A1) S is a commutative semiring; ≥ is a preorder; x ≥ y implies x + z ≥ y + z. *Proof.* Run programs side by side.
- (A2) x ≥ y implies xz ≥ yz. *Proof.* By distributivity it suffices to take z = [k]. A program P for x ≥ y is correct
  for every f, in particular at column count kf. Regroup k consecutive columns of width w into one column of width kw:
  an oracle [r] at count kf on a width-w array is the oracle [rk] at count f on a width-kw array; a task [m] at count kf
  is the task [mk] at count f; P's free letters, acting identically on the k sub-columns, are free letters of width kw
  (widened game: clearly; engine game: block-diagonal affine maps and Paulis are engine letters; that the engine's phase
  letters stay phase letters is [argued]). So P is a program for x[k] ≥ y[k]. ∎ (This is where "correct for every f"
  is used, hypothesis H3c of wave 17.)
- (A3) 1 ≥ 0, and N → S is an order embedding: n[1] ≥ n'[1] iff n ≥ n'. *Proof.* "If": drop calls. "Only if": take
  any Theorem-R certificate at width 1 with f(D) > 0 (engine game: the permutation module on frames, Y_0 = the line of
  the frame L_0, which every free letter fixes: f(g) = 1 if g moves L_0, else 0; widened game: a hook certificate). Theorem R gives n f(D) ≥ n' f(D). ∎
- (A4) Upper boundedness: [m] ≤ m[1] (the plain schedule K = dir e_1 ... dir e_m). ∎
- (A5) Lower boundedness *inside one width*: on an array of width w ≥ m, one direction costs at most 2 kernel calls:
      dir_e1 = S_1 · K · S_1 · K^-1 · S_1^3   exactly,   K^-1 = K · X_all,
  with S_1 = diag(1, i) on bit 1 and X_all the Pauli on all m bits. [computed exactly in Q(i), m = 1, 2, 3:
  `archimedes.py`, `out_archimedes.txt`; proved for every m: K S_1 K^-1 = (√X S √X^-1) ⊗ I on bits 2..m, so the
  identity is the m = 1 identity on bit 1 tensored with the identity.]
- (A6) Lower boundedness *across widths fails*: no program serves a task [1] on an array of width 1 from oracles [m],
  m ≥ 2, because an oracle of block m cannot act on a narrower array. Padding the array (x ↦ x ⊗ δ_0 on m − 1 extra bits
  per column) would make the identity of A5 usable, but it multiplies the array by 2^((m-1) f), which is not free in
  OpenAI's cost model (work per entry, wave 17 section 3). [argued: this is a statement about the cost model, not a
  theorem.]
- Consequence: S is **not Strassen preordered** in the sense of Definition 2.4 of Fritz, arXiv:1810.08667 (which asks
  for ℓ x ≥ 1 and ℓ ≥ x for every nonzero x), and it is **not of polynomial growth** (Definition 2.8 there asks for
  p(u) x ≥ 1), so neither Strassen's theorem (Theorem 2.6 there, due to Strassen 1988 and Zuiddam) nor Fritz's
  Theorem 2.12 applies to S as stated [read: Fritz 1810.08667, pp. 5-8]. Only (A6) fails.

## 2. Every spectral point of S is a power, or is degenerate

**Theorem C** [proved]. Let φ: S → R_+ be a monotone semiring homomorphism (φ([1]) = 1). Put t_n = φ([n]). Then
either t_n = n^θ for all n ≥ 1, for one θ ∈ [0, 1], or t_n ≤ 1 for all n.

*Proof.* Two inequalities of S: (i) [a + b] ≤ [a] + [b] (Kronecker: H on a + b bits of a column = H on a bits, then H on
the other b bits; one oracle of each); iterating, [Σ x_i] ≤ Σ [x_i]. (ii) [n] ≤ n[1] (A4), so 0 ≤ t_n ≤ n.
Fix a ≥ 2, n ≥ 2, k ≥ 1 and write n^k = Σ_{i=0..L} d_i a^i in base a (0 ≤ d_i ≤ a − 1, L = ⌊k log_a n⌋). By (i),
[n^k] ≤ Σ_i [d_i a^i] = Σ_i [d_i][a]^i (terms with d_i = 0 omitted), so, φ being a monotone homomorphism,
t_n^k ≤ Σ_i t_{d_i} t_a^i ≤ (a − 1) Σ_{i=0..L} t_a^i.
If t_a ≥ 1 this is ≤ (a − 1)(L + 1) t_a^L ≤ (a − 1)(k log_a n + 1) n^(k log_a t_a); taking k-th roots and k → ∞,
t_n ≤ n^(log_a t_a). If t_a < 1 it is ≤ (a − 1)/(1 − t_a) for every k, so t_n ≤ 1.
So if some t_a < 1, all t_n ≤ 1. Otherwise every t_n ≥ 1 and log_n t_n ≤ log_a t_a for all a, n ≥ 2, so
θ := log_a t_a does not depend on a, t_n = n^θ, and θ ∈ [0, 1] by (ii). ∎

**Both kinds occur** [proved].
- θ = 0 (t_n = 1, "count the arrays"): monotone because arrays of different widths never mix, and inside one width w
  the frame-permutation certificate of (A3) has f = 1 on every oracle and every task, so Theorem R gives
  (oracle calls on width-w arrays) ≥ (tasks of width w). [proved in the engine game, where free letters permute frames.
  Widened game: any certificate whose f-value on a rank-r block is at most its demand does the same (true for the hook
  values computed in b56-compute); in general [argued].]
- The width point t_1 = 1, t_n = 0 for n ≥ 2: monotone by the same certificate at width 1 (oracles of block ≥ 2 cannot
  act on width-1 arrays). It is exactly the obstruction to (A6).

**Corollary C1 (one exponent).** Put θ* = sup{θ : n ↦ n^θ is monotone on S}. A unit at block m serving W arrays with
oracle calls [r_1], ..., [r_M] (a block of rank r is the oracle [r]) is the inequality W[m] ≤ Σ_i [r_i]; n ↦ n^θ is
monotone on it iff Σ_i (r_i/m)^θ ≥ W, i.e. iff θ ≤ z_unit (the root of Σ (r_i/m)^z = W, which is the engine's unit
condition, `w-engine.md` 4.1). Hence θ* ≤ z_unit for every unit at every block size [proved]. That conversely every
violated inequality of S yields units with z close to θ* is the analogue of Schönhage's asymptotic sum inequality
(recursion with mixed sizes; Akra-Bazzi) [argued, standard, from memory]. So the spectrum has no other content than
the exponent: there are no further spectral points to discover in S, and the saving is one number s* = 1 − θ*, the
same whether one uses block sizes 2^N, 3^N or anything else.

## 3. What a per-block certificate says about the exponent: its θ-shadow

**Lemma D** [proved]. Let ψ be an additive monotone of A_m that depends on a block only through its rank (every
Theorem-R certificate is one: ψ_r = f(dir e_1 ... dir e_r), ψ_m = f(K) = the demand; the double-coset reduction of
pf-duality Lemma 1.1). Put
      θ(ψ) = min over 1 ≤ j < m of ln(ψ_m/ψ_j) / ln(m/j).
Then every unit at block m has z_unit ≥ θ(ψ), i.e. saving ≤ 1 − θ(ψ) at engine prices.
*Proof.* Theorem R: Σ_i ψ_{r_i} ≥ W ψ_m. By definition ψ_r/ψ_m ≤ (r/m)^θ(ψ). So Σ_i (r_i/m)^θ(ψ) ≥ W. ∎
(This is b56-compute's `pricecap.py` formula; the point here is that it is exactly the restriction of the power point
n ↦ n^θ to block m.)

**Corollary D1.** If for every block m some certificate ψ^(m) has θ(ψ^(m)) ≥ θ_0, then s* ≤ 1 − θ_0. A family that
does this is the only way a certificate method can bound s*. Hooks and coordinate certificates have ψ_m/ψ_1 below
4(n−1)/n < 4 (b56-note), so θ(ψ) ≤ log_m 4 → 0: they can never bound s* [proved, given the b56 ceilings, which are
proved for hooks at every m and for coordinate certificates at m = 3..10].

**Lemma E (mixtures)** [proved; Farkas]. Direct sums of certificates are certificates and their profile is the sum
(pf-ladder 2.1). For a finite set C of profiles, the shapes (n_1, ..., n_{m-1}, W) ≥ 0 that no c ∈ C excludes form the
cone {Σ_r n_r ψ^c_r ≥ W ψ^c_m for all c}. The largest z with Σ_r n_r (r/m)^z ≥ W on that whole cone equals
max over convex weights α of θ(Σ_c α_c ψ^c / ψ^c_m). *Proof.* For fixed z, Farkas' lemma: the inequality is implied by
the system and n ≥ 0 iff its coefficient vector is a nonnegative combination of theirs, i.e. iff some mixture has
ψ_r/ψ_m ≤ (r/m)^z for all r. ∎ So "choose the best certificate for each unit" never beats the best single mixture.

## 4. The only θ-consistent point we know is 2-adic

**Proposition G** [proved; w-bounds section 4 and pf-duality Corollary C(a), assembled]. Restrict S to programs whose
mixing coefficients have odd denominators. Then n ↦ n (θ = 1) is monotone: no such unit saves, at any block, in any
composition.
*Proof.* Theorem R over GF(2) with V = F_2^(2m) (Pauli labels, Γ acting by conjugation), Y_0 = L_0. A rank-r block is a
product of r commuting transvections along (e_i|0) and moves L_0 to distance r; K moves L_0 to the opposite frame
{(b|b)}; so f(rank-r block) = r and f(K) = m at every m, and Σ r_i ≥ W m. Per width in a mixed program. ∎
Over Q this point is not monotone: positive WHT savings exist (about 1.2e-3; HANDOFF line 338 [read]), so θ* < 1.
Lifts (w-bounds 5, pf-duality Lemma I2, b56-methods) are exactly the characteristic-0 deformations of this point:
δ ≥ distance, with equality on the pair frames of the (2m+1)-point model.

**No single characteristic-0 certificate attains θ = 1 at any block m ≥ 3** [argued for m ≥ 4; proved at m = 3,
w-bounds Theorem F]: a perfect family (δ ≤ c on edges, δ ≥ mc on opposite pairs) is tight on every geodesic, hence
perfect on every length-3 sub-interval (isomorphic to the block-3 interval, pf-duality Lemma I1), where the Fano
argument of Theorem F applies. The step "every pair at distance 3 lies on a geodesic between opposite frames" is a
standard property of dual polar graphs that I did not re-prove.
