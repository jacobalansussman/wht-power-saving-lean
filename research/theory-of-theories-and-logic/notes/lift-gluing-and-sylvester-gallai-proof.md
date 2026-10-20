> **Working note, kept as written.** An AI agent wrote this during the second abstract wave of 14 October 2026 (the lifts track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts.
>
> **Refuted part.** Proposition C (every counting argument of this kind is capped at 1/3 per edge, so the bar's
> threshold overtakes it from block 18) is false as written; both skeptics found it. What survives: such counting is
> capped by some constant, which is at least 0.38 per edge for a sharper version and at most about 1 for the whole
> class. Everything else used in the README (Lemma E3, Theorem G3, b_4 >= 18, b_5 >= 141, Theorem SG) was upheld by
> both skeptics.

# PROOF: a2-sg-lifts (14 October 2026)

Labels: [proved] (proof here), [computed] (program and output named), [read] (source read), [argued]. No skeptic has
read this file. Notation and Facts 0, T, U, Lemma LG, Lemma S, Corollary M, Lemma F2, Theorem B3, Lemma Dec, Lemma P
are those of `../../meta/c2-lifts/PROOF.md` (cited below as c2).

Additional notation for the block-3 interval [0, F] (c2 section 3): p ranges over the 7 points, t over the 7 lines;
end edges are (0,p) ("bottom") and (t,F) ("top"), middle edges are the 21 flags (p,t). p_0 and p_F are the
projections of Q^6 = Y_0 ⊕ Y_F (c2 Fact 0: Y_0 ∩ Y_F = 0). As in c2 Theorem B3:
A_p = Y_p ∩ Y_0, C_p = Y_p ∩ Y_F, B_t = Y_t ∩ Y_0, D_t = Y_t ∩ Y_F.

## 0. Three small facts used below

**(I) Triple intersections.** For frames U, V, W, dim(Y_U ∩ Y_V ∩ Y_W) ≤ dim(U ∩ V ∩ W). *Proof:* a primitive integer
vector of the intersection lies in each saturated lattice Λ_U, Λ_V, Λ_W, so its (nonzero) reduction lies in U ∩ V ∩ W;
a saturated lattice of rank r reduces to an r-dimensional space. ∎ Consequences in the 3-interval: the B_t (when
1-dimensional) are pairwise distinct, the C_p are pairwise distinct, and D_t ≠ D_{t'} for t ≠ t' (N_t ∩ N_{t'} ∩ M_0
= (t∩t'|0) is 1-dimensional).

**(II) Fano minus one line.** Let points X_q (q ∈ PG(2,2)) and lines Λ_ℓ (ℓ ≠ ℓ_0) be given in P^2(Q) with X_q ∈ Λ_ℓ
when q ∈ ℓ and X_q ∉ Λ_ℓ when q ∉ ℓ (all ℓ ≠ ℓ_0). Then the three points X_q, q ∈ ℓ_0, are not collinear.
*Proof:* the coordinate argument of c2 Lemma F2(ii), with a basis e_1, e_2, e_3 chosen so that ℓ_0 = {e_1+e_2,
e_1+e_3, e_2+e_3}; it uses only lines ≠ ℓ_0 and ends with det[[1,1,0],[1,0,1],[0,1,1]] = −2 ≠ 0. ∎
The same holds with the roles of points and lines exchanged (PG(2,2) is self-dual).

**(III) The bottom-defect map.** Suppose dim A_{p0} = 1 (so e(0,p0) = 1) and dim C_{p0} = 1. Then ℓ_0 := p_0(Y_{p0}) is
2-dimensional, A_{p0} ⊂ ℓ_0, and there is a linear map φ: ℓ_0 → Y_F/C_{p0} with
Y_{p0} = {y + w : y ∈ ℓ_0, w ∈ φ(y)}, ker φ = A_{p0}, rank φ = 1. Let t ∋ p0 be a line with (p0,t) clean and
Y_t = B_t ⊕ D_t (dim B_t = 1, dim D_t = 2). Then C_{p0} ⊂ D_t, B_t ⊂ ℓ_0, and either B_t = A_{p0} or
im φ = D_t/C_{p0}. *Proof:* Y_{p0} ∩ Y_t has dimension 2; p_0 maps it into p_0(Y_t) = B_t with kernel
Y_{p0} ∩ Y_t ∩ Y_F ⊂ C_{p0}; so both are equalities: C_{p0} ⊂ Y_t ∩ Y_F = D_t, and some u = b_t + d ∈ Y_{p0} with
b_t spanning B_t and d ∈ D_t. Then B_t ⊂ ℓ_0 and φ(b_t) = d + C_{p0} ∈ D_t/C_{p0}; if b_t ∉ ker φ = A_{p0} this
spans im φ, and D_t/C_{p0} is 1-dimensional. ∎
**Corollary (III′).** At most two lines t ∋ p0 satisfy the hypotheses of (III), and if two do, exactly one of them has
B_t = A_{p0}. *Proof:* two lines with B_t ≠ A_{p0} would give D_t/C_{p0} = D_{t'}/C_{p0}, hence D_t = D_{t'}, against
(I); and the B_t are distinct, so at most one equals A_{p0}. With three lines, two have B_t ≠ A_{p0}. ∎

## 1. Lemma E3 (end edges cost extra) [proved]

**Lemma E3.** Every lift of the 3-interval in which some end edge has positive excess has total excess E ≥ 3.

*Proof.* The integer matrix swapping x and z is symplectic and reverses the interval (c2 F2(i)), so take the bad end
edge to be (0,p0). Suppose E ≤ 2. By c2 Lemma F2(ii) some edge (p,t) or (t,F) is bad, so e(0,p0) = 1 and exactly one
other edge X ∈ {(p,t), (t,F)} is bad, with e(X) = 1. Everything else is clean. As in c2 Theorem B3 (Fact T along clean
chains): A_p is 2-dimensional for p ≠ p0, A_{p0} is 1-dimensional; B_t is 1-dimensional for every t and B_t ⊂ A_p
whenever p ∈ t, p ≠ p0 and (p,t) is clean (every line has such a p: at most p0 and one endpoint of X are excluded);
D_t is 2-dimensional when (t,F) is clean; C_p is 1-dimensional for every p (some line t ∋ p has (p,t) and (t,F)
clean) with C_p ⊂ D_t on clean chains. Non-adjacent p, t (p ∉ t) are opposite: Y_p ∩ Y_t = 0.

*Case 1: X = (t1,F).* All flags are clean, so Y_t = B_t ⊕ D_t for t ≠ t1 and Y_p = A_p ⊕ C_p for p ≠ p0. Every line
t ∋ p0 other than t1 satisfies the hypotheses of (III). By (III′), p0 ∈ t1; write t2, t3 for the other two lines
through p0, with B_{t2} = A_{p0}, so ℓ_0 = B_{t2} + B_{t3}.
The lower configuration (points B_t, lines A_p for p ≠ p0) has every incidence and non-incidence of PG(2,2)
except the missing line p0, so by (II) B_{t1}, B_{t2}, B_{t3} are not collinear: B_{t1} ⊄ ℓ_0.
For p ∈ t1 \ {p0}: Y_p ∩ Y_{t1} is 2-dimensional and contains B_{t1} (⊂ A_p), and p_F maps it into p_F(Y_p) = C_p with
kernel A_p ∩ B_{t1} = B_{t1}; so C_p ⊂ p_F(Y_{t1}), which is 2-dimensional (= 3 − dim B_{t1}). Hence
p_F(Y_{t1}) = C_{p'} + C_{p''} (the two points of t1 other than p0, distinct by (I)). The upper configuration (points
C_p, lines D_t for t ≠ t1) misses only the line t1, so by (II) C_{p0} ⊄ C_{p'} + C_{p''}; hence C_{p0} ⊄ Y_{t1}.
Now take W = Y_{p0} ∩ Y_{t1}. W ∩ Y_F ⊂ C_{p0} ∩ Y_{t1} = 0 and W ∩ Y_0 = A_{p0} ∩ B_{t1} = B_{t2} ∩ B_{t1} = 0. So p_0
is injective on W with image in ℓ_0 ∩ p_0(Y_{t1}); p_0(Y_{t1}) is 2-dimensional, contains B_{t1} ⊄ ℓ_0, so
ℓ_0 ∩ p_0(Y_{t1}) is 1-dimensional and dim W ≤ 1: the flag (p0,t1) is not clean. Contradiction.

*Case 2: X = (p1,t1), a flag.* All (t,F) are clean, so every Y_t = B_t ⊕ D_t, and Y_p = A_p ⊕ C_p for p ≠ p0.
(2a) p1 ≠ p0: all three lines through p0 satisfy (III), against (III′).
(2b) p1 = p0: t2, t3 (the lines through p0 other than t1) satisfy (III); by (III′) B_{t2} = A_{p0} and
ℓ_0 = B_{t2} + B_{t3}. All flags (p,t) with p ≠ p0 are clean, so the lower configuration misses only the line p0 and
B_{t1} ⊄ ℓ_0 (II). The upper configuration has every incidence except (p0,t1), so C_{p0} ⊄ D_{t1} (otherwise it
would realise PG(2,2) over Q, which c2 Lemma F2's argument excludes). For u ∈ Y_{p0} ∩ Y_{t1}: p_0(u) ∈ ℓ_0 ∩ B_{t1}
= 0, so u ∈ Y_{p0} ∩ Y_F = C_{p0} and u ∈ Y_{t1} ∩ Y_F = D_{t1}; so u = 0, dim(Y_{p0} ∩ Y_{t1}) = 0 and
e(p0,t1) = 2, giving E = 3. Contradiction.

There are no other cases (X lies in the upper two levels). ∎

**Remark (tightness) [computed].** Case 2b with e(p0,t1) = 2 can be completed to a lift with E = 3: lower
configuration = the complete quadrangle with its three diagonal points (PG(2,2) minus a line), upper configuration =
the same with one diagonal point moved off one line (PG(2,2) minus a flag), and
Y_{p0} = A_{p0} ⊕ C_{p0} ⊕ ⟨b_{t3} + 2w⟩ with A_{p0} = B_{t2} and w an integer vector of D_{t3} outside C_{p0}; all other
Y's decomposed. An explicit instance (both configurations from 0/1 lifts, one transported by an integer matrix so that
the failing flag sits at p0 = 111) has exactly E = 3, with e(0,111) = 1 and e(111,{001,110,111}) = 2, every Y reducing
to a basis of its frame (exact fractions: sgl-e3/e3_tight.py, out_e3_tight.txt). So the "3" of Lemma E3 is best
possible. (The annealer endedge3.c only found 5 with a bad end edge; it is weak.)

## 2. Lemma E21 (one more forbidden pattern) [proved]

**Lemma E21.** No lift of the 3-interval has excess 2 on the bottom end edges, 1 on the top end edges and 0 on the
middle flags (nor the mirror pattern).

*Proof.* All flags clean, exactly one top edge (t1,F) bad with e = 1. Bottom: either (i) e(0,p0) = e(0,p1) = 1, or
(ii) e(0,p0) = 2.
(i) A_p is 2-dimensional for p ∉ {p0,p1}; every line has a point outside {p0,p1}, so every B_t is 1-dimensional and
Y_t = B_t ⊕ D_t for t ≠ t1; every C_p is 1-dimensional. Through p0 the lines other than t1 satisfy (III), so by (III′)
p0 ∈ t1, and likewise p1 ∈ t1: t1 = {p0, p1, p2}. Let t3 ∋ p0 and t' ∋ p1 be the lines (≠ t1) with B ≠ A; then
p_F(Y_{p0}) = C_{p0} + D_{t3} = D_{t3} (im φ = D_{t3}/C_{p0}) and p_F(Y_{p1}) = D_{t'}. B_{t1} is 1-dimensional (via p2)
and differs from A_{p0} = B_{t2}, so p_F is injective on Y_{p0} ∩ Y_{t1}, which is 2-dimensional (clean flag): hence
p_F(Y_{t1}) (2-dimensional) = D_{t3}. In the same way p_F(Y_{t1}) = D_{t'}. So D_{t3} = D_{t'}, t3 = t' (by (I)), a line
through p0 and p1 other than t1: impossible.
(ii) A_{p0} = 0. For t ∋ p0, t ≠ t1, the proof of (III) still gives B_t ⊂ ℓ_0 := p_0(Y_{p0}) (2-dimensional) and
C_{p0} ⊂ D_t. The lower configuration misses only the line p0, so the three B_t (t ∋ p0) are not collinear (II):
hence p0 ∈ t1, ℓ_0 = B_{t2} + B_{t3}, B_{t1} ⊄ ℓ_0. W = Y_{p0} ∩ Y_{t1} is 2-dimensional; p_0(W) ⊂ ℓ_0 ∩ p_0(Y_{t1}) is at
most 1-dimensional, so W ∩ Y_F ≠ 0, i.e. C_{p0} ⊂ Y_{t1} ∩ Y_F, which is 1-dimensional: Y_{t1} ∩ Y_F = C_{p0}. As in
Case 1 of E3, p_F(Y_{t1}) = C_{p'} + C_{p''} (p', p'' the other points of t1), and Y_{t1} ∩ Y_F ⊂ p_F(Y_{t1}); so
C_{p0} ∈ C_{p'} + C_{p''}, against (II) for the upper configuration (which misses only t1). ∎

## 3. Theorem G3 (the gluing gap at block 3 is real): every whole-graph block-3 lift has E ≥ 68 [proved]

*Proof.* Restrict a whole-graph lift to each of the 4320 intervals (c2 Lemma LG: the restriction is an interval lift
with the same edge excesses; no slicing). Each of the 945 edges lies in 160 intervals: as an end edge in 64 and as a
middle flag in 96 [computed: glue/glue3_ilp.py prints "end-edge incidences per edge: [64]"; 4320·14/945 = 64]. For an
interval I let k_b, k_t be its excess on bottom and top end edges (k = k_b + k_t) and μ_I its excess on middle flags.
The local constraints give μ_I ≥ g(k_b, k_t) with
- g(0,0) = 2 (Theorem B3), g(1,0) = 2 (Lemma E3: E_I ≥ 3);
- g(2,0) = 1 and g(k_b,0) ≥ 1 for all k_b (F2(ii): the upper two levels need excess, and the top is clean);
- g(1,1) = 1 (Lemma E3); g(2,1) ≥ 1 (Lemma E21); symmetric in k_b, k_t;
- g ≥ 0 otherwise (k ≥ 4 with both sides positive).
Every case satisfies g ≥ 2 − k/2. Summing over intervals: 96E = Σ_I μ_I ≥ 2·4320 − (1/2)·Σ_I k_I = 8640 − 32E.
So E ≥ 8640/128 = 67.5, i.e. E ≥ 68. ∎

**Consequences.**
- λ_3 ≥ 68/945 = 0.0720 (was 2/35 = 0.0571), and by c2 Corollary M, λ_m ≥ 68/945 for every m ≥ 3. (At block 4 the
  interval bound b_4 ≥ 18 of Corollary G3′ is stronger: λ_m ≥ 18/240 = 0.075 for every m ≥ 4, so every whole-graph
  block-4 lift has E ≥ 2,582.)
- Every block-3 lift proves at most 3/(1 + 68/945) = 2835/1013 = 2.7986 moves per array (was 105/37 = 2.838); growth
  exponent ≤ ln(2.7986)/ln 3 = 0.9367.
- **No lift can exclude 5 in 14 at block 3** (that needs E ≤ 67). This settles c2-lifts' guess P3 (gluing gap
  g_3 ≥ 68/54 = 1.26).
- Čech reading: on the cover of the frame graph by intervals, a family of locally optimal interval lifts (E_I = 2)
  never glues: by Lemma E3 each such local optimum vanishes on its 14 end edges, every edge is an end edge of 64
  intervals, so a compatible family would vanish everywhere, against b_3 > 0. Theorem G3 is the quantitative form.

**Corollary G3′ (the gluing row survives slicing) [proved].** Let [a,b] be a 3-dimensional sub-interval of a block-m
interval (or any 3-interval of the block-m frame graph), μ its excess on its 21 middle flags and k its excess on its 14
end edges, all measured in the original lift. Then μ + k/2 ≥ 2. *Proof:* slicing (c2 Corollary S1 / Lemma S) gives a
lift of the block-3 interval with excesses μ′ ≤ μ, k′ ≤ k edge by edge; the table above gives μ′ ≥ g(k′_b, k′_t) ≥
2 − k′/2 ≥ 2 − k/2. ∎ With the Fano rows of c2 section 6 this gives, by integer programming,
**b_4 ≥ 18 for every lift** (LP 17.14, ILP optimal; sgl-cover/cover2_ilp.py, out_cover2_m4.txt), up from 14.
On the whole block-m graph the same row, summed over all 3-intervals (each edge is an end edge in 14/35 of those
through it), gives only λ_m ≥ 1/14 (the uniform LP optimum), weaker than 68/945 via c2 Corollary M.

## 4. The Sylvester–Gallai transfer for decomposed lifts [proved, asymptotics via DSW [read]]

An *A-part* on F_2^m (c2 Lemma Dec) gives every subspace s a rational subspace A_s ⊂ Q^m, dim A_s = dim s, with an
integer basis reducing to a basis of s; f(s,t) = dim s − dim(A_s ∩ A_t) on covering pairs s ⋖ t; E_A = Σ f. A
decomposed interval lift has E = E_A + E_{A'} for two A-parts (c2 Lemma Dec; the second is s ↦ B_s^⊥, which has an
integer basis reducing to a basis of s by the annihilator argument of c2 section 4), and b_m^dec = 2 b_m^A.
Let n_k be the least number of bad lines (three points not linearly dependent) over all configurations of the
2^k − 1 points of PG(k−1,2) in Q^k with rank_Q(S) ≥ rank_2(S) for every set S of points (integer lifts P_x ≡ x mod 2
are such configurations, c2 Fact 0; c2 Lemma N works in this generality). L_k = (2^k−1)(2^{k−1}−1)/3 is the number of
lines and ν_k = n_k/L_k. Known: ν_2 = 0, ν_3 = 1/7, ν_4 = 9/35 (c2 Proposition N4 and the E_A = 18 configuration),
ν_5 ≤ 55/155, ν_6 ≤ 270/651 [computed: goodlines.c, integer lifts].

**Lemma M (monotone).** ν_{k+1} ≥ ν_k. *Proof:* in a configuration for PG(k,2) the 2^k − 1 points of a hyperplane H
have rank_Q ≥ rank_2 on all their subsets; a generic linear map Q^{k+1} → Q^k sends every set of rank r to rank
min(r, k), so it keeps all these rank conditions (rank_2 ≤ k inside H) and keeps every 3-set's rank (k ≥ 3): the image
is a configuration for PG(k−1,2) with the same bad lines, at least n_k of them. Every line lies in the same number of
hyperplanes; average. ∎

**Lemma L (link lemma).** For an A-part and 1 ≤ j ≤ m − 1 let φ_i be the fraction of failing covering pairs at level i
(dim s = i). Then
  φ_j + φ_{j−1} ≥ ν_{m−j+1}/3   and   φ_j + φ_{j+1} ≥ ν_{j+2}/3,
with φ_0 = φ_{m−1} = 0 (A_0 = 0 and A_{F} = Q^m make the end pairs clean).
*Proof.* Fix a of dimension j − 1 and put V = Q^m/A_a with lattice Z^m/Λ_a (Λ_a = A_a ∩ Z^m, saturated), which is
free of rank k = m − j + 1 and reduces onto F_2^m/a. For s ⊃ a of dimension j with A_a ⊂ A_s, Λ_s/Λ_a is a saturated
rank-1 lattice of V reducing to the point s/a; call its line Q_s. For any set S of such s, rank_Q(Σ Q_s) =
dim(Σ A_s) − dim A_a ≥ dim_2(Σ s) − dim a (stack integer bases). For the remaining s (A_a ⊄ A_s) choose Q_s generic.
This is a configuration of PG(k−1,2) with rank_Q ≥ rank_2 (adding generic points keeps the rank conditions), so at
least n_k of its lines are bad. A line t/a
(t ⊃ a, dim t = j+1) all of whose three points s satisfy A_a ⊂ A_s ⊂ A_t is good (its points lie in the 2-space
A_t/A_a). So every bad line has a failing pair (s,t) (level j) or a failing pair (a,s) (level j − 1, and then all
2^{k−1} − 1 lines through s/a are charged to it). Summing over the [m, j−1] choices of a: a failing level-j pair is
counted for the 2^j − 1 hyperplanes a of s, a failing level-(j−1) pair (a,s) for 2^{m−j} − 1 lines, and
(2^j − 1)e_{m,j} = (2^{m−j} − 1)e_{m,j−1} = 3[m, j−1]L_k (e_{m,i} = [m,i](2^{m−i} − 1) pairs at level i). Dividing:
3(φ_j + φ_{j−1}) ≥ ν_k. The second inequality is the first for the dual A-part u ↦ A_{u^0}^⊥ (u = s^0), which has the
same failure on the reversed pair (dim(A_t^⊥ ∩ A_s^⊥) = m − dim(A_s + A_t)), at level m − j − 1. ∎
At j = 1 this is c2 Lemma P; at m = 4 it gives E_A ≥ 9 + 9 = 18 = b_4^A (c2 Theorem D4).

**Theorem SG (asymptotics).** 1 − ν_k ≤ (24 + o(1))/k. Hence, for decomposed lifts,
b_m^dec/e_m ≥ (2/3)κ_m − o(1), where κ_m is the least weight fraction of a set of levels meeting every adjacent pair
(levels 1..m−2; weights e_{m,i}/e_m). κ_m → 1/2 for even m (levels j and m−1−j have equal weight, so odd and even
levels each carry half; the forced end levels 1 and m−2 add o(1)) and κ_m → 0.4716 for odd m [computed: sgl-vc/vc.py, out_vc.txt]. So the asymptotic
SG bound is 1/3 per edge (even m) and 0.314 per edge (odd m).
*Proof.* Dvir–Saraf–Wigderson, arXiv:1211.0330, Corollary 1.11 [read: statement on p. 5, combining BDWY Theorem 1.10
with their Theorem 1.9, affine-dim ≤ 12/δ for δ-SG configurations over C]: n distinct points of C^d with at least δn^2
unordered pairs on special lines contain V′ with |V′| ≥ (δ/6)n that is a δ-SG configuration, hence affine-dim(V′) ≤
12/δ. Put the 2^k − 1 points P_x into an affine chart (a rational functional nonzero on all of them): projective
collinearity becomes affine collinearity and points stay distinct (rank 2 on pairs). The good lines give at least
3(1 − ν)L_k = (1 − ν)n(n−1)/2 pairs on special lines, so δ = (1 − ν)(n−1)/(2n). Then rank_2(V′) ≤ rank_Q(V′) ≤
12/δ + 1, so |V′| ≤ 2^{12/δ+1}, while |V′| ≥ δn/6 with n = 2^k − 1. Taking logarithms, k ≤ 12/δ + log_2(6/δ) + O(1),
i.e. 1 − ν_k ≤ 24/(k − O(log k)). For the lifts: Lemma L says φ_i + φ_{i+1} ≥ c_i/3 for every adjacent pair of levels,
with c_i ≥ min(ν_{m−i}, ν_{i+2}); φ_0 = φ_{m−1} = 0. The levels with i or m − i below εm carry weight o(1), and
elsewhere c_i ≥ 1 − O(1/(εm)). The resulting linear program on a path is a fractional vertex cover, whose optimum is
integral (paths are bipartite): its value is (1/3)κ_m per A-part, up to o(1). ∎
The bound is ineffective at our sizes: 1 − ν_k ≤ 1/2 needs k ≳ 54, i.e. blocks m ≳ 110.

**Proposition C (ceiling of line-level SG certificates) [proved].** The constraints of Lemma L (for any values ν ≤ 1)
are all satisfied by φ_i ≡ 1/6, and their exact optimum for ν ≡ 1 is (1/3)κ_m ≤ 1/6 per A-part. So no argument that
uses only "every bad line in every link costs one failing pair" can prove more than 1/3 per edge for decomposed lifts
(1/6 per A-part). Since m^0.1 − 1 > 1/3 for m ≥ 18, such arguments can never show that decomposed lifts miss
growth exponent 0.9 at large blocks. ∎

## 5. General lifts: the escape [argued]

For a general (non-decomposed) lift the same link argument, run on a sub-interval [a,b] after slicing (c2
Corollary S1), charges a bad link line either to a failing flag or to a *top defect* of the line relative to b
(c2 Lemma X: a(u) ≤ a(w) + e(u,w)), which can be produced by excess on any later level of [a,b]. In per-level
densities x_i this gives only 3x_j + x_{j+1} + … + x_{j+k−2} ≥ ν_k, satisfied by putting excess 1/(k+1) on every level,
or 1 on a single level. So the SG transfer for general lifts is weaker still, and the only general-lift bounds that
survive are the Fano cover (c2 section 6) and the gluing count of section 3.

## 6. Block 5, decomposed lifts: link-by-link integer programs [proved rules; bounds computed]

For an A-part of F_2^5 let F1, F2, F3 be the failure indicators of the covering pairs (point, line), (line, plane),
(plane, hyperplane); E_A ≥ ΣF1 + ΣF2 + ΣF3. Valid rows (all from Lemma L applied to one link at a time, and its dual):
(a) at least n_5 ≥ 40 lines carry a failing point flag (n_5 ≥ 155·9/35 by Lemma M); (b) every plane has such a line
(Fano); (c) every hyperplane contains at least n_4 = 9 such lines; (d) in the link of each point x, at least 9 planes
π ∋ x are killed by some F1[x,ℓ] or F2[ℓ,π] with x ∈ ℓ ⊂ π; (f) in the link of each line, at least one hyperplane is
killed; (e), (g) and the starred versions of (a)–(c) are the duals.
**Flatness rules.** Call a plane flat if its 7 points span a 3-dimensional Q-space.
- NF: at most one plane through a bad line is flat. *Proof:* if π ≠ π′ ⊃ ℓ are flat, both Q-spaces contain the three
  non-collinear points of ℓ, so both equal their 3-dimensional span W; then the points of π ∪ π′, of F_2-rank 4, lie
  in W, against rank_Q ≥ rank_2. ∎ (At block 4 this sharpens c2 Lemma N's R3 from "one non-flat plane" to "two of the
  three".)
- BLOCK: a non-flat plane π has a point x with A_x ⊄ A_π (else all 7 points lie in the 3-space A_π), so for each of the
  three lines ℓ with x ∈ ℓ ⊂ π, F1[x,ℓ] or F2[ℓ,π] holds (A_x ⊂ A_ℓ ⊂ A_π is impossible). ∎ Duals NF*, BLOCK* hold for
  the dual A-part.
**Membership rules.** Put in[x,π] = [A_x ⊂ A_π] for x ∈ π.
- CHAIN: F1[x,ℓ] = F2[ℓ,π] = 0 implies in[x,π] = 1.
- OUT: a non-flat plane has a point with in = 0; a flat plane (points in a 3-space W) has either all 7 points in A_π
  (A_π = W) or at least 4 out (A_π ∩ W is at most 2-dimensional and holds at most the 3 points of one line: any 3 points
  of π not on a line are Q-independent by the rank conditions).
- PAIR: if F2[ℓ,π] = 1, at most one point y of ℓ has F1[y,ℓ] = 0 and in[y,π] = 1 (two such points y ≠ z give
  A_ℓ = A_y + A_z ⊂ A_π).
- GOOD LINE: a line whose three points are collinear has 0 or at least 2 failing point flags (if A_ℓ is not their
  2-space W, A_ℓ ∩ W holds at most one of them).
- and the duals (in*[H,ℓ] = [A_ℓ ⊂ A_H]; dual flatness of a line: the normals of the 7 hyperplanes through it span a
  3-space). ∎ (each line is a one-step consequence of dimension counting.)
Calibration: on the two annealed A-parts (E_A = 290 and 307, re-checked exactly) every row holds; they have 55 and 61 bad
lines and 100 and 98 non-flat planes [computed: sgl-link5/check_apart5.py, calib_flat5.py].
Results [computed, HiGHS on a Colab machine]: without flatness rules LP 80, ILP dual bound 86, feasible point 123 with only 9
middle failures (out_link5_n40.txt). With the flatness rules (out_link5b_n40.txt) and the membership rules as well
(sgl-link5c/out_link5c_n40.txt) the LP stays at 80 and the integer programs made no progress in 20–25 minutes on a
loaded Colab (dual bounds 82 and 80). So **b_5^A ≥ 86** and every
decomposed block-5 interval lift has E ≥ 172 (0.083 per edge); the best found has 580 (0.279 per edge).
