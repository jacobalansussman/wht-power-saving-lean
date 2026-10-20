> **Working note, kept as written.** An AI agent wrote this during the theory-of-theories round of 14 October 2026 (the cohomology track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts. One skeptic later upheld Lemmas 1 to 3, the Proposition, Lemma 5 and Theorem E+, with its own code
> (`programs/clifford-cohomology/indep.py`).

# meta-cohomology: the sign-defect question of ut-direct 7.5 is decided. The class vanishes; the transfer works.

14 October 2026, 06:12-06:50 EDT (times from the out_* files). One agent. Nothing posted, nothing in Lean, nothing outside this folder written.
Labels: [proved] = written out in full here, every step checkable; [computed] = script in this folder, output named;
[argued] = reasoning with a named gap; [read] = taken from the file named; [guess].
Jobs: all tiny (< 0.1 GB, < 2 min); run under a memory-capped job runner.

## 0. Results

1. **The sign defect is identically zero** [proved, s2 Lemma 2; computed on 320 random QCAs, m = 1, 2, 3].
   ut-direct PROOF.md 7.5 says the naive collapse fails to be a Clifford map by (-1)^gamma,
   gamma(v,w) = sum_{d>d'} <M_d v, M_d' w>, and that gamma is alternating. In fact gamma = 0: grouping the sum by the lag
   e = d - d' >= 1, each lag sum is sum_d <M_d v, M_(d-e) w> = delta_(e,0) <v,w> = 0 by the same twisted symplectic
   identities. So every single collapsed rule is already a Clifford map; no quadratic refinement is needed.
2. **The naive collapse is multiplicative** [proved, s2 Lemma 3; computed exactly, 320 random pairs, signs included].
   Pi_0(alpha beta) = Pi_0(alpha) Pi_0(beta) for all Clifford QCAs, where Pi_0 replaces the tensor product over sites by
   the product over sites in increasing order. The proof is the same lag identity applied to the composite.
3. **The cohomology question, stated and answered** [proved, s1 and s2]. Group: Gamma = the bulk symbols, a subgroup of
   the twisted symplectic group Sp_tw(2m, F_2[S, S^-1]) = {M(S): M(S^-1)^T J M(S) = J}. Module: V = F_2^(2m) with Gamma
   acting through pi_1: M -> M(1). Classes in H^2(Gamma; V): [B] (the extension 0 -> V -> {Clifford QCAs mod scalars}
   -> Gamma -> 1, kernel = translation-invariant Pauli conjugations) and pi_1^*[C_m] (pull-back of the one-column
   Clifford extension 0 -> V -> C_m/mu_8 -> Sp(2m,2) -> 1). The question "can the corrections be chosen
   multiplicatively" is exactly "[B] = pi_1^*[C_m]?". Answer: **yes, the difference class is 0**, and Pi_0 is an explicit
   morphism of extensions (identity on V, identity on the constants). The set of such morphisms is a torsor under the
   crossed homomorphisms Z^1(Gamma; V) vanishing on Sp(2m,2).
4. **Why** [computed + argued, s3]: Pi_0 is restriction to the translation-invariant sector of a ring. On a ring of odd
   length L, F_2[S]/(S^L - 1) = F_2 x (the rest), with idempotent e_1 = 1 + S + ... + S^(L-1); the S = 1 summand is an
   orthogonal symplectic summand, the logical Paulis P_u^(x)L generate a copy of the m-qubit Pauli group, every
   Clifford QCA preserves it, and Clifford groups are functorial for orthogonal sums. Exactly: for L = 1 mod 4 the
   restriction equals Pi_0, for L = 3 mod 4 it equals the complex conjugate of Pi_0 [computed, all 317 cases, L = 7..15].
5. **Theorem E+ (hook certificates transfer to locality-preserving Clifford letters)** [proved given Theorem
   R-infinity (ut-direct s2) and the cited certificates, s4]. Let every letter be a Clifford family in Q_loc (s4.1:
   locality-preserving with translation-invariant bulk; contains one-column Cliffords, directions, the kernel, every
   unimodular chain map such as row0 ^= row1 << 1 and its anti-causal mirror, column-uniform CZ/Z4 couplings, cyclic
   column rotations, and every boundary-local Clifford); let every free letter be monomial or boundary-local. Then every
   certificate N in V^(x)t0 that is kept by g^(x)t0 for all monomial one-column Cliffords g gives in this game the same
   bound as in the tensor-power game. In particular: **>= 2W at block 2, >= 14W/5 at block 3, >= 41W/13 at block 4**
   (the hook certificates of d-tof-cert and b4-depth4 are kept by every monomial letter [read]).
   This upgrades ut-direct Theorem E at block 3 from 105W/41 to 14W/5 and adds block 4.
6. **Where the obstruction lives (where a unit could still hide)** [proved for the two named letters; computed exactly].
   The dividing line among Clifford letters is locality: a Clifford family whose label action is banded is a QCA and
   transfers; a free letter g that carries a local diagonal Pauli Z_A to a Pauli string spread over all columns yields
   the free word w = g S_A g^-1 = sqrt of a global Pauli string, with |tr w_f|/N_f = 2^(-1/2) at every f, so every
   certificate satisfying (T) must send w and w^2 = (global Z string) to 1, and Prop C (ut-direct 4.1) then kills all
   column-uniform Paulis. Realising words: **S_par = chi S_last chi^-1** for chi = (row ^= row << 1) (ut-direct 4.1), and
   **k S_first k^-1 = i^(parity of row 0)** for the monomial letter k = (row 0 of column 0 ^= parity of row 0 over
   the other columns), the letter of v-ut-direct-logic's E2 counterexample [computed exact, f = 2..10,
   out_t4_parity_word.txt]. Non-Clifford letters (Toffoli, CS) together with chain maps are not covered; the ring
   mechanism of item 4 cannot reach Clifford + Toffoli: a necessary trace identity fails at every odd L from 3 to 41
   [computed, floats, s5.2].
7. **Correction to ut-direct 7.1/7.2** [read + proved]: Q as written there is not closed under products
   (v-ut-direct-logic/out_e2_counter.txt). Q_loc (s4.1) is closed (s4.2) and contains every letter Theorem E lists, so
   Theorem E and Theorem E+ should be read with Q_loc.

## 1. Setting and the question

One digit column = m qubits; Pauli labels v in V = F_2^(2m), symplectic form <,>; P_v the Hermitian Pauli with label v;
P_v P_w = i^(beta(v,w)) P_(v+w) with beta(v,w) in Z/4. A **Clifford QCA** alpha on the chain (sites j in Z) is an
automorphism of the Pauli group of the chain, commuting with the shift, with
    alpha(P_(v,j)) = eps(v) (x)_d P_(M_d v, j+d)        (finitely many d; eps(v) = +-1).
Its symbol is M(S) = sum_d M_d S^d. Scalars are ignored throughout (an automorphism of the Pauli group sees U only
through Ad U). Pi_0 (the *collapse*): Pi_0(alpha)(P_v) := eps(v) prod_(d increasing) P_(M_d v), a one-column operator.
ut-direct 7.5 asked whether Pi_0, after a sign correction, can be made a homomorphism into the one-column Clifford group
C_m (mod scalars), identity on the one-column letters, so that Corollary A's certificate sigma = Pi^(x)t applies.

## 2. The lag identity and its two consequences [proved]

**Lemma 1 (lag identities).** For every integer e and v, w in V: sum_d <M_d v, M_(d-e) w> = delta_(e,0) <v, w>.
*Proof.* alpha(P_(v,0)) and alpha(P_(w,e)) commute or anticommute exactly as P_(v,0) and P_(w,e), i.e. according to
delta_(e,0)<v,w>. Their symplectic product is the sum over sites s of <M_s v, M_(s-e) w>. QED

**Lemma 2 (gamma = 0; each collapse is a Clifford map).** Put A_d = P_(M_d v), B_d = P_(M_d w) (one-column operators).
(a) gamma(v,w) := sum_(d > d') <M_d v, M_(d') w> = sum_(e >= 1) sum_d <M_d v, M_(d-e) w> = 0 by Lemma 1.
(b) (prod_d A_d)(prod_d B_d) = (-1)^gamma prod_d (A_d B_d) = prod_d (A_d B_d): moving each B_d left past the A_(d')
with d' > d costs (-1)^(<M_d' v, M_d w>), and the total exponent is gamma(v,w).
(c) On the chain, alpha(P_v) alpha(P_w) = i^(beta(v,w)) alpha(P_(v+w)), i.e. eps(v) eps(w) (x)_d (A_d B_d) =
i^(beta(v,w)) eps(v+w) (x)_d P_(M_d(v+w)). Site by site A_d B_d = i^(beta(M_d v, M_d w)) P_(M_d(v+w)), so
eps(v) eps(w) i^(sum_d beta(M_d v, M_d w)) = i^(beta(v,w)) eps(v+w). Collapsing with (b):
Pi_0(alpha)(P_v) Pi_0(alpha)(P_w) = eps(v) eps(w) prod_d (A_d B_d) = eps(v)eps(w) i^(sum_d beta(M_d v, M_d w))
prod_d P_(M_d(v+w)) = i^(beta(v,w)) Pi_0(alpha)(P_(v+w)).
(d) Hence v -> Pi_0(alpha)(P_v) satisfies the Pauli multiplication table: with w = v it squares to 1 (so it is
Hermitian), and P_v, P_w commute iff their images do. So it extends to an automorphism of the m-qubit Pauli group
fixing i, i.e. an element Pi_0(alpha) of PC_m := C_m/mu_8, with label map M(1). QED

**Lemma 3 (Pi_0 is multiplicative).** For Clifford QCAs alpha, beta: Pi_0(alpha beta) = Pi_0(alpha) Pi_0(beta).
*Proof.* Fix u; put b_d = M^beta_d u and c_(d,d') = M^alpha_(d') b_d. On the chain,
(alpha beta)(P_(u,0)) = alpha(eps_beta(u) prod_d P_(b_d, d)) = eps_beta(u) prod_d eps_alpha(b_d) prod_(d incr)
[(x)_(d') P_(c_(d,d'), d+d')], an operator whose factor on site e is the ordered product over increasing d of
P_(c_(d, e-d)). The collapse map (x)_e X_e -> prod_(e incr) X_e is well defined on tensor products of one-site operators
(it is multilinear), so
    Pi_0(alpha beta)(P_u) = eps_beta(u) prod_d eps_alpha(b_d) * prod over (e, d) in lexicographic order of P_(c_(d,e-d)).
Since Pi_0(alpha) is an automorphism (Lemma 2), Pi_0(alpha)(Pi_0(beta)(P_u)) = eps_beta(u) prod_(d incr) Pi_0(alpha)(P_(b_d))
    = eps_beta(u) prod_d eps_alpha(b_d) * prod over (d, d') in lexicographic order of P_(c_(d,d')).
Same factors, two orders. Two factors (d1,d1'), (d2,d2') with d1 = d2 are in the same order in both. With d1 < d2 they
are in opposite orders iff e1 > e2, i.e. d1' - d2' > Delta := d2 - d1 >= 1. The sign between the two products is
(-1)^N with N = sum_(d1<d2) sum_(e' > Delta) [sum_(d2') <M^alpha_(d2'+e') b_(d1), M^alpha_(d2') b_(d2)>], and each
bracket is the lag-e' identity of Lemma 1 for alpha, = delta_(e',0) <b_d1, b_d2> = 0 because e' > Delta >= 1. So N = 0.
QED

**Lemma 4.** Pi_0(g^(x)infinity) = g for every one-column Clifford g (one site, nothing to reorder), and
Pi_0(Ad P_w^(x)infinity) = Ad P_w. QED

**Proposition (the class).** Let Gamma be a group of symbols, B the group of Clifford QCAs mod scalars over it.
(i) 0 -> V -> B -> Gamma -> 1 is exact with V = {Ad P_u^(x)infinity}; Gamma acts on V by pi_1: alpha Ad(P_u^(x)inf)
alpha^-1 = Ad(P_(M(1)u)^(x)inf) (the sign alpha Ad(P_u) alpha^-1 puts on P_v is (-1)^<u, M^-1(1)v> = (-1)^<M(1)u, v>).
(ii) b -> (image of b in Gamma, Pi_0(b)) is a homomorphism B -> pi_1^*(PC_m) that is the identity on V (Lemmas 3, 4).
A morphism of extensions with identity on the kernel exists iff the classes agree, so [B] = pi_1^*[PC_m] in
H^2(Gamma; V_(pi_1)), and the difference class is 0. Two such morphisms differ by z in Z^1(Gamma; V); requiring the
identity on constants forces z = 0 on Sp(2m,2). QED

The sign-defect question was therefore posed one level too low: the defect gamma is 0 element by element, and the
2-cocycle that could have obstructed multiplicativity is the coboundary of nothing: it is 0 on the nose.

## 3. The ring picture [computed; the explanation argued]

Wrap a rule of radius r on a ring of L > 2r + 1 sites, L odd. Logical generators X_q^(x)L, Z_q^(x)L anticommute in
pairs (L odd) and generate a copy of the m-qubit Pauli group (identify i^k X^a Z^b with i^k X^(a..a) Z^(b..b)).
Translation-invariant automorphisms preserve it (the product over j of the images of P_(u,j) is translation
invariant with label M(1)u on every site). Restriction gives Psi_L: B -> PC_m, a homomorphism wherever wrapping is.
Computed (t1_bulk.py, 317 elements with a valid L in 7..15, two seeds): Psi_L depends only on L mod 4;
Psi_L = Pi_0 for L = 1 mod 4 and Psi_L = c Pi_0 c (c: i -> -i, X, Z fixed) for L = 3 mod 4; both are multiplicative
(320 pairs). The predicted correction chi(u) = <M_odd u, M_even u> between Psi and Pi_0 is identically 0 (proof: with
A_e = sum_(d odd) <M_d u, M_(d+e) u>, Lemma 1 at lag e gives A_e = A_(-e), and <M_odd u, M_even u> = sum_(e odd)
A_e = sum_(e odd > 0) (A_e + A_(-e)) = 0).
Why it had to work [argued]: for L odd, e_1 = 1 + S + ... + S^(L-1) is a self-adjoint idempotent of F_2[S]/(S^L - 1),
so the label space of the ring is the orthogonal sum of the S = 1 sector (= V) and the rest; a symplectic (Clifford)
automorphism commuting with S respects the sum, and the Clifford/Pauli extension restricts along orthogonal summands.
For L even the S = 1 sector is not a summand ((S^L - 1) = (S^(L/2) - 1)^2 over F_2) and the logical X^(x)L, Z^(x)L
commute: no copy of the Pauli group.

## 4. The family level and Theorem E+

### 4.1 The class Q_loc
A family k = (k_f) of Clifford unitaries is in **Q_loc** if there are r >= 0 and a Clifford QCA kappa = kappa(k) of
radius <= r such that for all large f in F:
(i) for every Pauli P supported on columns [r, f - r): k_f P k_f^* = kappa(P);
(ii) every Pauli supported on B_r := [0, r) u [f - r, f) is mapped to a Pauli supported on B_(2r).
Members (by inspection of the conjugation rules; computed for every letter in t2_family.py): one-column Cliffords
(r = 0), directions, the kernel, boundary-local Cliffords (kappa = id), causal and anti-causal address maps
x -> G(S)x with G and G^-1 polynomial (X labels by T_f(G), Z labels by T_f(G^-1)^T, no signs), column-uniform CZ and Z4
couplings of bounded range, cyclic rotations of a row (r = 1; column f-1 goes to column 0, inside B_2).
Not members: chi = row ^= row << 1 (a deep Z goes to a string over all earlier columns; computed: 4, 6, 8, 10 columns at
f = 6, 10, 14, 18) and the parity letter k of v-ut-direct-logic.

### 4.2 Lemma 5: Q_loc is a group and k -> kappa(k) is a homomorphism [proved]
Products: with (r_k, kappa), (r_l, lambda) take R = r_k + r_l. (i) P on [R, f-R): l(P) = lambda(P) lies on
[r_k, f - r_k), so k(l(P)) = kappa(lambda(P)); kappa lambda has radius <= R. (ii) P on B_R: split P = P_1 P_2 with P_1 on
B_(r_l), P_2 on the rest; l(P_1) on B_(2 r_l), l(P_2) = lambda(P_2) on B_(R + r_l); then split again for k: the part
on B_(r_k) goes to B_(2 r_k), the rest to kappa(.) on B_(R + r_l + r_k) = B_(2R).
Inverses: let r' be the radius of kappa^-1 and R = 2r + r'. (i) For P on [R, f-R), Q = kappa^-1(P) lies on [r, f - r),
so k(Q) = P, i.e. k^-1(P) = kappa^-1(P). (ii) For P on B_R write k^-1(P) = Q_1 Q_2 (Q_1 on B_r, Q_2 on the rest);
k(Q_1) lies on B_(2r), so kappa(Q_2) = k(Q_2) lies on B_(max(R,2r)) = B_R and Q_2 on B_(R + r') within B_(2R).
kappa(k) is read off from k_f at any single large f, so it depends only on the family. QED

### 4.3 Theorem E+ [proved, given Theorem R-infinity and the cited certificates]
*Hypotheses.* Every letter is a family in Q_loc; every free letter is monomial (an address permutation with Z4 phases)
or boundary-local; paid letters are directions and rank-2 blocks (one-column Cliffords); the target is K = H^(x)m per
column. N in V^(x)t0 is a Theorem-R certificate kept by g^(x)t0 for every monomial m-qubit Clifford g.
*Claim.* sum over paid steps of f_N(p) >= W f_N(K). With the hook certificates: >= 2W (m = 2, V (x) Lambda^2 V),
>= 14W/5 (m = 3, f = 5, V (x) Lambda^4 V: cap 80, dem 224) and >= 41W/13 (m = 4, f = 9, V (x) Lambda^8 V: cap 27456,
dem 86592) [read: d-tof-cert/STATUS.txt 10:17:35 and out_verify_hook.txt; b4-depth4/STATUS.txt 19:01]. Those N are
span{e_x (x) e_T : x in T}, kept by every permutation and every diagonal phase g (e_x (x) e_T -> phase * e_(gx) (x)
e_(gT), and x in T gives gx in gT).
*Proof.* Take t >= t0 a multiple of 8, Y_0 = N (x) V^(x)(t - t0), and sigma(k) := (any lift of Pi_0(kappa(k)) to C_m)^(x)t
(C_m = the group generated by H, S, CNOT; its scalars are mu_8, so the t-th tensor power of a lift is well defined).
- sigma is a representation of the family group G: kappa is a homomorphism (Lemma 5), Pi_0 is (Lemma 3).
- Free letters keep Y_0: a monomial letter maps Z-type Paulis to +-Z-type Paulis, so every M_d of its bulk does, so
  Pi_0 maps Z-type Paulis to +-Z-type Paulis: Pi_0(kappa(k)) normalises the diagonal group, i.e. it is a monomial
  Clifford, and N is kept by its t0-th power. Boundary-local letters have kappa = id, sigma = 1.
- Paid letters and K are one-column: Pi_0 gives them back (Lemma 4); f_sigma = 2^(m(t - t0)) f_N on both sides.
- Trace condition (T): let sigma(k) != 1, so Pi_0(kappa(k)) is not a scalar.
  (a) Bulk symbol M_k(S) != I: columns of M(k_f) for deep generators are those of T_f(M_k) by (i), the other 4rm columns
  are arbitrary, so ut-direct E5 applies (rank(M(k_f) - I) -> infinity) and E4 gives |tr k_f|/N_f -> 0.
  (b) M_k(S) = I: then kappa(k) = Ad(P_u^(x)inf) and Pi_0(kappa(k)) = Ad P_u, so u != 0. By (i) every single-column
  Pauli P_(v,j) with r <= j < f - r is fixed by k_f up to the sign (-1)^<u,v>; so these labels lie in ker(M(k_f) - I)
  and the sign character of E4 is nontrivial on that kernel: tr k_f = 0 exactly.
Theorem R-infinity (ut-direct s2) now gives sum_p f_sigma(p) >= W f_sigma(K). QED

## 5. Where the obstruction lives

### 5.1 Non-banded Clifford letters [proved for chi and k; computed exact]
For a Clifford family with label matrix M_f, "every local Pauli goes to a local Pauli" for g implies it for g^-1
(M^-1 = J M^T J, and the transpose of a banded matrix is banded). So a free Clifford letter outside Q_loc maps some
local Pauli to a Pauli whose support grows with f. For a monomial letter with address matrix G (label matrix
diag(G, G^-T)) the spreading can always be put on a Z-type Pauli: if G^-T is not banded, g spreads some Z; if G is not
banded, g^-1 (Z action G^T) does. Suppose the spread Z-type Pauli Z_A can be taken on boundary columns (true for chi:
Z_last; for k: Z_first). Then w = g S_A g^-1 with S_A = (1+i)/2 + (1-i)/2 Z_A (a boundary Z4 phase, free) is free, has
trace ratio (1+i)/2 at every f, and w^2 = g Z_A g^-1. Every sigma with (T) has sigma(w) = sigma(g Z_A g^-1) = 1. When g Z_A g^-1 is a column-uniform
Z string (both known cases), ut-direct Prop C steps 3-5 kill every column-uniform Pauli and i^f, so no Pauli-faithful
certificate (hook family) survives. Realising words: chi S_last chi^-1 = i^(parity of row 1) (ut-direct 4.1, re-checked
by v-ut-direct-compute c1), and k S_first k^-1 = i^(parity of row 0) with square Z(row 0)^(x)f (t4, f = 2..10, exact).
Non-banded letters of other shapes (a long-range swap of two columns, non-translation-invariant maps) are neither in
Q_loc nor of this kind; not analysed.

### 5.2 Non-Clifford letters with chain maps [computed (floats) + argued]
The ring mechanism needs a translation-invariant subalgebra B = M_n (n = 2^m) of M_n^(x)L on which g^(x)L acts by
Ad g. Then V^(x)L = V (x) W with the shift acting on W only, so tr(g^L) = tr(g^(x)L tau) = tr(g) chi_W(g): tr g = 0
must force tr g^L = 0. For the Clifford group this holds for every odd L (E4: tr U = 0 iff the sign character is
nontrivial on ker(M - I), and eps_(U^L) = eps_U^L = eps_U on that kernel) [proved; sampled: 2,433 zero-trace words,
m = 2, 3]. For Clifford + Toffoli (m = 3) it fails at every odd L in 3..41 [computed, out_t3_ring_limit.txt]: the
ring construction cannot exist for letter sets containing Toffoli. For Clifford + CS and Clifford + T (m = 2), 3,875
zero-trace words found no witness (out_t3b_cs_search.txt), so there the limit rests only on an argument [argued, from
memory: Clifford + any non-Clifford gate is dense (Nebe-Rains-Sloane 2001 maximality of the Clifford group), the
identity extends to the Zariski closure, and (x^L + y^L + z^L)/(x + y + z) is not a polynomial]. This rules out one
mechanism only; whether a Pauli-faithful sigma exists for {chain maps, Toffoli or CS, directions} is open.

## 6. Calculations (all exact unless marked)

| script | output | what | result |
|---|---|---|---|
| cqca.py | (library) | Pauli/Clifford/QCA arithmetic with phases mod 4; collapse; ring restriction | |
| t1_bulk.py [seed] | out_t1_bulk.txt, out_t1_bulk_seed2.txt | 320 random QCAs (m = 1, 2, 3; 9/28/57 elementary layers incl. E-type, CZ couplings, row shifts, directions), 320 pairs | gamma = 0, Pi_0 valid and multiplicative, Psi_L = Pi_0 (L = 1 mod 4) / conjugate (L = 3 mod 4), chi = 0, ut-direct 7.5 relation: ALL PASS |
| t2_family.py [seed] | out_t2_family.txt, out_t2_family_seed11.txt | 104 random words of finite-f families (open chains, f up to 84 columns, m = 2, 3), incl. commutator words with Pauli bulk and boundary letters | deep action = composed bulk QCA (0 mismatches), corner condition, (T): 75 words with symbol != I have rank(M - I) strictly growing and |tr|/N <= 2^(-rank/2); 19 with symbol I and Pauli bulk have tr = 0; monomial letters collapse to monomial Cliffords; E^-1 Dcz E Dcz^-1 = CZ chain at f = 1..8; chi not a QCA: ALL PASS |
| t3_ring_limit.py | out_t3_ring_limit.txt | floats: tr g = 0 => tr g^L = 0 ? | Clifford m = 2, 3: holds; Clifford + Toffoli m = 3: fails at all odd L 3..41 |
| t3b_cs_search.py | out_t3b_cs_search.txt | floats, longer words | Clifford + CS, Clifford + T (m = 2): no witness in 3,875 words each |
| t4_parity_word.py | out_t4_parity_word.txt | k S_first k^-1 for the parity letter | = i^(parity row 0), square Z(row 0)^(x)f, |tr|/N = 2^(-1/2), f = 2..10 |

## 7. In the brief's format

(a) **Abstract statement.** Pauli-faithful certificates transfer along a letter class exactly when the class acts by
locality-preserving automorphisms with translation-invariant bulk: then evaluation at S = 1 lifts to the Pauli level
(the extension class of Clifford QCAs equals the pull-back of the one-column Clifford class; the splitting is the
ordered collapse = restriction to the S = 1 summand of an odd ring). Outside that class the obstruction is a free
square root of a global Pauli string.
(b) **Prediction.** gamma = 0; the naive collapse is multiplicative; the ring restriction depends only on L mod 4 and
equals the collapse for L = 1 mod 4; therefore 14W/5 at block 3 (and 41W/13 at block 4) hold for every unit built from
unimodular chain maps and the other Q_loc letters; and every non-banded monomial letter yields an S_par-type word.
(c) **Calculation.** t1-t4 above, at the smallest sizes (m = 1..3, f up to 84 columns): every prediction held.
For the open cases: 5-in-14 at block 3 stays exactly tight and cannot be beaten with unimodular chain maps; at block 2
the only Clifford door left is non-locality-preserving cheap letters (chi, parity-to-boundary), where every
Pauli-faithful certificate is dead (Prop C) and Prop D empties the Pauli-blind ones that treat E as one-column.
