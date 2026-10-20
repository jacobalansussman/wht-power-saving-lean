> **Working note, kept as written.** An AI agent wrote this during the logic round of 14-15 October 2026 (the ml-logic track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts. Its skeptic's verdict is [logic-classification-skeptic.md](logic-classification-skeptic.md): the core holds; the
> transfer criterion for limits is false; item 6 of Theorem A is argued in one direction only; the dated record is not genuine.

# ml-logic: proofs

meta3, 14 October 2026, written 21:15 to 22:05 EDT. No skeptic has read this. Nothing was posted.
Labels: [proved] = proved here (standard facts from logic are cited and marked [memory] or [read]); [computed] = a
program in this folder and its output; [argued] = reasoning short of proof standard; [guess].
Facts quoted from the literature carry how they were read: [read: first-hand], [read: abstract], [memory].

## 0. Setting

This follows `../../moonshot/theory-guide/RULES.md` §0 and `../../wave17/xc-legit/every_f.py`.

- A **game** Γ fixes a block m, a computable menu of free letters, the paid moves, and scalar mixing of registers.
  V = C^(2^m). At f digit columns a letter acts on V^(⊗f).
- A letter is **finite-state** if its action at every f is a matrix product operator (a weighted automaton over the
  alphabet of pairs (output digit, input digit)) of one fixed bond dimension b:
  `<o|L_f|i> = α^T A[o_1,i_1] ⋯ A[o_f,i_f] β`. Column-uniform letters have b = 1. Chain maps row_a ^= row_b << e have
  b = 2^e. Carries, neighbour phases (CZ chains) and column-0 boundary letters are finite-state
  (`every_f.py` docstring). A letter's entries may be unknowns.
- The paid move in direction z acts as D_z^(⊗f), D_z = aI + bX^z with a = (1+i)/2 and b = (1−i)/2 (b = 1).
- A **shape** S is everything except the numbers: block, W, M, the registers, the sequence of steps, which letter each
  step uses and which register it touches, which entries are unknown. The unknowns are x ∈ C^n.
- The **unit condition**: for every delivering register l and input j, the transfer T_lj(x)_f from input j to
  register l equals δ_lj K^(⊗f) for every f ≥ 1, where K = D_{e_1} ⋯ D_{e_m}. A mixed target replaces K by a product
  of some of the D_{e_k}. Scratch registers are unconstrained.

## 1. "Correct for every f" is a finite condition in every game we use

### Theorem F (one bound per shape) [proved]

For every shape S whose letters are finite-state there is a number N(S), computable from S alone, such that for every
x ∈ C^n: the unit condition holds for all f ≥ 1 iff it holds for f = 1, …, N(S).

*Proof.* Track the bond dimension n_rj of each transfer T_rj as an automaton. Initially T_jj = I^(⊗f) (n = 1) and
T_rj = 0. A mixing step r ← Σ_s μ_s r_s takes the direct sum of the automata (dimensions add, boundary vector scaled by
μ_s). A letter of bond b on register r takes the product automaton (n_rj multiplied by b). A paid move has b = 1. So
every n_rj is a function of S only, and the entries of the automaton matrices are polynomials in x.

The difference Δ = T_lj − δ_lj K^(⊗f) is an automaton of dimension n = n_lj + 1. Its value at f columns is zero iff
the series of Δ vanishes on all words of length f. Let R_k = span{α^T A_w : 1 ≤ |w| ≤ k} ⊆ C^n. Since
R_{k+1} = R_1 + R_k·A, the chain stops growing once two consecutive terms agree. If R_1 = 0 the series is zero on all
nonempty words. Otherwise dim R_k ≥ k until the chain stops, so R_n = R_∞. The series vanishes on every nonempty word
iff β ⊥ R_∞ = R_n, that is iff it vanishes on all words of length 1..n. (This is Schützenberger's equivalence bound
for weighted automata; Tzeng 1992 gives the polynomial-time version [memory]; `every_f.py` uses it with exact
arithmetic.) For each fixed x this is linear algebra over C, and n does not depend on x. Take
N(S) = 1 + max_lj n_lj(S). ∎

### Corollary F1 (engine game: one bound per block, the formal identity) [proved; numbers computed]

Suppose every letter is column-uniform and all letters, paid moves included, lie in a finite group G ⊂ GL(V). In the
engine game they do: Paulis, affine maps and Clifford Z4 phases are monomial Clifford elements, and D_z is Clifford
with entries in Z[i, 1/2]. Let r'(m) be the number of distinct **nonzero** values of tr(g) on G. Then for every
content c ∈ C[G] these are equivalent:
(i) Σ_g c_g g^(⊗f) = 0 for every f ≥ 1; (ii) the same for f = 1, …, r'(m); (iii) c = 0 in C[G].
So a column-uniform word is right at every f iff its formal identity holds in C[G] iff it is right at f = 1..r'(m).
The bound depends on the block, not on the word.

*Proof.* (iii) ⇒ (i) ⇒ (ii) are clear. For (ii) ⇒ (iii), the Burnside–Brauer argument: let χ = tr, let a_1 = dim V,
and let a_2, …, a_r' be the other nonzero values. Let ψ be an irreducible character of G, and suppose ψ occurs in no
V^(⊗f) with 1 ≤ f ≤ r'. Then for those f, 0 = |G|⟨χ^f, ψ⟩ = Σ_i a_i^f S_i, where S_i = Σ_{g: χ(g)=a_i} conj ψ(g). The
terms with χ(g) = 0 drop out. The matrix (a_i^f), f = 1..r', is a Vandermonde matrix times diag(a_i), so it is
invertible and every S_i = 0. But χ(g) = dim V forces g = I, because g has finite order and so all its eigenvalues
are 1. So S_1 = ψ(1) ≠ 0, a contradiction. Hence every irreducible representation of G occurs in some V^(⊗f) with
1 ≤ f ≤ r'. If c acts as 0 there, it acts as 0 on every irreducible representation, so c = 0 by Wedderburn. ∎

**Computed** (`fbound.py`, `out_fbound.txt`; letters D_j, X_j, Z_j, S_j, CZ, CNOT and the scalar i):
- m = 1: |G| = 96, 13 distinct traces, r' = 12.
- m = 2: |G| = 46080, 21 distinct traces, r' = 20.
- m = 3: a sample of 200,000 random-walk elements shows 29 distinct traces, so r' ≥ 28.
- In all three, |tr|² takes the values 0, 1, 2, 4, …, 4^m.
- Each nonzero modulus carries exactly 4 trace values. That fits the traces lying in Z[i, 1/2], where the only
  elements of norm 2^k are units times (1+i)^k.
- If |tr|² ∈ {0} ∪ {2^k : 0 ≤ k ≤ 2m} for Clifford elements [memory], then r'(m) ≤ 4(2m+1). That is 12 and 20 at
  m = 1 and 2, as computed.

**Sharpness at m = 1** [computed]: the evaluation map C[G] → ⊕_{f ≤ F} End(V^(⊗f)) has rank 4, 14, 34, 57, 77, 91,
95, 96 for F = 1..8. Its Gram matrix is Σ_{f ≤ F} tr(g^{-1}h)^f. So F* = 8 ≤ r' = 12. The script finds a content that
vanishes at f = 1..7 (norms ≤ 1.1e-14) and has norm 37.9 at f = 8. A finite test with f ≤ 7 can therefore be fooled,
at least on general contents.

### Corollary F1' (widened game: one bound per word) [proved]

If the letters are column-uniform but generate an infinite group, take the support instead. Let g_1, …, g_k be the
distinct group elements in the support of c − δK. Then the unit condition at f = 1..k already forces the formal
identity. Choose v, w with the numbers λ_t = w^T g_t v distinct and nonzero. They exist because g_t − g_s ≠ 0 and the
g_t are invertible. Then Σ_t c_t λ_t^f = 0 for f = 1..k, and Vandermonde gives c = 0. This is the in-house "by
Vandermonde" argument with its bound made explicit.

The widened game needs it. Adding the widened letter CS = diag(1,1,1,i) to the engine letters at block 2 makes the BFS
pass 300,000 elements by depth 8, against 46,080 without it (`widened_check.py`, `out_widened_check.txt`).
Clifford + CS is universal [memory], so that group is infinite and no block-uniform bound of the F1 kind exists there.

### Corollary F2 (per-shape decidability) [proved]

For every shape S in E, XQ, Xchi, Wd and XW, in the mixed-target and spectator variants, and in OpenAI's finite-win
setting, "a unit of shape S exists" is one existential sentence over C with Gaussian-rational constants: ∃x (finitely
many polynomial equations = 0, and finitely many ≠ 0 where the shape demands them). The equations are the unit
condition at f ≤ N(S), or the formal identity.

Such a sentence is decidable. ACF_0 is complete (Tarski–Chevalley), and concretely the system has a complex solution
iff 1 is not in the ideal (Nullstellensatz; Rabinowitsch's trick handles ≠), which Buchberger's algorithm decides.
OpenAI's finite win is the case of Definition 2.1 in their preprint [read: first-hand]: A nonmonomial (a ≠ condition),
b, and a word with g < bq^(b−1).

The same holds for **border units** of ε-degree ≤ d. The unknowns become the coefficients of Laurent polynomials in ε.
Theorem F holds over the field C(ε), because the bond dimensions do not change.

### Where the finite reduction would fail [computed / memory]

- **Letters that are not finite-state.** A Cooley–Tukey twiddle diag(exp(2πi j·val(w))) on the digit columns has
  Hankel rank L + 1 on words of length ≤ L (L = 2..6). A column-uniform phase has rank 1 and a CZ chain has rank 2
  (`hankel.py`, `out_hankel.txt`). So twiddles are not weighted automata (Fliess / Carlyle–Paz rank criterion
  [memory]). A game with twiddles as free letters would need another every-f argument. The prime-factor (CRT) layout,
  OpenAI's route, avoids them. None of our WHT games has such letters.
- **A demand for rational numbers.** Each shape would then be a Diophantine question over Q, whose general
  decidability is open (Hilbert's tenth problem over Q) [memory]. Any complex numbers are acceptable
  (`../../meta/ex-note/NOTE.md` §1), so this does not arise.
- **No uniform bound in the chain-map games.** E_e = (row0 ^= row1 << e) agrees with the identity for f ≤ e and
  differs at f = e+1. Its automaton has dimension 2^e, so the bound of Theorem F is 2^e + 2 ≥ e + 1
  (`every_f.py` Part B: first difference at f = 2, 3, 4, 5 for e = 1..4). The needed f grows with the shape. In the
  engine game it does not (F1).

## 2. The arithmetic classification

Let Sh_Γ be the decidable set of shapes of game Γ. Let U_Γ ⊆ Sh_Γ be those for which a unit exists, and let s(S) be
a shape's saving: ln(mW/M)/ln m, or for a mixed target the τ-transfer value. Let s*_Γ be the exponent saving of
programs built from units. By Theorem C (spectrum, unchecked) together with ex-laser Theorem 3 (unchecked),
s*_Γ = sup_{S ∈ U_Γ} s(S).

### Theorem A [proved, given Theorem C for the identification of s* with the supremum]

1. **U_Γ is decidable** (Corollary F2), and the comparison s(S) > q is decidable for rational q = a/b.
   For a full unit it reads (mW)^b > m^a M^b; for q = 0.1, (mW)^10 > m M^10. For a mixed unit (saving 1 − τ with
   Σ_k W_k k^τ = M) it reads Σ_k W_k k^(1−q) > M, a comparison of real algebraic numbers.
2. **"Some unit has saving > q" is Σ_1**: ∃S (S ∈ U_Γ ∧ s(S) > q). Its witness is a shape plus a consistency
   certificate: a reduced Gröbner basis without 1, with the cofactors that prove it generates the same ideal, or an
   exact point in a number field. Either can be checked in time polynomial in its size, so the sentence is ∃w φ(w)
   with φ bounded (Δ_0).
3. **"No unit has saving > q", that is s*_Γ ≤ q, is Π_1.**
4. **s*_Γ is left-c.e.** (computable from below): run through the shapes, decide each, and keep the running maximum
   of rational lower approximations of s(S). This sequence rises to s*_Γ.
5. **s*_Γ is computable iff it is also right-c.e.**: iff the set {q ∈ Q : s*_Γ < q} is c.e.
6. **The full circuit model L0 is a level higher.** "Complex linear circuits for H_(2^n) of size ≤ c·2^n·n^(1−s)
   exist for every n" is Π_1, since each n is a finite decision over C as in F2. So "s*_L0 > q", which reads
   ∃ rational s > q ∃c ∀n (…), is **Σ_2**, and "s*_L0 ≤ q" is Π_2. With a "for infinitely many n" definition it
   moves up again. The unit framework is what lowers existence to Σ_1. One finite circuit for H_(2^n), used in every
   direction, saves only the constant factor n^(−s): cost k·2^(n(k−1))·S for H_(2^(nk)). A growing saving needs the
   recursive amortisation of a unit. So L0 collapses to Σ_1 exactly if every circuit saving is realised by a finite
   unit, which is not known.

**Calibration.** OpenAI's "a finite win exists" (their Theorem, `02-generation.tex` Def. 2.1) is a Σ_1 sentence of
this kind. They proved it by contradiction, through a price function obtained by "separation and compactness"
[read: first-hand, introduction §1]. Theorem A then predicts that explicit finite wins exist and can be found by
search. They do: OpenAI's companion "explicit power saving" manuscript, the hub's units, and the 5/6 result. The plain
schedule is the trivial member of U_Γ at every block.

## 3. Where independence can and cannot live

### Theorem I [proved, from Theorem A and standard facts]

1. **The existence side has no independence.** If s* > q then "s* > q" is a true Σ_1 sentence. Robinson's Q proves it
   (Q is Σ_1-complete [memory]), so every theory that contains Q proves it. A proof is: exhibit the shape and its
   certificate, then verify them. Its length is about the size of the witness plus the verification run.
2. **Bounded impossibility has none either.** "No unit with at most B steps and R registers has saving > q" is
   decidable. When it is true, EA proves it with a proof about as long as the decision procedure's run.
3. **Unbounded impossibility: independence happens in every sound theory iff s* is not computable.** For a sound c.e.
   theory T, let q_T = inf{q : T ⊢ "s* ≤ q"}. Then q_T ≥ s*.
   - If q_T = s* for some sound c.e. T, then s* is right-c.e., because one can enumerate T's proofs. With item 4 of
     Theorem A, s* is then computable.
   - Conversely, suppose s* is computable, by a program producing rationals a_k with |a_k − s*| ≤ 2^(−k). Take
     T = EA + the single true Π_1 sentence "for all k, t and S: if the program outputs a_k within t steps, then no
     unit of shape S has saving > a_k + 2^(−k)". Then q_T = s*.
   - So: **if s* is not computable, every sound c.e. theory (ZFC, ZFC with large cardinals, Lean's) leaves a gap**.
     There is a q > s* for which "no unit has saving > q" is true and unprovable.
   - If s* is computable, some sound theory proves every true ceiling, though not necessarily ZFC.
   - Neither case is known. The same reasoning applies to each block separately (s*_Γ,m).
4. **The 0.1 bar.** "s* ≤ 0.1" is one Π_1 sentence. If it is true, it may or may not be provable in ZFC. We know no
   mechanism that would make it independent. Every ceiling we have (Theorem R and its relatives: rank potentials,
   induction on word length, finite group computations) is finitary. It is the sort of argument that formalises in
   IΣ_1, the first-order part of RCA_0 [read: first-hand, Wikipedia "Reverse mathematics"] [argued]. The Lean proofs
   use a stronger theory, but their content is finite linear algebra.

**Comparison with the literature.**
- **Matrix multiplication ω** is the mirror image. Constructions approximate it from above: rank is decidable over C,
  so ω = inf_n log_n R(⟨n,n,n⟩) is right-c.e. [memory]. Whether asymptotic tensor rank is computable "is unknown".
  Christandl–Hoeberechts–Nieuwboer–Vrana–Zuiddam (arXiv:2411.15789) prove it "computable from above" in a strong
  sense: the sublevel sets are Zariski-closed and decidable. They also prove discreteness from above: upper bounds on
  ω from bilinear maps "snap" to ω [read: abstract].
- **The Shannon capacity of graphs** has our shape exactly: a supremum over finite objects (independent sets of
  powers), so left-c.e., with Lovász/Haemers-type certificates from above. "The Borel–Turing computability of the
  Shannon capacity of graphs is completely open" (Boche–Deppe, arXiv:2010.06873) [read: abstract].
- **Chaitin's Ω** shows that a left-c.e. real need not be computable [memory].
- **MIP* = RE** (Ji–Natarajan–Vidick–Wright–Yuen 2020) makes the entangled value of nonlocal games, a supremum over
  finite-dimensional strategies just as s* is a supremum over finite units, uncomputable as a function of the game
  [memory].

Our s* is a single number, though, and uncomputability results speak about families. A single attained supremum is
computable: if some unit attains s*, then s* = ln(mW/M)/ln m.

## 4. Provability logic: no self-referential or reflection shortcut to existence

Let σ be Σ_1, for instance "some unit has saving > 0.1". Let S and T be c.e. theories containing Q, with T ⊇ EA, and
write Pr_S for S's provability predicate.

### Theorem P [proved from standard facts]

1. **Formalized Σ_1-completeness**: EA ⊢ σ → Pr_S(σ) [memory; Hájek–Pudlák].
2. **Collapse.** If S and T are Σ_1-sound, the following are equivalent:
   - σ is true;
   - Q ⊢ σ;
   - S ⊢ σ;
   - Pr_S(σ) is true;
   - T ⊢ Pr_S(σ);
   - T ⊢ σ.

   (Item 1 and Q's Σ_1-completeness give "σ true ⇒ everything". Σ_1-soundness of T, then of S, gives
   "T ⊢ Pr_S(σ) ⇒ S ⊢ σ ⇒ σ".) "Proving that a proof of existence exists" is therefore exactly as strong as proving
   existence, and it is true exactly when existence is true.
3. **Löb.** T ⊢ (Pr_T(σ) → σ) iff T ⊢ σ ["PA proves Prov(⌜B⌝) → B only when PA already proves B"; read: first-hand,
   SEP "Provability logic"]. An argument inside T of the form "if units are provable then they exist" is already a
   proof of existence.
4. **Modal reasoning alone proves nothing.** Every Σ_1 sentence p satisfies p → □p. Any argument that uses only this,
   the derivability conditions and Löb (GL) is valid for every Σ_1 sentence, ⊥ included. ⊥ is Σ_1, satisfies
   ⊥ → □⊥, and is not a theorem. So no such argument concludes p. An existence proof has to use the unit equations
   themselves. (By Solovay's theorem GL captures all that PA can say about its own provability in modal terms
   [read: first-hand, SEP].)
5. **Reflection and stronger axioms add no Σ_1 truths.** T + RFN_Σ1(T), T + Con(T), ZFC + large cardinals: each that
   is Σ_1-sound proves exactly the true Σ_1 sentences, the same set as Q. A stronger theory can change only the
   **length** of the proof of σ, never whether σ holds or is provable.
6. **This is the one exception, and it is about length.** Strong theories can shorten proofs of true bounded
   sentences enormously (Gödel's speed-up, Parikh; e.g. T + Con(T) over T) [memory]. So a short abstract proof that a
   unit exists, whose smallest unit is too large to write down, is logically possible, and it would count for the
   challenge. Logic gives no way to find such a proof. It would be ordinary mathematics about the unit equations, a
   "lemma" in the sense of §5.
7. **Nonstandard and limit witnesses.**
   - (a) A model of T in which σ holds with a nonstandard witness shows only Con(T + σ), that is T ⊬ ¬σ. It is
     evidence that impossibility is unprovable, not that a unit exists.
   - (b) **Transfer criterion for a fixed shape.** If a standard shape S has an exact solution in some algebraically
     closed field of characteristic 0 (an ultrapower of C, Puiseux series), it has one in C (ACF_0 is complete). An
     approximate family whose residual tends to 0 transfers only if the residual is exactly 0 in the ultrapower. A
     border family with infinitesimal nonzero residual does not transfer, and the plain schedule trivially does. This
     is the precise form of what ex2-limits must check.
   - (c) **Primes (Lefschetz principle)** [memory]: a fixed shape S has a complex unit iff it has a solution over the
     algebraic closure of F_p for infinitely many primes p, iff it has one for all sufficiently large p.
     - Effectively, if no complex unit exists then N = Σ g_i P_i over Z for some nonzero integer N, and then no
       solution exists over F̄_p for any p not dividing N.
     - So a solution over F_p, or over F_(p^k), for a single prime p above an explicit arithmetic-Nullstellensatz
       bound on N proves a complex unit exists [memory: Krick–Pardo–Sombra 2001 give such bounds].
     - This is the exact logical content of ex2-modp's "solutions mod infinitely many p" route. The bound is huge in
       general.

## 5. Proof mining: what an existence proof says about the unit's size, and what searches say about proofs

Write σ = ∃w φ(w), with φ bounded and w = (shape, certificate).

### Proposition M [the growth classes are memory, standard proof theory; the converse is derived here]

Let B_T(ℓ) be the largest "smallest witness" over all T-proofs of length ≤ ℓ of Σ_1 sentences. This is computable
when T is Σ_1-sound: run the search. B_T eventually dominates every function that T proves total. The reason is that
"∃y y = F(n)" has a T-proof of length O(log n) when T proves F total. Orders of growth [memory]:

| Existence proof of length ℓ in | smallest witness bounded by about |
|---|---|
| a universal (open) theory, cut-free | ℓ: the witness is an explicit term in the Herbrand disjunction |
| first-order logic with cuts (lemmas); EA = IΔ_0 + exp | a tower of exponentials of height growing with ℓ (level F_3). Non-elementary cut elimination is unavoidable in general (Statman, Orevkov) |
| IΣ_1, PRA, WKL_0 | Ackermann-like (level F_ω). WKL_0 is Π^0_2-conservative over PRA and covers Heine–Borel compactness and separable Hahn–Banach [read: first-hand] |
| PA, ACA_0 | level F_(ε_0) |
| ZFC | computable, but eventually above every function ZFC proves total |

**Converse** [derived here]. Suppose an exhaustive search shows that no unit of shape size ≤ B exists, in some family.
Then no existence proof for that family of the forms above is shorter than the inverse of the bound at B:
- **cut-free (explicit) proofs**: shorter than about B. This is meaningful: an explicit construction would have been
  found;
- **proofs with lemmas**: shorter than an iterated logarithm of B. For any practical B (10^6 skeletons) that is a
  handful of symbols, so the search says nothing about abstract existence proofs.

**The algebraic part is bounded, the combinatorial part is not** [memory: Bézout / Heintz; arithmetic Nullstellensatz,
Krick–Pardo–Sombra]. For a fixed shape with n unknowns and equations of degree ≤ d:
- if a complex unit exists, one exists with algebraic coordinates of degree ≤ d^n, and with explicitly bounded
  heights;
- the numbers are therefore never the obstacle once the shape is known;
- only the shape's size escapes every logical bound. That is the quantity searches explore and abstract proofs hide.

**Calibration.** OpenAI's finite-win proof uses separation and compactness. If it formalises in WKL_0 [guess: their
packing step needs only countably many matrix types at a time], proof mining predicts a finite win of
primitive-recursively bounded size. One exists, explicitly, in their companion manuscript.

## 6. Calibration summary

| Case | What this track says | Agrees with |
|---|---|---|
| plain schedule W in mW | in U_Γ at every block; formal identity holds; consistent | trivially |
| 8-move control word (`../../moonshot/o3/o3-build/`) | column-uniform, so F1 applies: exact rational solution + formal identity ⇒ right at every f | o3-build's own check (formal identity; raw arrays f = 1..5) |
| block-2 stand-in 2 in 3 (ex-dimension) | a fixed finite-dimensional system, decidable; their exact solutions are witnesses | `../../meta/ex-v-ex-dimension-logic/solve23.py` |
| E_e chain shift (XW-type letter) | needs f = e+1 > any fixed test; Theorem F bound 2^e + 2 catches it | `every_f.py` Part B |
| m = 1 content passing f ≤ 7 | refused by F1 (r' = 12 ≥ F* = 8) | `out_fbound.txt` |
| engine 2 in 5 at block 3 | a Π_1 statement; each size bound is decidable (refuted by finite search); the unbounded statement needs a uniform proof (Theorem R), finitary, formalisable in IΣ_1 [argued] | ladder / Lean |
| the 755 Gröbner-{1} skeletons (w3-asym) | each is a bounded instance; GB = {1} is an EA-checkable refutation (a Nullstellensatz certificate) | `../../moonshot/w3-asym/out_singular*.txt` |
| OpenAI's finite win | Σ_1, proved by contradiction ⇒ an explicit witness must exist and search finds it | companion explicit manuscript; hub and our units |

## Sources

- In-house: `../BRIEF.md`; `../../meta/MAP.md` §0–2; `../../meta/ex-note/NOTE.md`; `../../exist2/BRIEF.md`;
  `../../iso/NOTE.md`; `../../meta/literature/NOTE.md` T2; `../../moonshot/BRIEF2.md`;
  `../../moonshot/theory-guide/RULES.md` §0; `../../wave17/xc-legit/every_f.py` and `out_every_f.txt`;
  `../../moonshot/w3-asym/NOTE.md`; `../../moonshot/o3/o3-build/verify_ctl.py`.
- OpenAI, "Finite tensor savings and exact Fourier circuits", 25 Sep 2026, `build/sections/01-introduction.tex`
  §1.1–1.2 and `02-generation.tex` Def. 2.1 (local copy under `<private path>`)
  [read: first-hand].
- Christandl, Hoeberechts, Nieuwboer, Vrana, Zuiddam, "Asymptotic tensor rank is characterized by polynomials",
  arXiv:2411.15789 [read: abstract].
- Boche, Deppe, "Computability of the zero-error capacity of noisy channels", arXiv:2010.06873 [read: abstract].
- Stanford Encyclopedia of Philosophy, "Provability logic" (Löb, Solovay) [read: first-hand].
- Wikipedia, "Reverse mathematics" (WKL_0 conservativity, RCA_0's first-order part IΣ_1, Heine–Borel and Hahn–Banach
  at WKL_0) [read: first-hand].
- From memory, not re-read: Schützenberger 1961, Tzeng 1992 (automaton equivalence); Burnside–Brauer (tensor powers
  of a faithful character); Fliess / Carlyle–Paz (Hankel rank); Hájek–Pudlák (formalized Σ_1-completeness); Gödel
  1936 and Parikh 1973 (speed-up); Herbrand; Statman 1979 and Orevkov 1979 (non-elementary cut elimination);
  Kohlenbach, Applied Proof Theory (2008); Krick–Pardo–Sombra 2001 (arithmetic Nullstellensatz); Heintz 1983;
  Ji–Natarajan–Vidick–Wright–Yuen 2020 (MIP* = RE); Coppersmith–Winograd 1982 (ω is not attained); Lefschetz
  principle (e.g. Marker, Model Theory, §2.2); Glaudell–Ross–Taylor (Clifford + CS universal).
