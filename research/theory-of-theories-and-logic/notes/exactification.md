> **Working note, kept as written.** An AI agent wrote this during the theory-of-theories round of 14 October 2026 (the border track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts.
>
> **Refuted part.** The claim of item 3 and section 4, that both known repairs make a relation-using approximate unit
> cost at least a factor 2 per level, was refuted by a skeptic (`programs/exactification/skeptic_exactify_count.py`):
> the factor can be close to 1, so such units are worth whatever their exact version costs. Item 4's conclusion
> ("worthless at blocks 2 to 4") falls with it. Theorem X (formally polynomial units count at full value) stands.

# What does it cost to turn an approximate (border) unit into an exact algorithm?

meta-border, 14 October 2026, about 06:10-07:00 EDT. Written by one agent of the meta round (`../BRIEF.md`). No person
has read it, nothing is in Lean, nothing was posted. Labels: [proved] (full argument here), [computed] (program in
this folder, output named), [argued], [guess], [read] (read today at the place cited).

## 0. The decision, in five lines

1. **A border unit counts at full value if, and only if (for every exactification we know), its cancellations of
   negative powers of eps are formal**: they hold with the paid blocks treated as free symbols, using no relation
   of the paid letters (such as dir z · dir z = shift z, or two directions commuting). Call such a unit
   *formally polynomial*. [definition; d-border's NOTE.txt 3.5 had the idea]
2. **Formally polynomial units: the cost is a factor O(1 + a·L) = O(log log N)** (a = largest total negative
   eps-power along a path through one copy, L = recursion depth = Theta(log log N)). The saving exponent
   x = ln(mW/M)/ln m is unchanged: the work is N (log N)^(1-x) · O(log log N). [proved, section 3]
3. **Units that use a relation: one shared eps gives a wrong answer from depth 2 on** [proved by example; computed,
   `out_exact_toy.txt`]. Both known repairs (separate scales per level; exactify each copy by interpolation or by
   truncated arithmetic) multiply the work by at least 2 per level, so the unit counts as if it had 2M moves or more.
   [proved for those two methods; that no cleverer repair exists is OPEN]
4. **A relation-using unit is then worthless wherever the proved bound exceeds m/2 moves per array**: blocks 2, 3, 4
   (bounds 2, 2.8, 3.154 against 1, 1.5, 2) and block 5 in the engine game (2.665 against 2.5). It could still matter
   at block 5 in the widened game (31/14 = 2.21) and at block 6 (engine 2.641, widened 2.032, against 3). [proved
   arithmetic, given 3]
5. **So the forecast in the brief needs one more distinction.** Rank methods (Theorem R, RF, D) bound the border
   complexity in the group G (ladder B1). The exponent counts the *formal* border complexity. Between the two sit
   relation-using border units, which rank methods cannot exclude and which are worth nothing at blocks 2-4.
   (Section 6.)

## 1. What Bini and Schoenhage do [read: Blaeser, "Fast Matrix Multiplication", Theory of Computing Library,
Graduate Surveys 5 (2013), section 6, pp. 25-28; copy `lit/blaeser_gs005.pdf`. Originals not re-read today:
Bini, Capovani, Lotti, Romani, Inform. Process. Lett. 8 (1979) 234-235; Bini, Calcolo 17 (1980) 87-97;
Schoenhage, SIAM J. Comput. 10 (1981) 434-455.]

- Definition 6.1: R_h(t) = least r with r products of vectors over K[eps] equal to eps^h t + O(eps^(h+1)).
  Border rank = min over h.
- Theorem 6.3(3): R_(h+h')(t (x) t') <= R_h(t) R_h'(t'). **Under the tensor (Kronecker) product the approximation
  degrees ADD.**
- Lemma 6.4: R(t) <= c_h R_h(t) with c_h <= C(h+2, 2) (expand in eps and keep the eps^h part). Remark 6.5: over
  infinite fields c_h = 2h + 1 (interpolation in eps).
- Theorem 6.6 (Bini): if R_h(<k,m,n>) <= r then omega <= 3 log_kmn r. Proof: tensor up s times first, so the degree
  is 3hs, then make it exact once: R(<N^s,N^s,N^s>) <= c_(3hs) r^(3s). The factor c_(3hs) is polynomial in s while
  the gain is exponential in s, so it vanishes in the exponent.

Two facts make this work. (i) Recursion in matrix multiplication IS the tensor power, and along any path an entry
passes through one copy of the base algorithm per level, so degrees add. (ii) The cancellations in a bilinear
approximate algorithm are identities of multilinear algebra: they survive when the inner products are replaced by
approximate inner products (the inner algorithm computes eps^h' · product + O(eps^(h'+1)), and that is again a
bilinear map). In our language: every bilinear border algorithm is formally polynomial (after the usual
normalisation by eps^h).

## 2. Our recursion is substitution, not a tensor power

[read: d-border NOTE.txt sections 3.1-3.4, which I re-derived; the engine's recursion facts there are from
E/Recursion.lean, which I did not re-read.]
- **Depth.** A paid block of rank k at block size m on f = n/m digit columns is a sub-call on k·f bits. So the
  depth is L = Theta(log log N), and the saving is (log N)^x = n^x with n = log2 N. A factor c paid per level costs
  c^L = n^(log_m c): it lowers x by log_m c. A factor polynomial in L costs nothing in x.
- **Paths.** Inside one unit a path through the data crosses several paid moves one after the other (an array's
  m directions in the plain schedule). Each is a sub-call, a whole copy of the unit at the next level. So a path
  crosses up to p^j copies at depth j (p = most paid moves on a path in one unit, p >= m). The total eps-degree of
  the output is therefore up to about a·p^L >= a·log N, not h·L as in Bini's tensor powers. **Global interpolation
  in eps (Lemma 6.4 / Remark 6.5 applied once at the top) costs a factor of about log N and kills every saving.**
  [argued: an upper bound on the degree that is attained unless there is cancellation]
- **Substitution.** The sub-call returns B + eps E (B the exact paid block, E some error), not B. The outer unit's
  eps^-a terms cancel at X = B. Whether they still cancel at X = B + eps E is exactly the question of whether the
  cancellation is formal.

## 3. Formally polynomial units count at full value [proved]

**Setting.** A unit U at block size m serves W arrays with M paid moves. Write it as a linear straight-line program
whose free steps (mixing with coefficients in C[eps, 1/eps], free letters) are exact and whose paid steps call black
boxes X_1, X_2, ... (one symbol per paid block kind; sub-calls of conjugate directions may be written as one symbol
conjugated by exact free letters). Its formal output is Phi(eps; X), an element of (A ⊗ C[eps, 1/eps])^W where A is the
free product of the free letters' group algebra with the free algebra on the symbols. Let a be the largest total
negative eps-power along a path through one copy, and t the largest number of terms of a coefficient.
- U is a **border unit** if Phi(eps; B) = [K] e_l + O(eps) for each live output (B = the exact paid blocks).
- U is **formally polynomial** if Phi(eps; X) has no negative power of eps as an element of A ⊗ C[eps, 1/eps], and
  Phi(0; B) = [K] e_l.

**Theorem X.** If U is a formally polynomial border unit, then the recursive algorithm built from U, with one
global eps and every value stored as a truncated Laurent series, computes the transform exactly, and its operation
count is at most t·(1 + 2a(L+1)) times that of the same recursion run with exact scalars. Hence its work is
N (log N)^(1-x) · O(log log N) with x = ln(mW/M)/ln m.

*Proof.* (1) The composite map is a power series. Let Y(eps) = sum_(k>=0) eps^k Y_k be what a depth-(j+1) call
computes, with Y_0 = B (exact block). Then Phi(eps; Y(eps)) = sum_k eps^k Phi_k(Y(eps)), with Phi_k = 0 for k < 0
because U is formally polynomial, so it is a power series in eps whose constant term is Phi_0(B) = [K]. The f-fold
action is a representation of the algebra, so the same holds with every letter acting on f digit columns. By
induction on the depth (base case: exact small transforms), every call computes a power series in eps whose constant
term is exact. (2) Precision. In one copy, every coefficient applied along a path from a value to the output has
eps-valuation >= -a. So to know the output mod eps^q it suffices to know every sub-call output mod eps^(q+a), and a
sub-call's output mod eps^(q+a) depends only on its input mod eps^(q+a), because the sub-call map is a power series
(step 1). The top level needs mod eps^1. So depth j needs precision q_j = 1 + j·a. Inputs to a depth-j call have
lowest order >= -j·a (each level lowers it by at most a). So every value at depth j is a window of at most
1 + 2(j+1)a coefficients. (3) Cost. Adding two windows costs one operation per coefficient; multiplying by eps^h is
an index shift; multiplying by a coefficient with t terms costs t operations per coefficient; exact free letters and
base cases act coefficient by coefficient. So every operation at depth j costs at most t(1 + 2(j+1)a) exact
operations. (4) The top-level answer is the eps^0 coefficient: read it off, no interpolation. All of this uses only
additions and scalings by complex constants, which is what OpenAI's model allows [read: d-border NOTE 3.5 quotes
UniformFourier.lean; I did not re-read it]. ∎

So "border units count with a log log overhead" (d-border) is **right for formally polynomial units**. The
overhead is Theta(log log N) for this method (the deepest levels, where the work is, need the largest windows); it
is (log N)^o(1), so it costs nothing in the exponent. Whether it fits the slack in the record's exact statement is a
bookkeeping question (d-border read a (log log n)^(4-theta) slack in paperTime; not re-checked here).

**Relative version.** The proof only needs Phi_k(Y) = 0 for k < 0 at every Y that a sub-call can return. So the
cancellations may use any relation that the recursion keeps exactly: linearity; conjugation by exact free letters
(if the sub-call for dir z is the sub-call for dir e1 conjugated by an address map); commutation with a subgroup H of
free letters, if every sub-call is built H-equivariantly. They may not use relations of degree >= 2 in the paid
symbols (dir z · dir z = shift z; dir z · dir w = dir w · dir z), since an error E would have to satisfy
D E + E D = 0 or similar.

## 4. Units that use a relation

**Example (shared eps fails)** [proved; computed `exact_toy.py` -> `out_exact_toy.txt`]. Phi_r(eps; X) =
X + eps N + eps^-1 (X X - S) Z is a border unit for D (D D = S, a = 1), but Phi_r(eps; D + eps E) =
D + (D E + E D) Z + O(eps): the eps^0 coefficient is wrong at depth 2. In exact series arithmetic, the constant term
is off by 6.7, 2.6e2, 9.5e4, 2.0e10, 2.8e20 at depths 2-6, whatever the precision (1 + j or 1 + 3j). The formally
polynomial unit Phi_f(eps; X) = eps^-1 [X (I + eps Y + eps^2 N) - X] Y^-1 (same a = 1) has the right constant term at
every depth (errors 2e-16 to 6e-16).

**The two repairs and their price** [proved for these methods].
- *Separate scales.* Use eps_j at depth j. Depth j+1's error eps_(j+1) is multiplied by eps_j^(-a), so
  eps_(j+1) <= eps_j^(a+1) is needed: eps_j = eps^((a+1)^j), and the bottom precision is (a+1)^L =
  n^(log_m(a+1)). In floating point with a = 1, eps_j = 1e-2^(2^j) already needs 16 digits at depth 4: the toy's error
  is 0.03, 0.10, 0.16, 0.72 at depths 1-4 [computed].
- *Exactify every copy.* With exact sub-calls a copy's output is K + eps E_1 + ... + eps^h E_h; it needs h + 1 >= 2
  evaluations (interpolation, Remark 6.5's idea per copy), or a window of a + 1 >= 2 coefficients that each
  sub-call must then process separately (it is linear, and its own output must be exact). Either way the factor is
  c >= 2 per level, i.e. the unit counts as c·M moves.
- *The one exception:* with exact sub-calls a copy's output Phi(eps; B) - K is a polynomial in eps; if it is
  identically 0, every fixed value of eps already gives an exact unit with the same M moves, and nothing needs
  exactifying. Otherwise it is nonzero at all but finitely many eps, and the factor is c >= 2.
- *What is not proved:* that every exactification of a relation-using unit costs a constant factor per level. A
  proof would need a lower bound against all algorithms, which we do not have. [open]
- *Observed once* (NOTE.md 3.4): in the relaxed sub-game (a) sector the border-only schemes need W + 1 moves (W = 2, 3; W = 3
  probably) where the exact ones found need 2W, and they are relation-using with a = 1, so their known exactifications cost a factor
  a + 1 = 2: 2(W + 1) >= 2W, no gain left. Whether exactification always eats the border gain is a [guess].

**Consequence** [proved arithmetic]. With c >= 2, a relation-using unit saves only if c·M < m·W, i.e. M/W < m/2.
The certified lower bounds (b56 note; ladder Theorem D) are above m/2 at blocks 2, 3, 4 (2 > 1, 2.8 > 1.5,
3.154 > 2) and at block 5 in the engine game (2.665 > 2.5). So there, relation-using border units are worthless, and
only formally polynomial ones (or exact ones) matter. At block 5 widened (2.214) and block 6 (2.641 engine, 2.032
widened, against 3) a relation-using unit with c = 2 could still count.

## 5. A test every numerical hit must pass

[proposed by d-border NOTE 3.6; sharpened here.] Let x_k be points of a border family (loss -> 0, coefficients
growing). Replace every paid block B by B + delta·E with E random and delta FIXED (independent of eps) and evaluate
the word at x_k. Formally polynomial: the outputs converge to Phi_0(B + delta E) = K + O(delta) as k grows.
Relation-using: they blow up like eps^(-a)·delta. (With delta tied to eps, e.g. delta = eps^c, c > a, every border
unit passes, so the test is useless in that form.) Implemented in `bsearch.py`'s companion `perturb.py` for any hit;
no hit has needed it yet (section 7 of NOTE.md).

## 6. What this does to the forecast

- Rank methods (Theorem R, RF, D; any certificate) bound tau_border,G, the border complexity in the group G (ladder
  B1).
- The exponent counts tau_formal, the least moves of a formally polynomial border unit (or an exact one):
  tau_border,G <= tau_formal <= tau_exact. [proved: an exact unit is formally polynomial with a = 0]
- No rank method can see the difference between tau_border,G and tau_formal: a certificate is a representation of
  G, and a formal identity is only known to evaluate to K at X = B. [argued: a representation of the free product
  in which X acts as anything other than B gives no information on Phi(0; X)]
- So for each open target there are three cases:
  (a) no border unit in G: a rank proof may exist (not guaranteed: rank methods need not be complete);
  (b) a border unit in G but none formally polynomial: rank methods are hopeless AND the target is worthless at
      blocks 2-4 (section 4); a proof would need an exactness-sensitive method, but there is nothing to gain;
  (c) a formally polynomial border unit: rank methods are hopeless and construction is the path, at O(log log N)
      cost (Theorem X), with new Lean (d-border 3.7: a border certificate plus a recursion over truncated Laurent
      data; its 2.5k-5.5k-line estimate is a guess I did not redo).

## 7. Lean consequence (unchanged from d-border 3.7, re-read, not re-estimated)

The engine's GCert asks for an exact identity with rational mixing matrices; nothing in it can state "correct mod
eps". A formally polynomial unit needs (a) a certificate stating the identity in the free product over Q[eps, 1/eps]
(no negative powers) and its eps^0 evaluation, and (b) a recursion over truncated Laurent data with precision growing
by a per level, with the charging layer paying the window size.
