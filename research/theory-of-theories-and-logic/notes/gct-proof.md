> **Working note, kept as written.** An AI agent wrote this during the second abstract wave of 14 October 2026 (the geometric complexity theory track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts. No skeptic has read it.

# a2-gct: proofs (14 October 2026)

Written by agent a2-gct (meta2 round). Nothing posted; no skeptic has read this. Labels: [proved] (proof here),
[read] (read today at the place cited), [cited from memory].

## Setting [read: ../../proofs/ladder/PROOF.md, Lemma 0, Rung 1, B1]

R = C[G], S = C[F] (F the free subgroup; finite in the engine game). Contents are rows in R^W. A scheme with M paid
moves (each a direction D, all conjugate to one D_0 by free letters, so we may take every paid move to be D_0 = h
after absorbing free letters into the S-coefficients) is a choice of coefficients

    r_i = sum_l beta_il e_l + sum_{j<i} gamma_ij u_j,   u_i = h r_i,   beta_il, gamma_ij in S,

and it serves W arrays iff K e_l lies in M_M = S e_1 + ... + S e_W + S u_1 + ... + S u_M for every l (Lemma 0).
Write Out_M for the set of W-tuples (P_1, ..., P_W) in (R^W)^W with every P_l in M_M for some choice of coefficients.
Out_M is the image of a polynomial map (coefficients beta, gamma and output coefficients c in S), so it is a
constructible subset of the affine space A = (R^W)^W, and its Zariski closure cl(Out_M) is an affine variety. The
target is T = (K e_1, ..., K e_W). A W-in-M unit exists iff T is in Out_M; a border unit (coefficients in C(eps),
right up to O(eps)) exists iff T is in cl(Out_M) [ladder B1 and the remark after it: a point of the closure of a
constructible set is reached along an algebraic curve].

Symmetry. H = GL_W(S) x GL_W(S) acts on A: (A, B).(P_l)_l = (sum_l' A_ll' P_l' B)... concretely, the left factor mixes
the W outputs with coefficients in S (free: outputs are S-combinations), the right factor applies free letters and
invertible mixing to the inputs. Out_M is H-stable, hence so is cl(Out_M) and its ideal I_M. S is semisimple (F finite),
so GL_W(S) = product over P_0-types tau of GL_(W d_tau)(C) is reductive. [proved: immediate]

## Proposition 1 (closed obstructions = border complexity) [proved]

Call a *closed obstruction* to a W-in-M unit any polynomial function f on A with f = 0 on Out_M and f(T) != 0. Every
occurrence or multiplicity obstruction of geometric complexity theory, every equation of an orbit closure, and every
rank (Theorem R) argument that holds for border units, is (or yields) a closed obstruction.

(a) A closed obstruction exists iff T is not in cl(Out_M), i.e. iff no border W-in-M unit exists.
(b) If one exists, one exists inside a single H-isotypic component of C[A]_d for some degree d, and a highest-weight
    vector of that component can be taken (the GCT form of an obstruction).

*Proof.* (a) f vanishes on Out_M iff it vanishes on cl(Out_M) (closure = zero set of I_M). So a closed obstruction is
an element of I_M not vanishing at T, which exists iff T is not in V(I_M) = cl(Out_M). (b) I_M is graded (Out_M is a
cone: scaling every output by a scalar is free) and H-stable; H is reductive, so each graded piece (I_M)_d splits into
H-isotypic components, f = sum of its components f_lambda, and f(T) != 0 forces f_lambda(T) != 0 for some lambda.
Inside that component, the H-translates of a highest-weight vector span it, and T is not a common zero of them (else
f_lambda(T) = 0), so some translate g.v with v a highest-weight vector satisfies v(g^-1 T) != 0; g^-1 T is in the
H-orbit of T, which is the GCT setting. QED

**Consequence.** The best possible closed obstruction (GCT included) proves exactly the border complexity
tau_border(W) = min{M : T in cl(Out_M)}. GCT-style obstructions beat Theorem R at W arrays iff tau_border(W) exceeds
the best Theorem-R bound at W. They never beat the border complexity: the border barrier B1 applies to them verbatim.

## Proposition 2 (relaxation duality) [proved]

Let Sigma be a finite set of irreducible representations of G, V_Sigma their sum, rho: R -> R_Sigma := image of R in
End(V_Sigma) (= the product of End(V_chi), chi in Sigma, by Burnside/Wedderburn, the chi being pairwise
non-isomorphic), S_Sigma = rho(S). The *relaxed game* is the same game with R, S, h, K replaced by R_Sigma, S_Sigma,
rho(h), rho(K); write Out_M^Sigma and T_Sigma for its output set and target.

(a) rho(Out_M) = Out_M^Sigma (rho applied entrywise).
(b) Every obstruction that factors through rho (a function f = f' o rho on A) and vanishes on Out_M vanishes at T
    whenever T_Sigma is in the closure of Out_M^Sigma; if f' is arbitrary (not polynomial), it vanishes at T whenever
    T_Sigma is in Out_M^Sigma.
(c) Every Theorem-R certificate (V, Y) whose representation V has all its irreducible constituents in Sigma, with any
    multiplicities (so every labelled certificate on copies), bounds the relaxed border complexity:
    tau_border^Sigma(W) >= W dem / f(D).

*Proof.* (a) rho is an algebra map with rho(S) = S_Sigma, so it maps every scheme of the game to a scheme of the
relaxed game with the same M and outputs rho(P_l): rho(Out_M) is contained in Out_M^Sigma. Conversely S -> S_Sigma is
surjective, so every relaxed coefficient lifts to S; the lifted scheme's outputs map onto the relaxed outputs.
(b) f' vanishes on rho(Out_M) = Out_M^Sigma, hence (if polynomial) on its closure, which contains T_Sigma, and
f(T) = f'(T_Sigma) = 0. The second claim is the same without closure. (c) The proof of Theorem R (ladder, Rung 1) uses
only the action of R and S on V, which factors through rho because the constituents of V are in Sigma (rho_V factors
through R_Sigma). Lemma 0 and Theorem R hold verbatim for the relaxed game (any algebra containing S_Sigma acting on
V), and B1 (run over C(eps)) gives the bound for border schemes. QED

**How it is used.** If a relaxed W-in-M scheme exists in R_Sigma (exactly or as a border scheme), no obstruction of any
kind that only sees the irreducibles in Sigma (rank, labels, GCT multiplicities, any polynomial; and, for exact relaxed
schemes, any function at all) can exclude a W-in-M unit. The relaxed complexity is a computable ceiling for every
certificate living on Sigma. Two remarks:
- If no P_0-type (irreducible of the free group's image) is shared between two parts Sigma_1, Sigma_2 of Sigma, then
  S_Sigma = S_1 x S_2 contains the central idempotents of the parts, the module splits, and
  tau^Sigma(W) = max(tau^Sigma_1(W), tau^Sigma_2(W)): coupling between irreducibles only happens through shared P_0-types.
  [proved: the idempotent (1, 0) in S_Sigma separates every content into its two parts; the paid moves are shared but
  each part can use them independently.]
- For one irreducible V with multiplicity-free P_0-types, S_V = sum over the types tau of End(tau). Writing U_tau for
  the row space of the tau-rows of the module, the target condition is U_tau contains the tau-rows of K e_l for all tau and l;
  this is the form run/relax2.py computes with.

## Lemma 3 (block 2, everything-free: finite-f relaxations are empty) [proved given the computation]

At a finite set Phi of digit-column counts, the relaxed game uses R_Phi = image of the letter algebra in
the product over f in Phi of End(C^(4^f)). If at each f in Phi the free letters act irreducibly (commutant = scalars),
then S_Phi = R_Phi (the C^(4^f) are pairwise non-isomorphic, having different dimensions, so by Burnside and the
Chinese remainder theorem for semisimple algebras the image of the free letters is the whole product), K lies in
S_Phi, and the relaxed game has a 0-move scheme for every W. Computation (xw/out_xw_commutant.txt): the commutant is
the scalars at f = 1, 2, 3 for the menus xw, xq and thmB. Hence no obstruction that factors through finitely many
column counts (rank, GCT, any) proves anything in these games, including the two that are proved closed (xq, thmB);
their lower bounds must use representations of the family group that factor through no finite f (as the QCA-symbol
representation sigma of meta/cohomology Theorem E+ does). [proved for f = 1, 2, 3; the statement for every finite f
is argued from the same letters (boundary Paulis plus chain or column-uniform letters generate every single-bit
Pauli), not computed beyond f = 3.]

## Literature (cited from memory; not re-read today)
- K. Mulmuley and M. Sohoni, "Geometric complexity theory I", SIAM J. Comput. 31 (2001); "GCT II", SIAM J. Comput. 38
  (2008).
- P. Buergisser, C. Ikenmeyer and G. Panova, "No occurrence obstructions in geometric complexity theory", J. Amer.
  Math. Soc. 32 (2019), arXiv:1604.06431.
- The point that every GCT obstruction is closed, hence bounds border complexity only, is standard (GCT works with orbit
  closures by design).
