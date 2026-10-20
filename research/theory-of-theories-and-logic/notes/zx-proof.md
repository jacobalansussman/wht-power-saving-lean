> **Working note, kept as written.** An AI agent wrote this during the second abstract wave of 14 October 2026 (the ZX track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts.
>
> **Corrected part.** The proof assumed that every free step only permutes and rephases entries, which is false in
> the games with Hadamards on end columns. Both skeptics supplied the repair: p*f + 2mWL >= mWf at every f forces
> p >= mW. Its "block-5 near miss of 0.008" (2 in 6) is not a near miss: the two-layer bound W + m excludes it.

# a2-zx: what ZX/ZXW rewriting can and cannot do for the paid-move count (14 Oct 2026)

Labels: [proved] = argument written here in full; [computed] = program + output in this folder; [argued]; [read].
No skeptic has read this. Nothing posted.

## 0. Dictionary [computed checks in zx2_slot.py, zx3_ancilla.py]

- One digit column holds m bits; a unit must work for every number f of columns.
- Paid move D_z = (1+i)/2 I + (1-i)/2 X^z in every column = the X-spider with phase pi/2 on the f wires of row z
  (sqrt X), a non-monomial Clifford. Up to free S gates it is a Hadamard (S sqrtX S = e^{i pi/4} H).
- Free letters of every game in the lattice (Paulis, affine maps, chain maps, Z4 phases, CZ, CS, Toffoli) are
  MONOMIAL matrices (one nonzero entry per row and column).
- Two encodings of W arrays:
  - tensor encoding: W registers of mf qubits; no linear combinations of arrays can be written (sum-free);
  - direct-sum (slot) encoding: one slot register (dimension W) next to one data register; a sum of arrays is a
    linear map on the slot; a paid move on one array (a combination of arrays) is the gadget
    G(Q, z) = I + Q (x) (D_z - I), Q a rank-1 idempotent on the slot. This is the ZXW reading (sums inside a diagram).
- Correct cost of a gadget exp(i pi/4 Q (x) (1 - X^z)) whose slot weight Q is diagonalisable with integer
  eigenvalues: the number of ODD eigenvalues (D^2 = X^z is a free Pauli). An uncontrolled sqrt-X on a data row has
  Q = I_W and costs W. PyZX counts it as one gadget (zx2 B: "one array" for W = 2).

## 1. Theorem 1 (sum-free routes are closed in every game) [proved]

**Statement.** Let C be any circuit on n = mWf qubits made of monomial gates (any arity) and single-qubit gates, and
suppose C equals K^{(x)W} = sqrtX^{(x) mWf} up to a scalar. Then C has at least mWf single-qubit gates that are not
monomial. Hence every sum-free scheme, and every circuit PyZX extracts from one (HAD, XPhase, ZPhase, CNOT, CZ, NOT,
CCZ, ...), has at least mW paid moves (each paid move is f such gates), in the engine, chi, widened and everything-free
games alike.

**Proof.** For a matrix A let s(A) be the largest number of nonzero entries in a column. A monomial matrix P
satisfies s(PA) = s(A) and s(AP) = s(A). A single-qubit gate g on wire q has at most two nonzero entries per
column of g (x) I, so s((g (x) I) A) <= 2 s(A). Starting from s(I) = 1 and multiplying gate by gate, s(C) <= 2^t
where t is the number of non-monomial single-qubit gates. Every entry of sqrtX^{(x) n} is nonzero, so
s = 2^n = 2^{mWf}, and t >= mWf. []

Remarks. The bound uses only the linear map, not the diagram, so no rewriting of a sum-free diagram (spider fusion,
local complementation, pivoting, gadget fusion, phase teleportation, extraction) can go below it. Calibration
(out_zx1.txt: 324 rows; out_zx4.txt: 54 rows with CS, CCZ and chain maps): 0 rows below mWf, as forced.

## 2. Lemma 2 (gadget fusion with slot weights never saves) [proved]

**Statement.** Let Q_1, Q_2 be rank-1 idempotents on the slot and A = Q_1 (x) (1 - P), B = Q_2 (x) (1 - P) with the
same data Pauli P. The fusion identity exp(i pi/4 A) exp(i pi/4 B) = exp(i pi/4 (A + B)) used by gadget_simp,
supplementarity and phase teleportation holds for all angles only if [Q_1, Q_2] = 0. Commuting rank-1 idempotents
are either equal or mutually annihilating. In the first case the fused weight 2Q_1 gives D^2 = P on one array (a
free Pauli: a redundant pair is removed, which a plain schedule does not contain). In the second case the fused
weight Q_1 + Q_2 is a rank-2 idempotent of cost 2, the same as before. So fusion never lowers the correct cost of a
scheme without redundant pairs.

**Proof.** Comparing the theta^2 terms of exp(theta A) exp(theta B) and exp(theta (A + B)) gives
theta^2 [A, B] / 2, so the identity for all angles forces [A, B] = 0. Since [A, B] = [Q_1, Q_2] (x) (1 - P)^2 and
(1 - P)^2 = 2(1 - P) is nonzero, [Q_1, Q_2] = 0. PyZX's gadget rules fuse only identical Pauli strings, which commute,
so they act within this case. Write Q_1 = a_1 b_1^T and Q_2 = a_2 b_2^T with b_k^T a_k = 1. Then
Q_1 Q_2 = (b_1^T a_2) a_1 b_2^T and Q_2 Q_1 = (b_2^T a_1) a_2 b_1^T. If they are equal and nonzero, a_1 is parallel to
a_2 and b_1 to b_2, and the normalisation gives Q_1 = Q_2. Otherwise both products vanish, Q_1 + Q_2 is idempotent of
rank 2, both eigenvalues are odd (1, 1), and the cost is 2. []

## 3. Proposition 3 (P4(a): commuting slot weights split the scheme) [proved, slot-uniform model]

Model: free steps are slot maps T (x) I and data words C applied to every array alike; paid moves are G(Q_j, P_j)
with the slot weights Q_j written in the input frame (conjugated by the slot maps before them) and the data Paulis
P_j written likewise (conjugated by the data words before them).

**Statement.** If the Q_j commute pairwise, the scheme costs at least mW.

**Proof.** By Lemma 2's argument, pairwise commuting rank-1 idempotents are equal or mutually annihilating, so they
lie among W' <= W mutually annihilating rank-1 idempotents E_1, ..., E_W' that sum to at most I. Complete these to a
basis of slot idempotents. In that basis every G(Q_j, P_j) is block diagonal, with D_{P_j} on block k(j) and I
elsewhere. The product of the gadgets is sum_k E_k (x) U_k, where U_k is the product of the D_{P_j} for j in J_k
(U_k = I for a block with no move). Correctness reads sum_k E_k (x) U_k = T^{-1} (x) C^{-1} K. Multiplying on the left
by E_k gives E_k (x) U_k = E_k T^{-1} (x) C^{-1} K with both slot factors nonzero, so U_k is a scalar times C^{-1} K
for every k. Each D_{P_j} is F^{-1} D_z F with F monomial, i.e. f non-monomial one-qubit gates. Theorem 1 with W = 1
then gives |J_k| >= m for each of the W blocks (a block with no move would make C^{-1} K a scalar, which is false).
The total is >= mW. []

This model leaves out per-array data letters and the free column-parity idempotents e_psi (sums of shifts). Those are
exactly what g-play's commutative slice exploits (section 6). The general case is not claimed here.

## 4. Remark 4: free states and effects would make paid moves free [computed]

zx3_ancilla.py: |+>^f on an ancilla row, CZ in every column, and <+|^f on the data row together give exactly H^{(x)f}
(f = 1, 2, 3). Every piece is column-uniform and phase-free. So a ZX model of our game must forbid states and effects
on data wires, i.e. post-selection. Inside the game the effect needs prod_c (I + X_c)/2, a sum of 2^f words. A single
sum of free images gives only (I + X^{(x) f})/2.

## 5. Remark 5: any search at one fixed f is degenerate [computed]

out_zx2.txt (C): at block 2 with f = 1 the 3072 engine letters span all 16 dimensions. So sqrtX (x) sqrtX is a free
combination, and two arrays cost 0 paid moves. At f = 2 the span of g (x) g has rank 88/256 and K (x) K is not in it
(relative residual 0.612). Solving f = 2 and f = 3 jointly leaves 0.822. All content of the problem sits in the
"same word at every f" requirement. Every lower bound we own uses it, and every ZX/ZXW search must impose it (several
f jointly, or SZX thick wires, Carette-Horsman-Perdrix arXiv:1905.00041).

## 6. Remark 6: ZX's strongest fragment is the commutative slice, already analysed [argued]

- X-phase gadgets of angle pi/2 on parities v are exactly dir(v) = exp(i pi/4 (1 - X^v)) in every column.
- Products of commuting gadgets with slot weights and sums are g-play's commutative slice: shifts, directions and
  additions (wave16/g-play/comm.py).
- The spider-nest identities that drive ZX T-count reduction (de Beaudrap-Bian-Wang arXiv:2004.05164), taken at
  angle pi/2, are the Fano closures prod_{v in U \ 0} dir(v) = 1 for dim U >= 3 (pf-logging F1, computed there).
- Inside the slice, sharing can only come from the free idempotents e_psi. These are combinations of shifts, which in
  the slot encoding are slot gadgets controlled by column parities. With them, one paid move serves different
  psi-sectors with independent vectors.
- The two-point test is the fixed-space argument applied to pairs of characters in one psi-coset; I re-derived this
  today. A product of k pseudo-reflections fixes a subspace of codimension at most k, so whenever rho(kappa) = -1,
  at least W moves must have rho(z_j) = -1.
- Status of the slice, read in proofs/logging/PROOF.md section 5 and moonshot/mw2-b5beyond/NOTE.md:
  - no saving at block 3;
  - at least 15W/4 moves at block 4 (closed for the 0.1 bar, which needs 3.48W or fewer);
  - 2 in 7, 2 in 8, 2 in 9 and 3 in 11 excluded at block 5 by exhaustive test.
- So ZX phase-gadget rewriting, the fragment where ZX is strongest, lands on ground that is already searched.

## 7. What ZX does not see

- Our cost is the rank of the slot weight of each D-gadget.
- ZX rewrites preserve the linear map. Their cost measures (Hadamard edges, T-count, the Hadamard count after
  extraction) over-count free CZ ladders (out_zx2 A: 8, 13 and 26 H-edges at paid count 0) and under-count
  uncontrolled gadgets (out_zx2 B: m f instead of 2m f/f per array).
- Completeness for Clifford+T (Jeandel-Perdrix-Vilmart arXiv:1705.11151) and for ZXW (Poor et al. arXiv:2302.12135)
  says every equation between diagrams is derivable. It says nothing about which diagram minimises our cost.
- With sums, our game is a ZXW-diagram optimisation whose cost is not local. A search inside it is the template
  search already run (c2-xw, moonshot).
