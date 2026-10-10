# An argument that two arrays need six moves at block size 3 (Clifford frames, characteristic 0)

Jacob Sussman, 10 October 2026. It belongs to [for-network-coding-readers.md](../for-network-coding-readers.md), whose statement C(ii) it proves and whose section 1 it restates in the words of pull request #288.

In the words of that note, the claim of the title reads: in the operator model at m = 3, two messages need 6 paid steps and not 5. This file also uses three words of this repository. An *array* is a register that starts with one input (a message), a *move* is one unit of paid cost, and the *block size* is m.

I put this note together with AI agents working under my direction. The argument has had no human review yet and it is not formalised in a proof assistant. A computer has checked two things. The first is the arithmetic of the table in section 3, exactly: the checker included here and a second program that shares no code with it were written by the same AI agent within the same hour, and programs by five other AI agents give the same two averages, 41/35 and 3. The second is a test of Lemma 1 and Lemma 2 on small exact cases at block sizes 1 and 2, by one AI agent that took no part in writing them and that also read both lemmas line by line. A test of that kind is not a proof. I would be glad to be told where the argument is wrong.

## 1. The claim

Pull request #288 of CrocSwap/integer-mult-bounds (eumemic, 10 October 2026) sets up a "register model" for lower bounds. For the complex primitive at m = 3 it lists this case as open: "the full register model at W = 2, s = 5 with two or more inserted free elements" (its README, lines 106 and 351). It reports as its own results that s = 4 is impossible, that s = 5 is impossible with at most one inserted free element, and the bound s ≥ 50W/21 for Clifford frames. I have not checked those reports.

This note argues that for Clifford frames s = 5 is impossible whatever is inserted. So, for Clifford frames, coefficients in characteristic 0 and a paid element charged at least its length, W = 2 needs s = 6, which the plain schedule attains.

**Words.**

- *Pauli matrix, label.* For 3-bit vectors x and z, X^x Z^z is the 8 × 8 matrix that sends the basis vector e_y to (−1)^(z·y) e_(y+x). Its *label* is the 6-bit vector (x|z). The *form* of two labels is x·z′ + x′·z modulo 2; it is 0 exactly when the two matrices commute.
- *Frame.* A set of 8 labels that is closed under addition and on which the form is always 0 (a Lagrangian subspace). There are 135. L₀ is the frame {(0|z)} of the diagonal Pauli matrices. The *distance* of two frames is 3 minus the dimension of their intersection: 0, 1, 2 or 3.
- *Clifford matrix.* An invertible 8 × 8 complex matrix g such that g P g⁻¹ is a multiple of a Pauli matrix for every Pauli matrix P. It then acts on labels, additively and keeping the form, so it carries frames to frames and keeps distances. Write gL₀ for the image of L₀, and call the distance from L₀ to gL₀ the *length* |g|. Length 0 means that g has one non-zero entry in every row and column. A Hadamard matrix on one qubit has length 1, and so has a square root of X^z; H⊗H⊗H has length 3.
- *The machine.* G is any group of Clifford matrices, and K is a field of characteristic 0 (the coefficients). There are R registers, R ≥ W. A register holds a formal sum of terms c·[g]·y_l, with c in K, g in G and y_l the l-th input; two terms are of the same kind only when their matrices are equal. At the start register l holds [1]·y_l for l = 1, …, W, and the other registers (scratch) hold 0. A step is one of three things.
  - *Mix* (free): replace all registers at once by combinations of the registers with coefficients in K. Any R × R matrix is allowed, invertible or not. Copying, erasing and overwriting are mixes.
  - *Free map* (free): multiply one register on the left by an element b of G of length 0. Each term c·[g]·y becomes c·[bg]·y.
  - *Paid step*: multiply one register on the left by any g in G. Its cost is at least |g|.

  The run is *correct* if at the end every input l has a register of its own that holds exactly [T_l]·y_l, where |T_l| = 3. The number of *moves* is the total cost of the paid steps. The plain schedule uses three paid steps of length 1 per array: 3W moves.

**Theorem.** Every correct run uses at least 105·W/41 = 2.561·W moves. For W = 2 that is more than 5, so two arrays need 6 moves.

In detail: if n_d is the number of paid steps of length d, then n₁·41/35 + n₂·1307/630 + n₃·3 ≥ 3W.

**In the notation of #288.** A factor I + ([g] − 1)N with N = x·yᵀ and yᵀx = 1 equals C·diag([g], 1, …, 1)·C⁻¹ for an invertible matrix C over K whose first column is x: a mix, then g on one register, then a mix. So a product of such factors that equals [T]·(permutation) is a correct run with R = W, no scratch and invertible mixes. Cost-zero factors are mixes and free maps, in any number and at any place. If a paid factor is charged at least the length of its g (my reading of "rank" in README lines 39-41), the theorem gives s ≥ 105W/41 for Clifford frames. The machine allows more than that identity does: scratch registers, mixes that cannot be undone, and, through scratch, mixes whose coefficients are combinations of free maps.

## 2. Proof

**Lemma 1 (the score).** Let ρ be a representation of G on a finite-dimensional space V over K, and Y₀ a subspace of V with ρ(b)Y₀ = Y₀ for every b of length 0. Put f(g) = dim(Y₀ + ρ(g)Y₀) − dim Y₀. Then in every correct run the sum of f(g) over the paid steps is at least f(T_1) + … + f(T_W).

This is the "FO lemma" of #288, stated in the header of its script `complex_m2_lambda2_certificate.py`. I write a proof out because the machine above has scratch and arbitrary mixes.

*Proof.* Apply ρ to a state. The content Σ c·[g]·y_l of one register becomes the linear map (v_1, …, v_W) ↦ Σ c·ρ(g)v_l from V^W to V, and the R registers together become a linear map M from V^W to V^R. Let A = M(Y₀^W): what the registers hold when every input is taken from Y₀. The *score* is Φ = dim(A + Y₀^R) − dim Y₀^R, the number of dimensions by which A sticks out of "Y₀ in every register".

1. *Start.* A is Y₀^W in the first W registers, so Φ = 0.
2. *A mix or a free map never raises Φ.* It is a linear map of V^R that sends Y₀^R into itself. So it acts on V^R/Y₀^R, and the new A, taken modulo Y₀^R, is the image of the old one.
3. *A paid g on register r raises Φ by at most f(g).* Let D be ρ(g) on coordinate r and the identity on the others. Then dim(DA + Y₀^R) = dim(A + D⁻¹Y₀^R). The space D⁻¹Y₀^R lies in the space obtained from Y₀^R by enlarging Y₀ to Y₀ + ρ(g)⁻¹Y₀ in coordinate r. That space has f(g) more dimensions than Y₀^R, because dim(Y₀ + ρ(g)⁻¹Y₀) = dim(ρ(g)Y₀ + Y₀).
4. *End.* Keep only the W registers that deliver. There A consists of all (ρ(T_1)v_1, …, ρ(T_W)v_W) with every v_l in Y₀, and this sticks out of Y₀^W by exactly f(T_1) + … + f(T_W) dimensions. Dropping the other registers cannot raise a dimension, so Φ is at least that sum. ∎

**Lemma 2 (averaging).** Give every frame U a subspace Y_U of K^n. Any table will do; no symmetry is asked. Put δ(U, V) = dim Y_U − dim(Y_U ∩ Y_V), and let ā_d be the average of δ(U, V) over all ordered pairs of frames at distance d (so ā_0 = 0). Then there are ρ, V and Y₀ as in Lemma 1 with f(g) = N·ā_|g| for every g in G, where N is a positive number that does not depend on g.

*Proof.* Let Γ be the group of all maps of the 64 labels that are additive and keep the form, and N the number of its elements. Every g in G acts on labels by an element of Γ, which I also call g; products go to products. Γ permutes the frames and keeps distances, and it carries every ordered pair of frames to every other ordered pair at the same distance (checker, step 2). Let V be the space of all lists v = (v_γ), one vector of K^n for each γ in Γ, and let ρ(g) move the entries: (ρ(g)v)_γ = v_(g⁻¹γ). Let Y₀ be the set of lists with v_γ in Y_(γ⁻¹L₀) for every γ.

- ρ(g)Y₀ is the set of lists with v_γ in Y_(γ⁻¹gL₀). Indeed, if w = ρ(g)v then w_γ = v_(g⁻¹γ), which lies in Y_((g⁻¹γ)⁻¹L₀) = Y_(γ⁻¹gL₀). So ρ(b)Y₀ = Y₀ whenever bL₀ = L₀, that is, whenever |b| = 0.
- Y₀ and ρ(g)Y₀ are both given position by position, so dimensions add up over the positions: f(g) = Σ_γ δ(γ⁻¹gL₀, γ⁻¹L₀).
- As γ runs through Γ, the pair (γ⁻¹gL₀, γ⁻¹L₀) runs through the ordered pairs at distance |g| and visits each of them equally often. So the sum is N times their average. ∎

**The table.** The file `cert_E162_bitmodel.json` gives each of the 135 frames three integer rows of length 7, with entries between −2 and 2; Y_U is their span in Q^7. A computer search found it; how it was found plays no part in the proof. For this table (checker, steps 3 and 4):

| distance d | ordered pairs of frames | values of δ (number of pairs) | average ā_d |
|---|---|---|---|
| 1 | 1,890 | 1 (1,626), 2 (204), 3 (60) | 41/35 |
| 2 | 7,560 | 2 (6,996), 3 (564) | 1307/630 |
| 3 | 8,640 | 3 (all) | 3 |

A rank over Q does not change when Q is enlarged to K.

**End of the proof.** By the two lemmas, N·(sum of ā_|g| over the paid steps) ≥ N·W·ā₃ = 3NW. This is the detailed inequality. Since 1307/630 ≤ 2·41/35 and 3 ≤ 3·41/35, a paid step of length d adds at most d·41/35 to the left side, and it costs at least d. So (41/35)·(moves) ≥ 3W. For W = 2: five moves give at most 5·41/35 = 41/7 = 5.857, and 6 is needed. ∎

*Remark.* A wrong table cannot prove a wrong statement, because Lemma 2 holds for every table. A poor table only gives a poor bound.

## 3. Checking the table

The program uses the Python standard library only. The time is for one core of a laptop.

```sh
python3 -B check_two_arrays.py cert_E162_bitmodel.json          # 1 second
```

`check_two_arrays.py` (151 lines, integers and fractions only) is the check the proof rests on. It (1) lists the 135 frames from their definition and compares them with the file, (2) checks that the form-keeping maps of labels have exactly one orbit on ordered pairs of frames at each distance, (3) computes δ for all 18,090 ordered pairs by exact integer elimination, and the averages, (4) checks 5·ā₁ < 2·ā₃. Its last lines should read:

```text
4. THE TWO NUMBERS:  score per move  abar_1 = 41/35   ;   score per array served  abar_3 = 3
   every correct run for W arrays uses at least  W * 105/41 = 2.560976 * W  moves
   two arrays: 5 moves supply at most 5 * 41/35 = 41/7 = 5.8571 ; two arrays need 2 * 3 = 6
PASS: two arrays cannot be served with 5 moves at block size 3; they need 6
```

The table file has sha256 `1a749c87e27bd04d43224ca9323dfab8620a73a3b856d47c4c4fca6ddd7f453c`.

Three further programs are not included here and are available on request: one that reaches the same numbers from the 8 x 8 matrices of section 1 of the note (X^z, the 168 linear maps of addresses, the diagonal phases i^(z·x + c), and the paid matrix ((1+i) + (1−i)·X^z)/2), proving every rank twice over; one that runs the checker on twelve control tables (seven spoiled copies of the table fail; an earlier table with 5·227/189 = 6.005 fails; a weaker table with 5·1133/945 = 5.9947 passes; at block size 2 a table on the 15 frames gives 2 moves per array, the answer pull request #288 reports there); and one that tests how Lemma 2 is written down at block size 2, where Γ has 720 elements. One AI agent wrote these three and the checker. The AI agents that wrote the programs ran them; I list their results as reported.

A further AI agent, which took no part in writing the proof or the checker, tried to break both. Its own program recomputes the table by another route (frames grown from smaller subspaces, distances by a search in the graph, ranks with fractions) and gives 41/35, 1307/630 and 3 with the same counts of each value of δ. It spoiled the table in three ways of its own, in files that declare no averages: one frame replaced by 8 labels that are closed under addition but on which the form is not always 0; the 135 sets of rows shuffled over the frames; every third frame short of one row. The checker included here rejects all three. Two rough edges remain, and neither lets a wrong table pass: a table whose rows have unequal lengths stops the checker with a Python error instead of a line that starts `FAIL`, and the averages that a file declares are compared as text, so a file that writes 3/1 for 3 is refused.

## 4. Where I would look first for a mistake

1. **Whether the machine is #288's model.** I assumed: the frames are Clifford matrices; the free elements are exactly those of length 0; a paid element costs at least its length; the target has length 3. Lines 37-49 of #288's README are all I went by.
2. **Lemma 2.** It is the part this note adds, and it is half a page long. The place to check is the direction in which ρ(g) moves the list entries.
3. **Lemma 1.** The statement is #288's; the proof written out here, which allows scratch, is mine.
4. **The arithmetic.** The included checker and a second program that shares no code with it agree on the averages and on how often each value of δ occurs; one AI agent wrote both within the same hour. A program by another AI agent, by another route, gives the same averages and counts.

## 5. What is not covered

- **Frames or free maps outside the Clifford group.** #288 mentions an "enlarged model" (its README, line 81) and "Clifford plus 2^k-th-root phases" (line 100). Lemma 2 needs every matrix of G to permute the 135 frames, so it says nothing there.
- **A cost rule that charges a paid element less than its length.**
- **Coefficient fields of characteristic 3 or 5.** The same table gives the same two numbers modulo every prime from 11 up, and it still excludes five moves modulo 2 and modulo 7. Modulo 3 and modulo 5 it does not. (These figures modulo primes come from programs of the AI agents that are not included here; I have not rechecked them.)
- **Steps that are not linear, and identities that hold only approximately or only for some sizes.** The statement is about the formal identity. That correctness for all column counts is exactly this identity is #288's statement (README line 49), which I use as given.
- **W = 3 with s = 8.** The bound gives 315/41 = 7.68, so 8 stays open. A table with ā₃ = 3 and ā₁ < 9/8 would close it. I do not have one.
- **Other block sizes.**
- **Review.** No person has read the proof yet, and none of it is formalised.

## 6. Credit, and what this note adds

The register model, the FO lemma and the list of open cases are from #288 (eumemic; prepared with AI assistance by its own statement). #288 credits the lemma to its authors' earlier exploration and points to papers P1 and P2 of Swapnil Jain for overlapping statements; I have not read those papers.

This note writes out a proof of the lemma in a form that allows scratch registers and mixes that cannot be undone; I do not know whether the proof behind #288's statement already covers them. It adds the averaging step (Lemma 2), which lets a table without any symmetry serve as a certificate, and one table that reaches 105/41. I have not searched the literature for the averaging step. It puts the translates of a family under a group side by side, which is a standard way to build a representation, so it may well be written down.
