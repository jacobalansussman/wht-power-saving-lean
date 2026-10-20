# Fast Walsh–Hadamard transforms and stabilizer states: a dictionary and some questions

I work on exact fast Walsh–Hadamard and discrete Fourier transforms, and I count their arithmetic operations. Recent work gives exact transforms that use slightly fewer operations than the textbook ones. The bookkeeping behind these savings is the stabilizer formalism, under other names. I wrote this note for three reasons: to explain the correspondence to people who know that formalism, to credit what is already known, and to ask questions I cannot answer.

This note makes no claim about the overall complexity of these transforms.

I did this work in October 2026 with a team of AI agents (Claude) that I directed: I set the goals, chose which ideas to pursue and decided at each stage what to do next; the agents ran the computations, wrote the arguments and programs, and tried to refute one another's work. Nothing here is checked in Lean, and no person has read the arguments. The programs are in [../programs/quantum-dictionary/](../programs/quantum-dictionary/), and the full account, with the status of every statement, is in sections 6 and 7 of [../README.md](../README.md).

## The model

There are W registers. Each holds a vector of length 2^m, indexed by m bits. Think of it as data on m qubits, not as a physical quantum state.

- **The goal.** Each working register must end up multiplied by the kernel K = (√X)^{⊗m}, where √X = ((1+i)I + (1−i)X)/2. Scratch registers start at zero and are not checked at the end.
- **Free:**
  - on one register, any Hadamard-free Clifford: Paulis, CNOT circuits, CZ, S, powers of i;
  - across registers, any rational linear mixing: copying, adding, erasing.
- **Paid:** one √X on one qubit of one register, after free gates. I call this a move. Up to free gates and a scalar e^{iπ/4}, a move is one Hadamard.
- **Correctness.** The same schedule must stay correct when each gate U is replaced by U^{⊗f}, for every f. Because of this, a register is best read as a formal rational combination of terms g·x_j, where g is a Clifford element and x_j an input. Terms cancel only when they are equal.

The plain schedule costs W·m moves. A unit *saves* if it uses fewer. The published savings are built from such units, and they depend on the free mixing between registers. That mixing has no counterpart in a quantum circuit: a sum of two Clifford matrices is in general not a multiple of a unitary.

## The dictionary

Each row was checked by exact computation at small sizes by one AI agent, then re-checked by a second agent with its own code, written from the definitions in this project's Lean engine. These are checks at small sizes, not proofs for every m.

| Transform side | Stabilizer side | Checked |
|---|---|---|
| The m address bits of a register | m qubits | Notation |
| Address shift; sign pattern (−1)^{a·x}; label (x\|z) in F_2^{2m} | Pauli X; Pauli Z; Pauli X^x Z^z with the symplectic form | m ≤ 3 |
| Free steps on one register | The Hadamard-free Clifford group ⟨X, CNOT, CZ, S, i⟩. On labels it is the stabiliser of L_0 = {(0\|z)} | All elements at m = 1, 2; label images at m ≤ 3 |
| One paid move | √X on one qubit, after a CNOT circuit. Equals H up to free S gates and the scalar e^{iπ/4} | m ≤ 3 |
| The kernel | (√X)^{⊗m}. The Walsh–Hadamard matrix H^{⊗m} differs from it by free S layers, up to scalars | m ≤ 3 |
| A *frame*: the coordinate system of a register, up to free steps | A Lagrangian L(T) = T^{-1}L_0: the stabilizer group, signs dropped, of T^{-1}\|0…0⟩. Equivalently an orthonormal stabilizer basis | 3, 15, 135, 2,295 frames and 6, 60, 1,080, 36,720 states at m ≤ 4 |
| Frames built from a subspace U. The published certificates use only these | CSS states twisted by S on every qubit: Σ_{y∈U} i^{\|y\|}\|y⟩. There are 2, 5, 16, 67 at m = 1 to 4 | m ≤ 4 |
| The other frames | The other quadratic forms in the Dehaene–De Moor description | m ≤ 4 |
| Frames one move apart | Neighbours in the dual polar graph of Sp(2m,2) | m ≤ 4 |
| Frame distance m − dim(A∩B) | The least Hadamard count of the transition Clifford, which is the rank of its lower-left block. Also \|⟨ψ_A\|ψ_B⟩\|² ∈ {0, 2^{-d}} | All of Sp(2m,2), m ≤ 3 |
| Distance between subspace frames | The Koetter–Kschischang subspace distance of random network coding | m ≤ 4 |
| Which ports of the best known unit can serve each other (120 E8 root lines; orthogonal = can serve) | The commutation graph of the 120 four-qubit Paulis with an odd number of Y's. E7's 63 lines are the 63 three-qubit Paulis | All pairs |
| E8 octads (2,025); E7 heptads (135) | Four-qubit Lagrangians that are not totally singular; three-qubit Lagrangians | Exhaustive |

The correspondence stops being exact in four places:
- **Scalars.** Move counts ignore them, but exact correctness does not.
- **Mixing across registers.** It has no unitary counterpart.
- **The other half of the construction.** It uses labels over Q and has no stabilizer reading.
- **Fourier twiddle factors.** They are non-Clifford from order 4 on, so the dictionary covers the kernels only.

The frame distance is the Hadamard count of the transition operator. It is not the cost of steering one state into another: S takes |+⟩ to |+i⟩ with no Hadamard, yet their distance is 1.

## What is already known

Every single fact in the table is standard in your field:
- counts and overlaps (Aaronson and Gottesman);
- the Hadamard-free group and the Bruhat form (Bravyi and Maslov);
- Hadamard-optimal normal forms, with the count given by the block rank (Maslov and Roetteler; the same rank argument in Vandaele, Martiel, Perdrix and Vuillot);
- nearest-neighbour stabilizer states (Garcia, Markov and Cross);
- the quadratic-form description (Dehaene and De Moor; Pllaha, Tirkkonen and Calderbank);
- H^{⊗n}|A⟩ = |A^⊥⟩ (Aaronson and Christiano);
- three qubits and E7 (Cerchiai and van Geemen);
- the 135 heptads (Lévay, Planat and Saniga);
- the 120 E8 rays with their 2,025 orthogonal bases in Kochen–Specker work (Lisoněk, Raussendorf and Singh);
- the subspace metric (Koetter and Kschischang).

On the transform side, the Clifford reading of frames is not mine either.
- The idea of a frame comes from OpenAI's preprint on a power saving for the exact discrete Fourier transform.
- icekylinx, working with GPT-6 Astra, wrote the general Clifford frame lemma in the community repository. It gives the frame label T^{-1}L_0, the transition cost as a block rank, and the subspace frames with their distances.
- eumemic's pull request there already speaks of a Clifford model, stabilizer-state modules and the dual polar graph.

What I add is the explicit match with the quantum literature, checked by computer, and the two observations below. Neither observation is new mathematics; I give them because they are what the dictionary did for my problem.

## Two observations

**1. Shortest paths and CSS states.** Take two Lagrangians X and Y at distance d. The Lagrangians Z on shortest paths between them are exactly those with Z = (Z∩X) + (Z∩Y), and they form a lattice isomorphic to the subspaces of F_2^d.
- Status: argued by hand, not checked. Two independent AI skeptics re-derived the argument and found no error; no person has checked it. A computer checked it on every pair at m = 2 and 3 and on 94,095 pairs at m = 4.
- For X = the Z-type group and Y = the X-type group, it says: a stabilizer group lies on a shortest path from |0…0⟩ to |+…+⟩ exactly when it is CSS. Every other stabilizer state costs at least one extra Hadamard on any such path.
- In our setting this is why frames beyond the subspace frames do not help on their own.

**2. Cancellation with fixed offsets, and a rank bound.** Two pieces of content that ride in one register keep a fixed Clifford offset between them.
- With the right offsets, a target register can cancel content whose Pauli label anticommutes with its own. Without offsets, only commuting pairs can cancel. I have an explicit example on three qubits.
- The paid cost of the shared register is at least the F_2-rank of a symplectic cross-Gram matrix with ones on the diagonal.
- So for a pair, nothing is saved: two separate registers cost the same.
- Status: argued by hand, not checked. Both skeptics rebuilt the three-qubit example with their own code, and one extended the rank bound to registers that pass through several frames.

## Questions for you

1. **Hadamard count with free mixing.** In this model:
   - Hadamard-free Cliffords on each register are free;
   - rational linear mixing across registers is free;
   - each register needs (√X)^{⊗m}.
   
   Savings below W·m exist. Is there a quantum or simulation setting where gate counts are taken with free linear combinations, and are there monotones that survive such mixing? Sum-over-Cliffords methods and linear combinations of unitaries are the nearest things I know, but there combinations are expensive.
2. **Minrank of a Pauli commutation graph.** In the best known construction, ports are the 120 four-qubit Paulis with an odd number of Y's. A port can be served for free by another exactly when the two Paulis commute. The cost is driven by the minrank of this side-information pattern. A computer check finds that the largest set of pairwise anticommuting Paulis among the 120 has 8 members, so the minrank is at least 8 over every field. What is the minrank over F_2, over F_3 and over Q? Do Kochen–Specker results on these 120 rays, or Lovász-type bounds, say anything?
3. **A reference for observation 1.** Is it the equality case of the discrete uncertainty principle on Z_2^m? Is the general form stated in the literature on dual polar spaces?
4. **A reference for the E8 map.** The 120 root lines of E8 map to the 120 odd-Y four-qubit Paulis, with orthogonal lines going to commuting pairs. The octads map to the non-totally-singular Lagrangians. I see this as classical O+(8,2) geometry made explicit. Is there a standard source?

## Sources

- S. Aaronson, D. Gottesman, Improved simulation of stabilizer circuits, quant-ph/0406196 (Prop. 2; Sec. III).
- S. Bravyi, D. Maslov, Hadamard-free circuits expose the structure of the Clifford group, arXiv:2003.09412 (Sec. 2, Theorem 1).
- D. Maslov, M. Roetteler, Shorter stabilizer circuits via Bruhat decomposition and quantum circuit transformations, arXiv:1705.09176 (abstract; Theorem 13).
- V. Vandaele, S. Martiel, S. Perdrix, C. Vuillot, Optimal Hadamard gate count for Clifford+T synthesis of Pauli rotations sequences, arXiv:2302.07040 (Prop. 1).
- J. Dehaene, B. De Moor, The Clifford group, stabilizer states, and linear and quadratic operations over GF(2), quant-ph/0304125.
- H. J. Garcia, I. L. Markov, A. W. Cross, arXiv:1210.6646 (abstract) and On the geometry of stabilizer states, arXiv:1711.07848 (Theorem 13).
- T. Pllaha, O. Tirkkonen, R. Calderbank, Binary subspace chirps, arXiv:2102.12384 (Lemma 1).
- S. Aaronson, P. Christiano, Quantum money from hidden subspaces, arXiv:1203.4740 (Sec. 5.1).
- S. L. Cerchiai, B. van Geemen, arXiv:1003.4255.
- P. Lévay, M. Planat, M. Saniga, arXiv:1305.5689.
- P. Lisoněk, R. Raussendorf, V. Singh, Generalized parity proofs of the Kochen–Specker theorem, arXiv:1401.3035 (Example 2.7).
- R. Koetter, F. R. Kschischang, Coding for errors and erasures in random network coding, arXiv:cs/0703061 (eq. 3).
- R. F. Bailey, P. Spiga, Metric dimension of dual polar graphs, arXiv:1711.04012 (Lemma 1.2, after Brouwer–Cohen–Neumaier Sec. 9.4).
- J. Alman, K. Rao, arXiv:2211.06459.
- OpenAI, An explicit power saving for the exact discrete Fourier transform (preprint, September 2026).
- Community repository CrocSwap/integer-mult-bounds: research/native-complex-moment552/notes/general-clifford-frames.tex (icekylinx, PR #130); research/hyperplane-end-lean/HyperplaneEnd/Clifford.lean; PR #288 (eumemic).