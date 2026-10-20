# famous-problems: small exact checks behind the famous-problems survey (14 October 2026)

These programs back the reshapings and reductions in sections 1 to 5 of [../../README.md](../../README.md). Each one
checks an identity or a reduction on small cases, exactly where it can (integers, fractions, Gaussian rationals,
finite fields) and in floating point otherwise. A check on small cases is not a proof of the general statement; the
README gives the argument and its status.

Run every command from the folder that holds the program, with Python 3.9 or later. Programs marked "numpy" need
numpy; the others use only the standard library. Every program here takes at most a few seconds and well under
200 MB. The outputs of the original runs are in [../../data/famous-problems/](../../data/famous-problems/), in a folder
of the same name. Every program was re-run for this revision on an Apple Mac with Python 3.9.6 and numpy 2.0.2, and
every output matched the saved one (apart from timing lines and numpy's harmless warnings about matrix products).

Programs whose name starts with `skeptic_` were written by an independent AI agent whose job was to break the
claim, with its own code.

## lower-bounds/

| program | what it checks | command | needs |
|---|---|---|---|
| `c1_wht_via_matmul.py` | the WHT on n² points is two products with a fixed n × n matrix (Alman 2021, Lemma 8.3 with k = 2), so a matrix-multiplication circuit gives a WHT circuit; and the linearisation step on small exact cases | `python3 c1_wht_via_matmul.py` | - |
| `c2_semigroup_ring.py` | why a unit is an identity in a semigroup ring and not a statement about one matrix: tensor powers of distinct products stay linearly independent once enough digit columns are used | `python3 c2_semigroup_ring.py` | numpy |
| `c3_numbers.py` | the arithmetic of Morgenstern's bound (how large some constant must be in a circuit of size C·N·(log2 N)^(1−s), for several sizes and savings; vacuous at every size anyone will run) and the other numbers of the lower-bound survey | `python3 c3_numbers.py` | - |
| `c4_dft_polymult_equivalence.py` | the DFT and polynomial multiplication in both directions (the convolution theorem and Bluestein's trick), on small lengths | `python3 c4_dft_polymult_equivalence.py` | numpy |
| `skeptic_v_checks.py` | the skeptic's own checks: the XOR form of Bluestein's trick for the WHT up to N = 128, Bluestein for n = 1 to 40, the superconcentrator remark and the potential argument | `python3 skeptic_v_checks.py` | numpy |

## kernels-and-matrix-multiplication/

| program | what it checks | command | needs |
|---|---|---|---|
| `fourstep.py` | the reduction from matrix multiplication to the WHT and to the DFT (the four-step method) on small cases, exact | `python3 fourstep.py` | - |
| `universal.py` | every invertible 2 × 2 kernel that is not a permutation times a diagonal gives the same exponent for its Kronecker powers: the explicit words H = U·diag(1,−2)·P·U·P and three calls to H giving U up to diagonals | `python3 universal.py` | - |
| `depth2.py` | Alman and Li's depth-2 rebalancing lemma can only ever give about N log N, for every fixed depth-2 exponent | `python3 depth2.py` | - |
| `skeptic_linearise.py` | Strassen's algorithm (n = 2 and 4) built as a straight-line circuit, one factor fixed to H, linearised gate by gate: the products come out exact within 6·C(n) gates; also a non-bilinear circuit with cancelling junk | `python3 skeptic_linearise.py` | - |
| `skeptic_xorbluestein.py` | (−1)^(x·y) = i^\|x\| · i^\|y\| · i^(−\|x XOR y\|), so the WHT is one XOR convolution between two diagonals, exact for n ≤ 6 | `python3 skeptic_xorbluestein.py` | - |

## network-coding-map/

| program | what it checks | command | needs |
|---|---|---|---|
| `checks.py` | the WHT, OpenAI's kernel C and XOR convolution; Bluestein and Rader; the shift identity H·diag(H e_s / N)·H = XOR shift by s that turns a fast circuit into a network carrying every shift; the edge-to-distance ratio of the plain FFT's own graph | `python3 checks.py` | - |
| `skeptic_vchecks.py` | the skeptic's checks: the shift identity for the cyclic DFT over the field with 97 elements, the WHT-to-C identity modulo 10009, the constant factor guaranteed for any degree-3 network (ratios 11.9, 9.6, 8.7, 8.4, 8.0 at log2 N = 8, 16, 32, 64, 1024), and the vertex counts of the frame graphs | `python3 skeptic_vchecks.py` | - |

## quantum-simulation/

| program | what it checks | command | needs |
|---|---|---|---|
| `stab_frames.py` | stabiliser states on m qubits (60 and 1,080 at m = 2, 3), their groups modulo sign (15 and 135 Lagrangians, the frames), and the frame distance as a Hadamard count | `python3 stab_frames.py` | numpy |
| `skeptic_vcheck.py` | one paid direction equals H·S·H exactly, and S·dir·S = e^(iπ/4)·H, so a paid move is one Hadamard up to free gates; the frame distance is additive along paths; the bound on nonzeros along paths that shows one array needs m paid moves | `python3 skeptic_vcheck.py` | numpy |

A third program of this survey (`f_uniform.py`, why "correct at every number of digit columns" makes Hadamards cost
something) is not included: it used about 8 GB and was stopped by the memory guard of the machine it ran on before it
finished.

## finite-geometry/

| program | what it checks | command | needs |
|---|---|---|---|
| `c1_e8_labels.py` | the community's E8 label family as a piece of the E8-mod-2 quadric: labels are the odd words of weight 3 or 7 in F_2^9, "knows" is even overlap | `python3 c1_e8_labels.py` | - |
| `c2_frustration_rm.py` | for relabellings that form a 2-group, the "frustrated" sign patterns of two-array separations are exactly the words outside the Reed-Muller code RM(1,k), as an exact linear-algebra statement | `python3 c2_frustration_rm.py` | - |
| `c3_frames_vs_quadric.py` | the 135 block-3 frames against the 135 singular points of the quadric O+(8,2): the two graphs compared | `python3 c3_frames_vs_quadric.py` | - |
| `skeptic_v1_iso.py` | the skeptic's explicit graph isomorphisms (E8 labels with root lines; block-3 frames with the quadric's points), and the two-piece separation by an order-3 letter | `python3 skeptic_v1_iso.py` | - |
| `skeptic_v2_fullpauli.py` | does the Reed-Muller separation criterion survive when the separating coefficients may be any combination of Pauli letters? | `python3 skeptic_v2_fullpauli.py` | - |

`out_c4_hirzebruch.txt` in the data folder is the output of a four-line evaluation of Hirzebruch's inequality whose
program was not kept; it is the arithmetic quoted in section 5 of the README (lhs 0 against rhs 16, 55, 352 at
multiplicity 1 for 4, 5 and 8 lines, so those configurations are not realisable over the complex numbers).

## quadratic-polynomial/

| program | what it checks | command | needs |
|---|---|---|---|
| `k1_baur_strassen.py` | the transform H and the single quadratic polynomial y^T·H·x have the same saving within a factor 15 (reverse-mode differentiation, Baur and Strassen, then linearisation), on H with up to 32 points and on random matrices | `python3 k1_baur_strassen.py` | - |
