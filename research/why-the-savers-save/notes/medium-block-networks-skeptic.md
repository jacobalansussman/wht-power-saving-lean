> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# v-cx-medium: skeptic's check of cx-medium (a crossing unit at blocks 13 to 25)

Anatomy round 3, 15 October 2026, 13:22 to about 13:55 EDT (Mac clock). I am the skeptic of
`../cx-medium/` and tried to refute it.

**Labels:**
- **[computed]**: a program of mine in `calc/` and its output.
- **[proved]**: a complete argument written here.
- **[read]**: taken from a named file.

All programs are my own. Where the thing under test *is* cx-medium's code (its network builder, its shrink layout), I
import it and check it with my own semantics. I reran our own checker. Two short jobs ran on Colab machines;
nothing of mine is still running. Nothing was posted.

## Verdict: refuted in part

**What stands.** I rebuilt every number and every proof that I could reach, with my own code. All of them reproduce:
- the prices;
- the clique floors;
- E7's rigidity;
- the exact raw-array calls;
- the checker rerun.

So the headline stands. No configuration that cx-medium built in blocks 13 to 25 beats the hub, and nothing built comes
near 1e-2.

**What fails.** Two of the claims the note rests its conclusions on are wrong.

1. **The cause it gives is not the cause.** cx-medium says the one design that saves is held to ×0.44 because one pair
   has no idle levels for hosting helpers, "the hub's only way to make junk cleaning free". In cx-medium's own
   calibrated model this is false:
   - Hosting does not make junk cleaning free. In the hub's own layout (four pairs, block 73, exact fill
     1.6647327e-3), the hosted helper pieces carry **66.1 %** of the narrow mass. On banks at one pair they carry
     68.8 %.
   - Hosting is worth only **×1.195** there: putting every helper on ideal banks instead gives 1.3929e-3.
   - At one pair, perfect hosting would give **exactly** the ideal-bank price, 7.3057e-4, because the histogram is the
     same.
   - Most of the ×0.44 comes from the centre's share of each call. Each call keeps 10.7 free directions at one pair and
     35.6 at four pairs, while the helper mass per call is the same.
2. **The list of template points in range is incomplete.** The note says the shrunk template has five points in 13 to
   25 and that only one of them saves. But sub-designs with an even window give template points at **block 16**
   (window 6) and **block 22** (window 8).
   - Block 16 closes: E7's 32-line section has a clique floor exactly at break-even.
   - Block 22 does not close. The 56 lines of E8 at 60° to one root have a clique floor of **28 per call**, against a
     one-pair break-even of **37.3**. That gap is larger than E7's (35 against 42).
   - No decoder below 60 per call is known for it, so this is not a candidate. It is an open point that the note's
     conclusion ("medium blocks closed for the record"; "the only openings left are an E7 or a 112-port decoder")
     leaves out.

## 1. What I rebuilt, and what matched

| item | cx-medium | mine | how mine differs |
|---|---|---|---|
| clique floor A4 / A5 / D5 / D6 (per call) | 8 / 12 / 16 / 24 | **8 / 12 / 16 / 24** | own root coordinates (undoubled); **no labels built**; ranks taken from rows of B = J + Gram mod 2; subspaces as label-spanned closures (`calc/vclique.py`) |
| clique floor E7 | 35 (135 cliques, best W of dim 5 meeting all 135) | **35** (135; dim 5, all 135) | E7 as the E8 roots orthogonal to (½)⁸ (cx-medium used e1 + e2) |
| E7 rigidity: tangent dimension at the root decoder | 111 | **111** (rank 771 of 882; singular values jump from 1.41 to 4e-15) | own embedding and Hadamard basis (`calc/vrigid.py`) |
| U8X0 at block 25, ideal banks | 7.3057223e-4 | **7.3057223e-4** (moment root the same; first order 7.313e-4) | fold re-typed from its statement; histogram rebuilt from the certificate's class split, not from cx-medium's constant (`calc/vprice.py`) |
| T74 gauge44 at block 25 | 7.3920199e-4 | **7.3920199e-4** | same tool |
| helpers free | 2.3338e-3 | **2.3338e-3** | |
| most generous ceiling (data in {8, 16}) | 8.2607e-3 | **8.2607e-3** | (cx-medium's `items.out` prints 8.290e-3, a first-order D/N of the same layout; its NOTE cites the fold's 8.26e-3, the right figure) |
| A5 call on raw arrays, block 5 | exact; 60/60/30/181 = 331 | **exact** (max error 1e-14) for seeds 0 to 2; seed 2 gives 60/60/30/181 = **331** | my own simulator (`calc/vraw.py`), see below |
| E7 call on raw arrays, block 7 | exact; 7,220 paid levels (best of 6 seeds) | **exact** (max error 8e-14) for seeds 0 and 3: 7,262 and 7,297 | a Colab machine, `vcxm-raw-e7/` |
| our checker on U8X0 | every check passes; certificate identical | **identical** `csx_checks.json`; certificate content sha256 ccbfe460…, the same in all three copies | rerun into `calc/csx-rerun2/` |
| U8X0 per-call classes | x 896, y 896, helpers 3,778, copies 64 | **the same**, read from the certificate | |

**My raw-array simulator** (`calc/vraw.py`) was written from the docstring's semantics, not from cx-medium's class:
- *The frame operator.* Φ_U = H diag(i^{wt(π_U x)}) H. The projection π_U is computed from an arbitrary basis through
  the inverse Gram matrix over F₂. It does not assume an orthonormal basis, and it fails on a degenerate frame.
- *Gates.* A gate compares frames as sets of vectors.
- *Moves.* Every move asserts nesting.
- *Cost.* A move is paid only when the array holds content.

Since Φ_U depends only on U, cx-medium's cache by subspace is sound.

## 2. The proofs, line by line

**The clique floor (cx-medium section 2): valid in its scope, and a little stronger than stated.**
- *The exact cliques.*
  - Ports with mutually orthogonal roots have non-orthogonal labels, so no side information passes between them (an2-smallnet 2.1).
  - Each port's own data also needs a descent.
  - So on a clique C the delivered map is diagonal and invertible, and it can come only from totals. That needs at least
    ω totals whose *gather* sets meet C.
- *The counting.* Σ_k T(S_k) ≥ ω·#cliques. A total costs at least dim span(labels of S_k), and T(S_k) ≤ T(W_k). Both
  steps are correct, and the minimum may be taken over spans of labels.
- *Robustness.* Only the gather side is counted, so the floor also survives two things the note did not discuss:
  - targets that mix with each other after the scatter, since the rank of (I + E)·R·G restricted to C is at most the
    number of totals meeting C;
  - mutable sources, since an array holding x_c content without a descent has a frame containing u_c.
- *Its real scope.* Totals that descend once to frame 0, which is the hub's accounting. The note states this.

**E7's symmetric-decoder bound (42): correct.**
- At most 36 lines lie in a hyperplane: the roots in a hyperplane form a root system of rank at most 6, and the largest
  is E6.
- So each total gathers at least 27 lines. A 5-dimensional space holds at most 16 odd vectors, so each total costs at
  least 6.
- There are at least 7 totals, so the bound is 42.
- *A tighter remark* [proved here]. Among coordinate decoders (orthonormal bases), at most one total can cost 6. A
  cost-6 total needs f orthogonal to an E6, so f is a minuscule weight of E7. Two such weights always have inner
  product ±1/2 or ±3/2, never 0, so no two can be orthogonal. Coordinate decoders therefore cost at least 48. This only
  strengthens the note.

**The free-meeting bound: correct in its stated scope.** Copies that never climb, and distinct labels.
- At dimension i, the meetings are Σ_F t_F·s_F, with t_F ≤ 2^{k−i−1}, s_F ≤ 2^{i−1}, Σ t_F ≤ v and Σ s_F ≤ v.
- Moves of rank above 1 skip dimensions, which only lowers the count.
- 3,360 (57.7 %), 882 (43.8 %) and 90 reproduce.

**The piece bound behind 8.26e-3 (argued): plausible.**
- A data array gates in its middle call, so one cut falls in the levels 8 to 16, inclusive.
- r·ln(m/r) is concave, so the best split is at the ends of that range: {8, 16}.

**A scope note on the off-template blocks.**
- The note excludes E7 at 13 and the 112-port design at 17 by citing an2-smallnet's window theorem. That theorem assumes
  condition A: each side-information copy is erased against its own source.
- cx-medium's own explicit calls do not satisfy it. They drop their copies and carriers non-zero: 121 to 126 scratch
  arrays in the A5 call and about 2,000 in the E7 call (`calc/vraw_*.json`).
- That is legal, since it is the hub's own accounting for centre copies. But it means the exclusion is cited for a
  class that does not contain cx-medium's own networks.
- I did not find a way to exploit this. In the three-shear swap, call 3's side information must carry −y and so is
  copied after call 2, which probably forces a fresh window anyway. So this is a limit of scope, not a hole.

## 3. Refuted: "no idle levels" as the cause of ×0.44

**What cx-medium claims** (NOTE "In short" 4, section 3, section 9.1; the claim's obstruction 1 and next step 4):
- One pair leaves no idle levels, so helpers cannot be hosted, "which is the hub's only way to make junk cleaning
  free".
- They go on banks and carry 69 % of the narrow mass.
- The one design that saves is "held to ×0.44 by junk cleaning without idle levels".

**What hosting actually does** [computed with cx-medium's own `shrink.py` model, which it calibrated on 1.6647327e-3;
`calc/vhost.py`, `calc/vhost.out`].
- A hosted helper's in-call pieces are paid pieces in the histogram, exactly as on a bank.
- Hosting only lets those pieces *replace* idle pieces of data arrays that would otherwise be paid in large blocks.
- Ideal banks are D-neutral, so at first order (a ≈ D/N) a bank costs exactly its helper pieces' narrow mass, which
  hosting would pay too.

| pairs, block | hosted (cx-medium's exact fill) | every helper on ideal banks | hosting worth | helper share of narrow mass (hosted) | helpers free |
|---|---:|---:|---:|---:|---:|
| 1, 25 | 7.3057e-4 | 7.3057e-4 | ×1.000 (nothing to replace) | 68.8 % | 2.334e-3 |
| 2, 41 | 1.4667e-3 | 1.3233e-3 | ×1.108 | 68.3 % | 3.440e-3 |
| 3, 57 | 1.6523e-3 | 1.4054e-3 | ×1.176 | 68.1 % | |
| 4, 73 (the hub's layout) | 1.6647e-3 | 1.3929e-3 | ×1.195 | **66.1 %** | 3.109e-3 |
| 5, 89 | 1.6493e-3 | 1.3559e-3 | ×1.216 | 64.9 % | |

**So:**
- The hub does not make junk cleaning free. Its hosted helpers carry two thirds of its narrow mass, just as the banks do
  at one pair.
- Hosting is worth ×1.2 at the hub's layout. At one pair, perfect hosting would change nothing.
- The ×0.435 gap splits into two factors:
  - ×0.84: the hosting gain the hub has (1/1.195);
  - ×0.52: the centre's share of each call. One pair keeps 224 − 192 = 32 free directions per three calls, 10.7 per
    call. Four pairs keep 320 per nine calls, 35.6 per call. The helper mass per call is the same in both.
- Even with junk cleaning free, one pair gives only ×0.75 of four pairs (2.33e-3 against 3.11e-3).
- The note's own "2.33e-3, ×1.39 over the hub" compares free junk cleaning at one pair with the hub's *paid* junk
  cleaning.

**Where this leaves the note's figures.**
- *What holds.*
  - The prices.
  - The requirement that the helpers' narrow mass fall to 17.9 % of U8X0's.
  - The conclusion that P112 at 25 does not beat the hub.
- *What does not hold.*
  - The diagnosis. The lever at medium blocks is the centre's share per call, which comes from too few pairs per centre.
    It is not the absence of idle levels.
  - Section 3 item 1, "each [bank] costs a whole array of narrow mass". A bank's narrow mass at one pair is its helper
    pieces' mass, which hosting would pay too.

**A smaller slip in obstruction 2.** It says that only the hub's 112-port decoder clears the one-pair break-even. By
cx-medium's own table (`calc/items.out`), the coordinate decoders also clear it:
- E8 with 120 ports: 72 against 80, D = +24;
- the 112-port design: 72 against 74.7, D = +8.

## 4. Refuted: "five template points, one saves", and "medium blocks closed"

**The missing points.**
- With one pair the block is 1 + 3(k − 1). cx-medium listed only odd windows (5, 7, 9).
- Root systems do have odd windows, but sub-designs need not. The hub's own 112-port design is a sub-design of E8.
- Hyperplane sections of the label space give windows 6 and 8, so template points at blocks **16** and **22** [computed,
  `calc/vsub.py`, `vsub.out`].

**Block 16 closes.** E7's 32-line section has k = 6, ω = 4 and 120 maximum cliques. Its clique floor is 21.33 per call,
exactly the one-pair break-even, so D ≤ 0 [computed, `calc/vsub2.py`].

**Block 22 does not close.** The design is E8's 56-line section: the 56 lines of E8 at 60° to one root ρ, that is, E8
minus E7 minus ρ (`calc/vsub_identify.out`). Its figures [computed, `calc/vsub3.out`, `calc/vfloor8.out`]:
- n = 56, k = 8, ω = 4, 630 maximum cliques;
- 28 junk partners and 27 orthogonal partners per port, σ = 7, so block 22 = 1 + 3·7;
- 1,568 junk entries per call, against 5,824 for the 112-port design and 2,016 for E7.

| | E8 56-line section, block 22 | E7, block 19 (for comparison) |
|---|---:|---:|
| one-pair break-even (per call) | 37.33 | 42 |
| clique floor (per call) | **28.0** (every one of the 417,198 subspaces of F₂⁸; best W of dim 6 meeting 540 of 630 cliques; calibration: the same code gives E7's 35) | 35 |
| best decoder known | 60 (coordinate decoders over 202 orthonormal bases) | 49 |
| D per replica at the floor | +28 (0.25 per data array; the 112-port design has 0.14) | +21 |
| ceiling, junk free, data in {σ, 2σ}, at the floor | 1.6e-2 | 1.2e-2 |
| the same with U8X0's helper mass per junk entry on ideal banks | 2.95e-3 | 1.90e-3 |

- In the same generous model the 112-port design gives 9.4e-4 against its real 7.3e-4, so these ceilings run about
  ×1.3 high (`calc/vceil22.py`, `vceil22.out`).
- This is not a candidate. No decoder below 60 is known, so D is −68 today, and the 28 is only a floor.
- But it is a template point inside the range that cx-medium's tools cannot exclude, and its margin is wider than
  E7's:
  - floor 25 % below break-even, against E7's 17 %;
  - a quarter of the 112-port design's junk.
- So these three statements are not supported:
  - "the hub's template has five points [in 13 to 25], and only one of them saves";
  - "Medium blocks inside the hub's template should be treated as closed for the record";
  - "The only openings left in this range are new centres [for E7 or the 112-port design]".

## 5. Things the brief asked for that the note does not have

- **The refund (BRIEF3 notice of 12:55).** The notice asked for refund through shear copies to be itemised separately
  in the price breakdown. The note's breakdown does not do so. The prices do not change, but the breakdown the notice
  asked for is missing.
- **The height filter (notice of 13:20).** It arrived as cx-medium was finishing. The only word priced, U8X0, is
  Lean-checked, so the filter does not matter here.

## 6. Predictions (STATUS.txt, written 13:24) against outcomes

| prediction | outcome |
|---|---|
| V1 (90 %): my labels give the same clique counts and floors | held; I built no labels, which is stronger |
| V2 (85 %): the clique floor is valid in scope | held; it is also robust to target mixing and mutable sources |
| V3 (80 %): the coordinate decoders and the 42 bound reproduce | held. I did not rebuild A5's and D6's 30 separately; the floors and rigidity took the time |
| V4 (70 %): the price lands within 2 % | held; it matches to every digit |
| V5 (70 %): an own-simulator replay of A5 is exact | held, for A5 and E7 |
| V6 (65 %): the free-meeting bound holds in scope | held |
| V7 (60 %): headline not refuted; a scope limit or slip found (35 %: material) | held for the headline. The material finding came in two parts: the wrong diagnosis of ×0.44, and the missed block-22 point |

## Files

All in `checks/wht43/anatomy/v-cx-medium/`:
- `STATUS.txt`, `RUNLOG.txt`, this `NOTE.md`.
- `calc/vclique.py` → `vclique.out`, `vclique.json`: label-free clique floor (A4, A5, D5, D6, E7).
- `calc/vprice.py` → `vprice.out`, `vprice.json`: fold re-typed; block-25 prices, helpers free, ceiling.
- `calc/vraw.py` → `vraw_A5_s{0,1,2}.json`: own raw-array simulator; E7 replay in `vcxm-raw-e7/vraw_E7.out`.
- `calc/vrigid.py` → `vrigid.out`: E7 tangent dimension.
- `calc/vhost.py` → `vhost.out`: hosting against banks, helper shares, helpers free.
- `calc/vsub.py` → `vsub.out`; `calc/vsub2.py` → `vsub2_E7-32.json`; `calc/vsub3.py` → `vsub3.out`;
  `calc/vfloor8.py` → `vfloor8.out`; `calc/vsub_identify.out`; `calc/vceil22.py` → `vceil22.out`: the even-window
  points.
- `calc/csx-rerun2/`: our checker rerun (`run.log`, `csx_checks.json`, certificate).
- `vcxm-sub/`: a Colab copy of the block-22 floor job. I stopped it once `vfloor8.py` finished the same computation in
  12 s on the Mac.

## Sources

- `../BRIEF3.md`, `../BRIEF.md`, `../BRIEF2.md`, `../NOTE.md` (rounds 1 and 2).
- `../an2-smallnet/NOTE.md` (sections 2.1 to 2.6, 3.2).
- `../cx-medium/` (NOTE, STATUS, RUNLOG, and calc: `designs_cx.py`, `cliquefloor.py`, `e7call.py`, `price_p112.py`,
  `price_t74_p1.py`, `ceilings.py`, `shrink.py`, `items.out`, `csx-rerun/`).
- `../an-ours/anatomy_ours.py`: the fold's statement and the Decimal root used for the cross-check.
- `bound-check/fsp-cc/tools/csx_check.py` and `bound-check/fsp-cc/out/`: our checker and the stored certificate.
