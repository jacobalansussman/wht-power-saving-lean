# hx-price NOTE: A-fc's full-cut price, re-derived independently

Round wht45/hubx, agent hx-price, 20 Oct 2026, 04:36-04:50 EDT. Not public text. All work was light Python on the Mac
(each run under 1 minute and 0.3 GB). No Lean, no Colab, no hub code. Nothing posted. Commands: `COMMANDS.txt`. Outputs: `out/`.
I did not use fc-reopt's pricing code (`tools/price_one.py`) and did not read it. Labels: [computed], [argued], [read].

## 1. Verdict: the price reproduces. STOP not written

**A-fc full cut b = 1.7887046e-3** [computed, three routes]:
- route A, `pricer.py`'s exact rational enclosures, certified on the 1e-24 grid: 0.00178870456430307919393;
- route B, `mixed.py` at 60 digits: 0.0017887045643030791939300419693164749977017732938491302815;
- route C, my own exact-rational interval arithmetic, 70 digits: 0.0017887045643030791939300419693164749977017732938491302814669639.

Routes B and C agree to 3e-60. All three equal fc-reopt's 0.00178870456430307919393.

**Grid numerator N = 17887045** (the root sits 0.643 of a grid step above it) [computed, certified by routes A and C].
- N/1e10 is below the root, and (N+1)/1e10 is above it.
- The Lean stack's convention (fc-rate NOTE section 2) also holds: the rate point (100N+2)/1e12 is below the root, the
  seam (100N+1)/1e12 is below it, and the rate point for N+1, (100(N+1)+2)/1e12, is above it. So a generator must
  refuse 17887046.

Margins [computed]:

| compared with | difference | relative |
|---|---|---|
| design B full cut, 1.7806935e-3 (the Lean result) | +8.011e-6 | +0.450 % |
| design A full cut, 1.7797095e-3 | +8.995e-6 | +0.505 % |
| #3945's hub b, 1.7306632e-3 | +5.804e-5 | +3.354 % |
| #2979's hub b, 1.7290340e-3 | +5.967e-5 | +3.451 % |
| A-fc's own uniform b, 1.6897030e-3 | +9.900e-5 | +5.859 % |

## 2. Heads [computed]

There are **16,576 literal heads**: 8,288 on child 8 and 8,288 on child 10. Three counters agree:
- My own code, `tools/hx_heads.py`, written for this check. It takes pairs and sides from the template's idle sweeps
  (checked against PLAN), role widths from #2979's published PROVIDERS.json (checked against my own GF(2) ranks of the
  words, and equal), and builds the chronology from one global event sequence. It also asserts:
  - the FILL schema and canonical row order;
  - the compact prefix and capacity of every bin;
  - every role comes before its bin's idle sweep;
  - at most one head per array.
- lever-reprice `heads_fill.py`, run on pkg2979 with A-fc's FILL: 16,576 HEAD-literal.
- lever-reprice `recount.py`: 16,576 fully idle rank-48 bins on children 8 and 10, and none on children 12 and 13.

Controls:
- The negative control (the bridge that reads pair 3 moved ahead of the sweep) gives 0 heads in my code and in heads_fill.
- On #2979's published fill, my code gives 16,086 heads (the known value) and a profile equal to its PROFILE.json.
- On A-fc, my profile equals A-fc's PROFILE.json digit for digit.
- The rank mass is 4,816,512 = 65,088 x 74.

A-fc's ORDER and DISPATCH have the same content as #2979's published ORDER and DISPATCH. The template, plan, words and
providers files match #2979's pins.

## 3. The full cut, built from the counts [computed]

The paid profile is 666 calls x the local word {1: 2,445, 2: 706, 3: 232, 4: 108, 5: 18, 8: 8}, plus the idle pieces
{1: 2,269, 6: 5, 7: 3, 17: 8, 19: 1, 23: 25,905, 48: 16,576}.

- **Heads.** All 16,576 rank-48 idle pieces are heads. They are removed, and their arrays get target 25 = 73 - 48.
- **Other arrays.** Every other array keeps target 73. Targets: **{73: 49,728, 25: 16,576}**, W = 66,304.
- **Split.** No other idle piece is above 24, so nothing is split.
- **Children** (2,370,513 paid blocks): {1: 1,630,639, 2: 470,196, 3: 154,512, 4: 71,928, 5: 11,988, 6: 5, 7: 3,
  8: 5,328, 17: 8, 19: 1, 23: 25,905}. The rank mass drops by exactly 48 x 16,576 (4,816,512 to 4,020,864).
- **Every paid rank is at most 23 <= 24 < 25.** So the rule is separated, and it fits the Lean rule shape
  `Rule.twoTarget 73 25 W H P 24` (paid ranks in [1, 24]).

## 4. Theorem 1G's hypotheses for this job model

| hypothesis (rl-prove/PROOF.md) | A-fc | status |
|---|---|---|
| Targets W_k > 0, paid ranks r_j >= 1, tau the root of sum r_j^tau = sum W_k k^tau | {73: 49,728, 25: 16,576}; ranks 1..23; the root above | [computed] |
| (S) separation r_max < k_min | 23 < 25 (also <= 24 < 25 under the Lean cap) | [computed] |
| theta = ln r_max / ln k_min < 1 (Lemma D exits) | 0.9741 (0.9873 with cap 24) | [computed] |
| kappa finite (log factor only; exponent tau + o(1)) | eta_m = 0.99059, kappa = 27 (56 with cap 24, same as design B) | [computed] |
| Each removed block is a literal head: rank 48, the first operation on its own array, acting on that array alone, at most one per array, unhosted | all 16,576, by three counters; asserted one per array; 0 under the negative control | [computed] |
| The truncated word is correct for the reduced targets (73 and 25) at every column count that occurs | Same template (21297b20) and head position (pair 3 idle sweep from level -1, read first by bridge 15) as design B. On design B, M1's `cut_heads_full` proves this in Lean; for A-fc it is not proved | [argued] (Lean on design B only) |
| Splitting idle pieces above 24 into pieces of 24 is a valid job split | not used: no piece above 24 | n/a |
| The unit is one whole word (the 74-replica T74 word, 66,304 arrays); hosted roles stay inside it; paid profile = the word's actual blocks | the profile is recomputed from the FILL and the words and equals PROFILE.json; the cut blocks have no hosted role | [computed] |
| Finite free passes F per unit | the word's translations, bridges and adapters; they only change Theorem 1G's constant | [argued] |

A caveat to state, not a failure: the hub engine asserts m^2 > 2 r_max^2. This holds with m = 73 (5,329 > 1,058). It
would fail with k_min = 25 in place of m (625 < 1,058). Theorem 1G does not use this condition, and design B has the
same situation (r_max 24 against target 25) yet is proved in Lean.

## 5. Controls [computed]

Each control was priced by all three routes, and the grid facts hold for each.

| control | expected | got |
|---|---|---|
| design B full cut (#3945 profile, 80,640 heads, 57,312 other 48s split 24 + 24) | 1.7806935e-3, N = 17806935 | 0.001780693540576287216656, N = 17806935; 11,641,032 paid blocks, as in fc-rate |
| design A full cut (#2979 profile, 16,576 heads, 11,614 other 48s split) | 1.7797095e-3 | 0.001779709547136452622863, N = 17797095 |
| A-fc uniform b (hub engine re-implementation, 1e-24) | 1.6897030e-3 | 0.001689703011137435213783 |
| #3945 uniform by the hub engine | hub-stated 0.001730663236534549309146 | equal, digit for digit |
| #2979 uniform by the hub engine | CLAIM complex_lower 1729033982646513517717/1e24 | equal, digit for digit |
| half cuts, design B / design A | 1.7581520e-3 / 1.7564708e-3 | 1.7581520168e-3 / 1.7564708007e-3 |

## 6. Not checked

- The upper bounds 1.7887182e-3 (T74 clock, Hall LP) and 1.7893658e-3 (any fill of #2979's words). They come from
  200 LP solves with Lagrangian certificates, so checking them is not quick. A-fc is below both, which is consistent.
- Hub admission and owner-engine checks of A-fc as a fill, including the 9 spare roles. These belong to the
  admission check (hx-admit), not to me.

## 7. Files

- Tools: `tools/hx_heads.py`, `tools/run_heads.sh`, `tools/hx_price.py`, `tools/hx_1g.py`.
- Outputs:
  - `out/hx_heads-Afc.json` (+ `-control`);
  - `out/hx_heads-2979pub.json` (+ `-control`);
  - `out/heads_fill-Afc.json` (+ `-control`);
  - `out/recount-Afc.txt`, `out/recount-2979pub.txt`;
  - `out/price-Afc.json`, `out/price-controls.json`, `out/thm1g-constants.json`.
