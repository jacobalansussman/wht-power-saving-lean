# Proofs of this section

These are the working proof files of the four rounds of 15 October 2026, as the AI agents wrote them, with private
file paths removed. Each file starts with a header that gives its status now, after the skeptics' checks. Where the
header and the file disagree, the header (the later verdict) wins. No person has read these proofs, and nothing in
them is checked in Lean. "I" inside a file is the agent that wrote it, not me.

| file | what it proves or checks | status now |
|---|---|---|
| [ledger-for-every-unit.md](ledger-for-every-unit.md) | the ledger M = W·m − 2N + X + 2D + Z + B − P for every unit; the anatomy's identity as a special case; parity; Theorem 4.3 (two arrays need packing, via Hu's theorem); Corollary 4.4 (clustered sources or sinks); the any-W bound; the network-coding reading | argued by hand, not checked; one skeptic upheld it with corrections |
| [common-frame-accounting.md](common-frame-accounting.md) | exact accounting for common-frame designs (Theorem 1); Corollaries 2 and 2′ (two arrays need 2m); the general bound (Theorem 5); the median bound | argued by hand, not checked; one skeptic upheld it with one correction |
| [slack-skeptic.md](slack-skeptic.md) | the skeptic's record of the slack round: definitions, the disagreement over the Lean savers, parity, U_m, ceilings | a skeptic's check |
| [packing-capacity-and-rigidity.md](packing-capacity-and-rigidity.md) | Theorems A and A′ (what packing buys), M ≥ 2m − P*, Proposition D, Lemma S, rigidity, the refund reading of the certificates | argued by hand, not checked; one skeptic upheld it with corrections (one of which is itself corrected in the header) |
| [frame-game-sandwich.md](frame-game-sandwich.md) | the frame game; the lift at f = 1 mod 4; the walker and landing lemmas; the position test | argued by hand, not checked; one skeptic upheld the theorems and refuted the "closest design" |
| [frame-game-sandwich-skeptic.md](frame-game-sandwich-skeptic.md) | the skeptic's record of the build track, with the f = 3 mod 4 lift | a skeptic's check |
| [frame-game-completion.md](frame-game-completion.md) | Theorem 1: every frame-game program completes to an engine unit at every block | argued by hand, not checked; one skeptic upheld it (wording only) |
| [closed-patterns-and-conjecture-tc.md](closed-patterns-and-conjecture-tc.md) | the parity law (Lemma 1), "one shared pair gives one ratio", why orbits multiply junk, Conjecture TC | argued by hand, not checked; one skeptic upheld it with six corrections; Conjecture TC open |
| [heights-theorem-c.md](heights-theorem-c.md) | Lemma M, Theorem A, Theorem C (heights) with Proposition D, Theorem B, Corollary T₀, Conjecture T; Propositions S and H3 refuted | argued by hand, not checked; one skeptic upheld it with four corrections |
| [single-term-game.md](single-term-game.md) | where a counterexample to Proposition 6 can live; the single-term game (section 4: Lemmas 4.1 to 4.7, Propositions 4.4 and 4.5) | section 4 upheld by one skeptic with corrections; sections 1 to 3 not checked by a skeptic |
| [single-term-game-skeptic.md](single-term-game-skeptic.md) | the skeptic's record of the single-term game, with the independent block-4 cube run | a skeptic's check |

## Working names used inside the files

The files cite each other and their programs by the rounds' working names:

| name inside the files | here |
|---|---|
| `sl-combinatorics`, "sl-comb" | [ledger-for-every-unit.md](ledger-for-every-unit.md), [../programs/ledger/](../programs/ledger/) |
| `sl-algebra`, "sl-alg" | [common-frame-accounting.md](common-frame-accounting.md), [../programs/common-frame-accounting/](../programs/common-frame-accounting/) |
| `v-slack` | [slack-skeptic.md](slack-skeptic.md), [../programs/skeptic-checks/slack/](../programs/skeptic-checks/slack/) |
| `fp-measure` | [../programs/packing-measurement/](../programs/packing-measurement/) (no proof file; the refund identity is written out in the section README, section 4) |
| `fp-capacity` | [packing-capacity-and-rigidity.md](packing-capacity-and-rigidity.md), [../programs/packing-capacity/](../programs/packing-capacity/) |
| `fp-build` | [frame-game-sandwich.md](frame-game-sandwich.md), [../programs/frame-game/](../programs/frame-game/) |
| `v-fp-build` | [frame-game-sandwich-skeptic.md](frame-game-sandwich-skeptic.md), [../programs/skeptic-checks/frame-game/](../programs/skeptic-checks/frame-game/) |
| `pk2-complete` | [frame-game-completion.md](frame-game-completion.md), [../programs/frame-game-completion/](../programs/frame-game-completion/) |
| `pk2-orbits` | [closed-patterns-and-conjecture-tc.md](closed-patterns-and-conjecture-tc.md), [../programs/orbit-designs/](../programs/orbit-designs/) |
| `p6-prove` | [heights-theorem-c.md](heights-theorem-c.md), [../programs/heights/](../programs/heights/) |
| `p6-refute` | [single-term-game.md](single-term-game.md), [../programs/single-term-game/](../programs/single-term-game/); its `job10` to `job19` are `../data/single-term-game/run10` to `run19` |
| `v-p6-refute` | [single-term-game-skeptic.md](single-term-game-skeptic.md), [../programs/skeptic-checks/single-term-game/](../programs/skeptic-checks/single-term-game/); `vp6r_b4` is `../data/skeptic-checks/single-term-game/block4-cube-2in7/` |
| `anatomy`, `an-openai`, `an-ours`, `an-synth` | the anatomy round, in [../../why-the-savers-save/](../../why-the-savers-save/) |
| `nc-gain` | the network-coding track, in [../../other-fields/](../../other-fields/) |
| `exist2`, `ex2-certnum`, `meta/MAP.md`, `moonshot`, `o3`, `theory-guide/RULES.md`, `proofs/ladder` | earlier rounds, in the other sections of [../../](../../) |

Some files mention outputs that are not shipped, such as run logs, the briefs and machine job folders. Their
results are in the shipped output files named next to them.
