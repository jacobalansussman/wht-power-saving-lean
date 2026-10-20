> Working note written by an AI agent during the round of 13 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# Open: does the full E8 unit need 72 copy steps? (parked, 13 October 2026)

Parked on 13 October 2026. Nothing is running on it.

## The statement

For the E8 complex unit with all 120 ports (P120), every decoder needs copy cost cst >= 72.

## Why it matters

- Monotonicity (agent d-m-hubopen-1, 13 Oct, own argument, unreviewed): adding one port raises the copy cost by
  at most 1, so cost(P') <= cost(P) + |P' \ P|. Hence cst(P120) >= 72 implies cst(P_v) >= v - 48 for every
  punctured unit. All known decoders sit exactly there: P104 at 56, P111 at 63, P112 at 64 (the hub's #799
  complex side), P120 at 72; bit side P111 at 63, P120 at 72.
- So a proof would close the hub's "open problem 2" (beating the star floor) for the whole puncture family, and
  tell everyone to stop searching decoders of this shape. It is a closure, not a gain.

## What is known

- Upper bound 72: attained (nine totals of span 8).
- Lower bound 55 for every decoder on P120 (ours: `wht35/e8/s-copies.md`); 53 on P112 (`u-ours` in the copy-wall
  workflow, octad count over all subspaces of F_2^9).
- Proved 72 (P120) and 64 (P112) for decoders whose totals are linear in the root vectors (the label family):
  exhaustive over subspaces (`u-ours`, copy-wall workflow run a workflow run).
- With the hub's helpers fixed, every star decoder with 9 or 10 totals on P112 costs at least 64
  (d-m-hubopen-1, exact runs; files `../bound-check/g-p112-decoder/`).
- SAT: P111 at budget 62 UNSAT in 85 s in d-m-hubopen-1's model (narrower than the general question); the hub's
  own SAT runs on the general case were killed, its flip-walks (~1.5e8 steps) found nothing below.
- The difficulty is the non-linear totals (totals outside the label row space, mixed ranks). Lower bounds of this
  kind are minrank-like, NP-hard in general.

## Two routes if anyone picks it up

1. Computer proof: SAT/ILP for "a decoder of cost <= 71 exists on P120", with symmetry breaking under W(E8)
   (order 696,729,600). A Colab with 8 cores / 50 GB could run it overnight.
2. Hand proof that lifts the octad count from 55 to 72 using E8 structure; needs a new idea.

## Sources

- Agent d-m-hubopen-1, workflow `build-on-community` (run a workflow run), files in `../bound-check/g-p112-decoder/`.
- Copy-wall workflow (run a workflow run): u-floor (`u-floor/copyfloor.py`, the hub's floors reproduced), u-ours.
- Hub: `research/format-limits-corrected/` (A.6 star floor, C.1-C.6), `research/e8-decoder-floors-and-cold14-notes/`,
  PR #717 (DreamingOfClouds; its D1 open item).
