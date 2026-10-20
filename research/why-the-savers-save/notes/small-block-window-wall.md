> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# an2-smallnet: proposal 3, the small-block crossing network: no network, and the bar not reached

Anatomy round 2, track an2-smallnet, 15 October 2026, 10:39 to 11:47 EDT (Mac clock). Task (BRIEF2.md): build or
exclude an-synth's proposal 3, a crossing network at block 8 to 16 whose stages share window coordinates, which an-synth
priced at a = 0.096 at block 16 if it existed (two free directions per array, free bridges, children of rank 4).

Labels: **[proved]** a complete argument written here (checked by me only; no skeptic yet); **[computed]** a program in
`calc/` and its output; **[argued]** reasoning with a gap that is named; **[read]** from a named file. All programs are
our own, run on the Mac. No Colab was used, no third-party code was run, nothing was posted.

## In short

1. **Proposal 3's premise fails where it matters: an array's source calls cannot share window coordinates.** In the
   model proposal 3 was posed in (frames only grow, except on descending bridges; one common frame per gate; no
   per-array address maps), I prove a per-array form of the stage–window wall:
   - a target must take up, in every call, the labels of all its junk partners (in mirrored form for calls x −= y)
     [proved];
   - an array that is a source in two amortised calls (calls whose bridge is shared by many pairs) must take up two
     windows that are orthogonal to each other. So the block is at least **1 + 2σ**, where σ is the window a call
     costs: σ = k − 1 for every design checked (8 for E8 and the hub's P112, 6 for E7, 4 for D5, D6, A4, A5, and 2 for
     D4, A3)
     [proved, given that each side-information copy is erased against its own source, which is how every network we
     know works; an exact count closes the gate-order gap for each design's natural code];
   - in any network of bank-to-bank calls, some bank is a source in at least two calls and some bank takes part in at
     least three [proved for every number of banks; also checked by exhaustive search up to 5 banks];
   - a backward call and the forward call after it cannot share a window either [proved, under condition A]; with
     a counting step that I outline but have not carried out, the swap's third call costs a fresh window, and the
     block is at least **1 + 3σ**: 25 for E8 and P112, 19 for E7 (sec. 2.6).
2. **Which designs fit in block 16.**
   - D4 and A3 save nothing under *any* centre [proved: a new floor on the centre charge, valid for every
     decomposition of the index code].
   - D5, A4 and A5 save nothing at the hub's own centre accounting (at least one level per port per call). D6 would
     save at the hub-type charge 0.8, which is unproved. For these four, the floor I can prove is only 0.4 to 0.6.
   - E7 (σ = 6), and also D8 and the non-saving E6, D7, A6, A7, fit at the proved floor 1 + 2σ = 13, but not at 19,
     the floor if the third call needs a fresh window. E8 and P112 need 17 and 25.
3. **Ceilings, solved exactly from the moment equation, under deliberately generous assumptions** (junk cleaning free,
   every in-call climb one child, the hub's centre charge) [computed, `calc/ceiling.out`]:

   | design, block | one free direction per array | 1.5 (the build's ceiling for a second line) | 2 (second line free) |
   |---|---:|---:|---:|
   | E7, block 13 (proved floor) | 0.034 | 0.086 | 0.140 |
   | E7, block 19 (third window fresh) | 0.015 | 0.039 | 0.063 |
   | P112, block 17 (proved floor) | 0.034 | 0.075 | 0.117 |
   | P112, block 25 (third window fresh) | 0.015 | 0.033 | 0.052 |
   | D6, block 13, hub-type charge assumed | 0.013 | 0.047 | 0.082 |
   | largest over D4..A5 at 1 + 3σ (block 7 to 14), proved centre floors only | 0.040 | 0.065 | 0.097 |
   | largest over D4..A5 at the proved floor (block 5 to 10), best charge known (proved for D4, A3; hub-type for the rest) | 0.028 | 0.079 | 0.166 |
   | largest over D4..A5 at the proved floor (block 5 to 10), proved centre floors only | 0.085 | 0.127 | 0.190 |

   - With one free direction per array (the only case anyone has built), nothing exceeds 0.04, except the small root
     designs at the proved floor (blocks 5 to 9) when only their weak proved centre floors are used (up to 0.085);
     with the best charge known for each design they stay below 0.03.
   - Every entry at or above 0.1 needs at least two of the following at once:
     - a third call that shares coordinates (excluded outright, and for partial overlaps up to a count, sec. 2.6);
     - a second free line per array at no cost (argued against, sec. 4; the build refuted it on our word);
     - a centre charge for D5, D6, A4, A5 at the bare proved floor, 0.4 to 0.6 instead of the 0.8 to 1.6 every known
       code pays.
   - Free junk cleaning at a block with no idle levels is assumed in every entry (argued against, sec. 5).
4. **Verdict: excluded for the bar, with a named gap.**
   - No network exists in the form proposal 3 assumed: windows shared between an array's source calls [proved].
   - Every reading that is not built on at least two unproved favourable assumptions stays below 0.09, and every
     one-sided reading stays below 0.09.
   - What a full proof still needs (sec. 9):
     - the count that finishes the third window (2.6), or a centre floor of about 0.8 per port per call for D5, D6,
       A4 and A5;
     - the cost of a second line (sec. 4);
     - the cost of junk cleaning without idle levels (sec. 5).
   - Together with the anatomy's other findings, no reading of the crossing class reaches the bar at any block
     without stacking such assumptions. The bar needs a mechanism that is not crossing slack.
5. **One loophole is real and named:** per-array address maps (gauges), which the hub's words use between calls,
   move frames for free. My lemmas hold only in a fixed gauge. The class in which proposal 3 was posed has no gauges,
   so the verdict on proposal 3 stands. But "the crossing class with gauges is closed" is not proved here (sec. 2.7).
6. **Raw-array checks pass** [computed]: the frame calculus identities at blocks 1 to 5; a lone pair costs exactly 2m;
   a one-pair three-shear swap through scratch bridges costs 2m + 4 and is exact; two pairs sharing one bridge, with
   their junk removed by side information, end exactly with y + x, with scratch erased, and each target's frame holds
   its junk partner's label (the lemma of sec. 2.2, seen on arrays).

## Words used here

- **Frame, label, port.** Each array carries a subspace of F₂^m, its frame; a paid move climbs or drops one level. A
  data array X_a starts on a line ⟨u_a⟩ (its port label) and ends on the whole space; its partner Y_a starts at 0 and
  ends on u_a^⊥. Its data must cross from frames that contain u_a to frames inside u_a^⊥, and only a descent can do
  that (an-openai 2.3).
- **Call.** One transfer from a source bank to a target bank (y += x or x −= y on every pair of the call), through a
  shared bridge (the "centre": t arrays that gather, descend, scatter). An **amortised** call has fewer bridge arrays
  than pairs; then each target also receives junk (other sources' data), which side information removes.
- **Side information, junk partner.** For target a, a junk partner b is a source whose data the bridge delivers to a by
  mistake. It is removed by a copy of x_b that reaches Y_a without any descent. The ports of a and b must have
  orthogonal labels (an-synth sec. 3).
- **σ (sigma).** The F₂ rank of the labels of a port's junk partners. For every design here, σ = k − 1, where k is the
  window (the label space): the junk partners' labels span the whole orthogonal complement of the port inside the
  window.
- **Gauge.** A free per-array address map that moves an array's frame. Not allowed in proposal 3's model.

## 1. What proposal 3 needed

an-synth (sec. 4.3) priced a network at block 16 with two free directions per array, bridges free and children of rank
4: a = 0.096. Our code reproduces this: 0.0963 [computed, `calc/ceiling.py`, calibration line]. Two things were
assumed, and each is tested below:
- that consecutive calls can reuse the same window coordinates, so that a call does not cost every array fresh levels
  (sec. 2);
- that the bridges and the junk cleaning cost little at such small blocks (sec. 3 and 5).

## 2. The stage–window wall, per array

All statements in this section are in a fixed gauge: frames only grow, except where an array descends, and every gate
joins arrays that share one frame. Data arrays never descend (a data array that descends is a bridge in its own right
and pays 2 moves per level). So X_a's frame always contains u_a, and Y_a's frame always lies inside u_a^⊥.

Sections 2.1 to 2.4 are written for calls y += x, whose sources hold their port lines. A call x −= y is the
time-reversed mirror of such a call: reverse time, invert every gate and replace every frame by its orthogonal
complement. That gives a legal network of the same cost, in which descents stay descents and X and Y swap roles. So
every statement holds for x −= y calls with the roles swapped, and every block floor is the same.

### 2.1 Side information leaves a source only while the source is orthogonal to the target's port [proved]

Take a copy of x_b that reaches target Y_a with no descent on the way. Along a route with no descent, frames only grow:
each array's frame grows between gates, and each gate equalises frames. So the source's frame at the moment it sends
the copy lies inside Y_a's frame when the copy arrives. That frame lies inside u_a^⊥ (Y_a's frame always does). So
**at the moment of sending, X_b's whole frame is orthogonal to u_a.** Since frames only grow, once X_b's frame contains
a vector that is not orthogonal to u_a, x_b can never again serve as side information for a without a new descent.

### 2.2 Every target takes up its junk partners' labels, in every call of this direction [proved]

By 2.1, Y_a's frame after the call contains every junk partner's frame at sending, and that frame contains the
partner's port label. So **Y_a's frame takes up span{u_b : b a junk partner of a}**, a space of dimension σ = k − 1.

### 2.3 Sources take up the window of their call [proved under condition A]

The copy of x_b that went to target a must itself be erased after use. **Condition A:** it is erased against X_b, the
only array that holds x_b alone during the call. That erasing gate is at a frame that contains Y_a's frame at the time
the copy arrived.
- Call the copy **late** if Y_a's frame then already contained S(a), the span of a's junk partners' labels (the
  orthogonal complement of u_a in the window). A late copy makes X_b take up S(a).
- Two late copies at targets with different ports make X_b take up S(a₁) + S(a₂), which is the whole window.
- An **early** copy arrives while Y_a's frame meets S(a) in a proper subspace. All earlier arrivals' labels lie in that
  subspace, so a target has at most f early arrivals, where f is the largest number of a port's partners whose labels
  lie in one hyperplane of S(a).
- Every source sends s copies (s = junk partners per port). So at most n·f/(s − 1) sources have one late copy or
  none.

`calc/early.py` computes f exactly [computed]:
- E8: f = 32 of s = 56, so at most 70 of 120 sources escape;
- P112: f = 32 of s = 52, at most 70 of 112;
- E7: f = 16 of s = 32, at most 33 of 63;
- A4, A5: at most 8 of 10 and 9 of 15;
- D5, D6: these have twin ports (e_i + e_j and e_i − e_j) with the same label, so a source needs three late copies to be
  sure of two different labels. The count n·f/(s − 2) gives at most 16 of 20 and 17 of 30.

E8, P112, E7, A4 and A5 have distinct labels [computed].

So **in every amortised call, at least one source (roughly half of them) takes up the whole window.** One array is
enough for 2.4, because the block is common to all arrays.

All three savers we know (OpenAI's, the hub's, ours) satisfy condition A: the hub's helpers are restored by running
the call backwards, which meets the sources again. A network that erased its copies only against other arrays would
move the taking-up onto those arrays. That is the remaining gap on the source side.

These counts are for each design's natural code (the line Gram, every partner at ±1/2). A code that kept each port's
junk partners inside a smaller subspace would cost a smaller window per call, but it would have fewer junk entries per
port and so need more bridges. That trade-off is not explored here.

### 2.4 Two amortised source calls need orthogonal windows: block ≥ 1 + 2σ [proved, given 2.3]

Let a bank be a source in calls s < s′, both amortised and of the same direction. 2.5 shows that such a bank exists in
every exchange between two banks, where all of a bank's source calls go the same way. Every array of the bank is
a source in both calls. By 2.3, at least one of them, X_b, takes up its whole window in call s. Let V_s be the span of
the ports of X_b's junk receivers in call s; V_s has dimension σ and lies inside u_b^⊥. Then:
- X_b's frame after call s contains V_s.
- By 2.1, every port of a call-s′ junk receiver is orthogonal to X_b's frame at sending, hence to V_s.
- So V_s′ ⊥ V_s.
- Both are non-degenerate: each is the orthogonal complement of u_b inside a non-degenerate window, which holds
  because σ = k − 1.
- Pairwise orthogonal non-degenerate subspaces are independent.

X_b's final frame, the whole space, therefore contains ⟨u_b⟩ ⊕ V_s ⊕ V_s′. So **m ≥ 1 + 2σ**, and with r amortised
source calls, m ≥ 1 + rσ.

This is the wall in the form proposal 3 had to break. A later call's ports must avoid what earlier calls left in the
frame, for every array, and not only through the template's layout.

### 2.5 Some bank is a source in two calls; some bank is in three [proved]

Banks exchange data by calls of the form "bank i += ±bank j". A bank saves only if it ends with exactly one other
bank's data, so its data must leave (it is a source at least once) and it must receive (it is a target at least once).

**Lemma (no exchange with two calls per bank).** Suppose every bank took part in at most two calls. Then each bank is a
source exactly once and a target exactly once.
- If its target call came first, the bank would hold its own data plus another's. Its later source call does not change
  it, so it would end holding its own data. That is not allowed.
- So every bank is a source first and a target second. But the very first call of the network has a target bank for
  which this is its first call: a contradiction.

So some bank takes part in three calls or more. Counting gives more:
- The calls have as many sources as targets in total.
- If every bank were a source exactly once, every bank would also be a target exactly once, and the lemma rules that out.
- So some bank is a source in at least two calls.

For networks of more than two banks with one-sided labels, the bank that is a source twice might send in both
directions; there I have not checked that 2.4 still applies [argued].

The search `calc/shears.py` confirms the lemma exhaustively for 2 to 5 banks [computed]. With three calls per bank,
2 and 4 banks can exchange (the three-shear swap); 3 banks cannot.

Together with 2.4, any layout in which that bank's two source calls are both amortised needs m ≥ 1 + 2σ.

### 2.6 The third call of a swap [proved for sharing a window outright; the general count outlined]

In the three-shear swap, X is a source in calls 1 and 3 (y += x) and a target in call 2 (x −= y), and Y the other way
round. 2.4 already separates calls 1 and 3. The question is whether call 2 can share its window with call 3.

At first I thought the lemmas allowed it, because during call 2 the X arrays take up only the Y-sources' frames. That
overlooked one fact.
- In call 2, the copies of y_b go up to their targets X_a and are erased against Y_b. That is condition A, mirrored:
  Y_b is the only array that holds y_b alone.
- The erasing gate is at a frame that contains X_a's frame, which holds u_a. So **after call 2, Y_b holds S₂(b)**, the
  span of its call-2 receivers' ports. This needs no gate-order count: every such copy carries its target's port
  line.
- The erasures must come before call 3 changes Y_b's content.

Now suppose call 3 uses the same window W.
- A call-3 bridge scatters first to some Y_{b₁}, at Y_{b₁}'s frame, which contains S₂(b₁) = W ∩ u_{b₁}^⊥. Every later
  target b of the same bridge is met at a larger frame, which must lie inside u_b^⊥.
- So u_b ⊥ S₂(b₁), which inside W leaves only u_b = u_{b₁}.
- With distinct labels (E7, E8, P112, A4, A5) each call-3 bridge serves one port, so call 3 is not amortised. With D4's
  repeated labels, it serves one label class, and the class floor then gives one level per port.

**So a backward call and the forward call after it cannot share a window** [proved, under condition A in both
directions].

For windows that only partly overlap, the same argument shows that every port a call-3 bridge serves after its first
target must be orthogonal to everything that target holds: S₁(b₁) from call 1 and S₂(b₁) from call 2.
- Take a Y_{b₁} that is also among the Y arrays whose frame starts call 2 inside W₂^⊥. That is the mirror of 2.3, which
  holds for most of them by the same count.
- Then S₁, S₂ and S₃ of that array are pairwise orthogonal, and the block is at least 1 + 3σ.

The step that remains is to count that such a Y_{b₁} exists. It needs (a) the mirrored count of 2.3 and (b) a count of
how many targets each bridge serves before the rest of its targets stop spanning the window. I have not carried it out
[outline]. If it holds, the floor is **1 + 3σ: 19 for E7, 25 for E8 and P112**. Then no design that saves at the hub's
own charge fits in block 16: only D6, whose hub-type charge 0.8 is unproved, would fit.

### 2.7 What a gauge would change

A free per-array address map moves an array's frame without paying. It could carry an absorbed window out of the way of
the next call's ports. The hub's gluing between calls (Φ_k(1) = Φ_{k+1}(t), our `HDef.lean` as read by an-build) may be
of this kind. Each array's frame dimension still only grows, so the dimension count of 2.4 survives a gauge. But the
orthogonality step V_s′ ⊥ V_s compares frames at different times, and a gauge breaks that comparison. So the wall is
proved here for the gauge-free class only. Proposal 3 was posed in that class (an-synth's own test: "monotone frames
except on centres and one common frame per gate").

## 3. The bridge charge for small designs

### 3.1 The index-code floor, in any network [proved]

For a set Q of ports whose labels are pairwise non-orthogonal, Q's data must leave the arrays that are high for every
port of Q, which are the arrays whose frame is not inside any u_q^⊥. Gates inside that set are invertible, so they
keep the rank of Q's data there. Only an array leaving the set (a descent) lowers that rank, and by at most one. So
the network has **at least |Q| descents**, for the largest such Q (ω). This is the brief's "t ≥ largest mutually
non-orthogonal port set". It holds for every gate order and every way of moving data, not only for the hub's calls.

### 3.2 A floor on the charge of one call, valid for every decomposition [proved]

Write the call's bridge as rank-one pieces: bridge k gathers Σ g_kb x_b and scatters r_ka times it. On any two ports
whose labels are not orthogonal, the delivered map must be exact, because no side information can pass between them
(2.1). Bridge k's frame before its descent contains the labels of the ports it gathered, and after its descent it lies
inside the orthogonal complement of the labels of the ports it scatters to. So **its drop is at least the F₂ rank of
the label pairing between its scatter set and its gather set.**

Now let a set Z of ports split into groups such that:
- labels inside a group are pairwise non-orthogonal;
- labels in different groups are orthogonal.

The delivered map is the identity on each group Q, so at least |Q| bridges touch Q on both sides. A bridge that touches
c groups has a pairing that contains a c × c identity, so it drops at least c. Summing over bridges, **the call drops
at least |Z| levels.**

`calc/classfloor.py` finds the best Z (exactly for at most 20 ports) [computed]:

| design | ports | best split set Z | floor per port per call | hub-type charge t(k−1)/n | σ |
|---|---:|---:|---:|---:|---:|
| D4 | 12 | 12 (three frames of 4) | **1.00: saves nothing** | 0.667 | 2 |
| A3 | 6 | 6 | **1.00: saves nothing** | 1.000 | 2 |
| D5 | 20 | 12 | 0.60 | 1.000 | 4 |
| D6 | 30 | ≥ 12 | ≥ 0.40 | 0.800 | 4 |
| A4 | 10 | 6 | 0.60 | 1.600 | 4 |
| A5 | 15 | ≥ 6 | ≥ 0.40 | 1.333 | 4 |

D4 is the instructive case. The hub's formula t(k − 1)/n would say 0.667, but D4's labels are three orthogonal classes
(one label per frame). Every bridge must touch all three, so every bridge drops 3 levels. The hub's formula holds for
the hub's P112 (64 per call, checked in Lean), but it is not a law for other designs.

### 3.3 The design table [computed, `calc/designs.py`, rebuilt from scratch]

| design | ports n | centres t | ω | window k | junk per port | σ | proved floor 1 + 2σ | floor 1 + 3σ (third window fresh) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| E8 | 120 | 8 | 8 | 9 | 56 | 8 | 17 | 25 |
| P112 | 112 | 8 | 8 | 9 | 52 | 8 | 17 | 25 |
| E7 | 63 | 7 | 7 | 7 | 32 | 6 | 13 | 19 |
| E6 | 36 | 6 | 4 | 7 | 20 | 6 | 13 | 19 |
| D4 / D5 / D6 / D7 / D8 | 12 / 20 / 30 / 42 / 56 | n | 4 / 4 / 6 / 6 / 8 | 3 / 5 / 5 / 7 / 7 | 8 to 24 | 2 / 4 / 4 / 6 / 6 | 5 to 13 | 7 to 19 |
| A3 / A4 / A5 / A6 / A7 | 6 / 10 / 15 / 21 / 28 | n | 2 / 2 / 3 / 3 / 4 | 3 / 5 / 5 / 7 / 7 | 4 to 12 | 2 / 4 / 4 / 6 / 6 | 5 to 13 | 7 to 19 |

The E8, P112 and E7 figures agree with an-synth and v-an-synth. Every design is tight (A² = (n/t)A), and σ = k − 1
for all of them. The Fano and Hamming structures the brief names are E7 and E8 here: E7's labels come from the Fano
plane and E8's from the extended Hamming code.

## 4. A second line per array is not free [argued]

With two-sided lines, X_a must also stay inside v_a^⊥ and Y_a must contain v_a.
- **Side information.** Take a junk pair (a, b) cleaned without descent, with the copy erased against X_b. The erasing
  gate is at a frame that contains Y_a's frame, which contains v_a. That frame is X_b's frame, which lies inside
  v_b^⊥. So **v_a ⊥ v_b for every junk pair of the first line** [proved, under condition A].
- **The bridge.** The bridge's own diagonal route (it gathers X_a and scatters to Y_a) makes it climb v_a. If it is
  then erased against X_a, it must drop v_a again, which costs at least the F₂ rank of the second labels over its
  scatter set.
- **The window.** The cheapest second labels that satisfy the first point are a mirrored copy of the window in fresh
  coordinates. That costs a fresh window per backward call, and so a block of about 2 + 3σ (E7: 20; E8: 26).

This is the general form of what an-build found on our word, where the Y line was forced to one shared coordinate s
and the dirt gates paid it back. The table's "second line free" column is therefore a reference, not a case anyone
can expect to reach [argued].

## 5. Junk cleaning has nowhere to hide at small blocks [argued]

At the floors above, the data arrays have no idle levels (13 = 1 + 2·6 for E7). The hub makes its helpers free in
rank by parking their leftover climbs in idle climbs (an-ours). Without idle levels, every helper level is paid.
- At the hub's own helper density (504 helpers for 112 × 52 junk entries per call, about 0.087 per entry), E7 has
  63 × 32 junk entries per call, which would need about 175 helpers.
- At even one paid level each, that is more than the 63 levels a call's pairs gain.
- The hub's direct delivery (sources meeting targets with no helper) avoids helpers only through ordering, which
  2.1 constrains.

So the generous ceilings of sec. 6 (junk cleaning free) are far above what such a network could reach.

## 6. Ceilings [computed, `calc/ceiling.out`]

The equation is Σ_children (r/m)^(1−a) = W, per pair (W = 2). The children are the data arrays' m − u levels in blocks
of σ, plus the bridge charge in blocks of σ. Junk cleaning is free, and the chain is long (two calls per pair).

| design | block | u = 1 | u = 1.5 | u = 2 |
|---|---|---:|---:|---:|
| E7 | 13 (proved floor) | 0.0336 | 0.0857 | 0.1400 |
| E7 | 19 (third window fresh) | 0.0154 | 0.0389 | 0.0631 |
| E7 two-sided, mirrored | 20 | | | 0.0573 |
| P112 | 17 (proved floor) | 0.0339 | 0.0745 | 0.1165 |
| P112 | 25 (third window fresh) | 0.0152 | 0.0332 | 0.0516 |
| E8 (120 ports) | 17 / 25 | 0.0369 / 0.0165 | 0.0777 / 0.0346 | 0.1197 / 0.0531 |
| D6, hub-type charge assumed | 13 | 0.0132 | 0.0470 | 0.0822 |
| D4..A5 at 1 + 3σ (third window fresh), proved centre floors | 7 to 14 | ≤ 0.0401 | ≤ 0.0653 | ≤ 0.0969 |
| D4..A5 at the proved floor 1 + 2σ, proved centre floors | 5 to 10 | ≤ 0.0851 | ≤ 0.1272 | ≤ 0.1903 |
| D4..A5 at the proved floor 1 + 2σ, best charge known (proved for D4 and A3; hub-type for D5, D6, A4, A5) | 5 to 10 | ≤ 0.0277 | ≤ 0.0792 | ≤ 0.1660 |

If the wall did not exist, the same generous model would give 0.2 or more at block 10 for P112 and E8 (`ceiling.out`,
"if stages could share windows"). That is why the wall, not the pricing, decides proposal 3.

The rows at the proved floor 1 + 2σ show the remaining gap. They use the weakest centre charges I can prove for the
small D and A designs (0.4 to 0.6 per port per call). Every code known for them pays 0.8 to 1.6, and D4's and A3's true
floor is 1. Closing either the third window or these centre floors brings every one-sided and u = 1.5 entry below
0.09.

A side remark, not a lead [argued, untested]: these one-sided ceilings at 1 + 3σ (about 0.015 for E7 and P112) are
about nine times the record. They assume that pairs reuse each other's coordinates and that junk cleaning is free. The
hub's block 1 + 8(2P + 1) is set by hosting, which needs exactly the idle levels this layout removes (an-ours, stage
knob). So the gap between 0.015 and 1.68e-3 is the price of junk cleaning, not a missed lever.

## 7. Raw-array checks [computed]

- `calc/frames.py`, at blocks 1 to 5:
  - Φ_0 = I and Φ_whole = C^⊗m;
  - Φ is multiplicative on orthogonal sums;
  - Φ_u² = R_u, and Φ_{u⊥}Φ_u⁻¹ = C^⊗m R_u;
  - a nested step costs its dimension difference in directional kernels;
  - a lone pair is swapped and fully transformed in exactly 2m paid moves.

  The first run failed at block 3 because the all-ones port has an alternating complement (no orthonormal basis, so
  no frame). The run now picks a port whose complement has one.
- `calc/net.py`, a frame-checked simulator on raw complex arrays (every gate asserts one common frame, every step
  asserts nesting):
  - **N2**, the one-pair three-shear swap through three scratch bridges: paid 2m + 4, exact outputs C^⊗m R_u x and
    −C^⊗m y, at blocks 2 to 5;
  - **N3**, two pairs with orthogonal ports and one shared bridge x₁ + x₂: the junk is removed by copies that leave
    each source while it is still on its own line, and is erased against the source afterwards. Y_a = y_a + x_a
    exactly, bridge and copies end at zero, and Y_a's frame holds its junk partner's label, at blocks 3 to 5.

No network that saves was found, so there was nothing further to verify on raw arrays or to price.

## 8. Predictions (STATUS.txt, written 10:58) against outcomes

| prediction | outcome |
|---|---|
| P1 (85 %): no network at block ≤ 16 with D > 0 in the stage-structured class | **partly held.** Proved: no saving for D4 and A3, and no network that shares windows between an array's source calls. Not proved: D5, D6, A4, A5 at blocks 9 to 10 (gap: their centre floor) and E7 at 13 to 16 (gap: the count that finishes the third window). No network was found, and the generous one-sided ceilings there are at most 0.085 |
| P2 (75 %): the wall provable as "later junk needs ports orthogonal to what was absorbed" | **held:** proved as 2.1, 2.2 and 2.4 (block ≥ 1 + 2σ), given condition A; a backward and the next forward call cannot share a window (proved); the general third window needs one count (2.6) |
| P3 (70 %): no port set in windows ≤ 6 below one level per port per call, hub accounting | **wrong as stated:** the hub-type formula gives 0.667 for D4 and 0.8 for D6. But the formula is not a law, and D4's true floor is 1.00 (proved). The bare-floor half held: D5, D6, A4, A5 have proved floors of only 0.4 to 0.6 |
| P4 (90 %): raw-array calibration passes | **held** (frames, lone pair, N2, N3) |
| P5 (90 %): nothing near 1e-1 at block ≤ 16 inside the class | **held for the best charges known:** one free direction ≤ 0.04, and ≤ 0.09 at u = 1.5. **Not proved:** with only the proved centre floors, the small D and A designs reach 0.085 one-sided and 0.127 at u = 1.5 under free junk cleaning. Outside the fixed gauge nothing is proved (the 60 % part held) |

## 9. What this leaves

1. **For the 1e-1 bar.** No reading of the crossing class reaches the bar at any block without at least two
   unproved favourable assumptions. The wall stops shared windows at small blocks, and the design and centre
   ceilings stop large blocks (an-synth, v-an-synth). The bar needs a different mechanism. That is anatomy item 1.
   The premise of its item 2 (shared windows) is now settled for the gauge-free class, apart from the named gaps.
2. **The cheapest remaining step.** Finish the count of 2.6. It is the same kind of count as `calc/early.py`, applied
   to the order in which each bridge serves its targets. With it, the floor 1 + 3σ is proved for swaps, and E7
   leaves block 16.
3. **To make the exclusion airtight** (cheap, in order):
   - carry out item 2;
   - prove condition A, or show that erasing copies elsewhere only moves the window cost onto other arrays;
   - prove a centre floor of about 0.8 per port per call for D5, D6, A4 and A5. The touch-count relaxation of
     sec. 3.2 cannot give it: its best symmetric weighting for A4 is again 0.6, because a centre on two disjoint
     edges touches five groups at F₂ rank 1. The off-diagonal zeros of the code must enter the argument;
   - redo 2.4 with gauges allowed, which is the hub's setting.
   A skeptic should check 2.1 to 2.5 and the class floor first.
4. **Not for the bar, but recorded.** The class floor of sec. 3.2 is a new lower bound on the centre charge that holds
   for every decoder. It does not move the hub's P112 floor of 64 (its split sets are smaller), but it could be run on
   any proposed port design in seconds before any helper search.

## Files

All in `checks/wht43/anatomy/an2-smallnet/`:
- `STATUS.txt` (predictions, then progress), `RUNLOG.txt`, this `NOTE.md`.
- `calc/shears.py` → `calc/shears.json` (no exchange with two calls per bank; up to 4 banks saved, 5 printed).
- `calc/frames.py` → `calc/frames.json` (raw-array frame calculus, lone pair).
- `calc/net.py` → `calc/net.json` (raw-array networks N2, N3).
- `calc/designs.py` → `calc/designs.json` (design table, σ, block floors).
- `calc/classfloor.py` → `calc/classfloor.json` (the centre floor of sec. 3.2).
- `calc/early.py` → `calc/early.json` (early and late side-information copies, sec. 2.3).
- `calc/ceiling.py` → `calc/ceiling.out`, `calc/ceiling.json` (ceilings; calibration against an-synth's 0.096).

The raw-array programs run with numpy from our own environment `bound-check/g-p112-decoder/venv`; the others need only
the standard library.

## Sources

- Briefs: `../BRIEF2.md` (section an2-smallnet), `../BRIEF.md`.
- `../NOTE.md` (the anatomy round), `../an-synth/NOTE.md` (sec. 3 walls; sec. 4.3 proposal 3 and its test; `calc/cap.py`
  table (d)), `../v-an-synth/NOTE.md` (the cap at 5.59e-3; Wall 3 not proved), `../an-build/NOTE.md` and
  `../v-an-build/NOTE.md` (the second line, the dirt gates, the forward calls), `../an-openai/NOTE.md` (sec. 2: the
  accounting identity, the bridge, carriers, restoration), `../an-ours/NOTE.md`, `stages.out`, and
  `../an-build/calc/calib.py` (the replica frame chain, block 1 + 8 × stages).
- `../../meta/MAP.md` sections 0 to 1.4 (necessary conditions; proved lower bounds).
- Our Lean, read as text: `Work/Reframe/SNet.lean` (a stage applies
  GF or GB on every bank pair, helper sets shared by all stages).
- The classification of line systems at 60° and 90° (Cameron, Goethals, Seidel and Shult, 1976) is used only as
  an-synth used it, cited from memory there; nothing here depends on it.
