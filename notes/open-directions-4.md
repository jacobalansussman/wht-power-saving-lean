# Open directions, part 4: the recursion lever (19 October 2026)

Each lead has its price, its status and a first experiment. A price is computed by the balance of the README
("The lever, with the numbers"); none of these is proved in Lean.

## 1. Rules that are not separated (the largest lead)

**Price on #3945:** 1.8033130e-3 (+4.20 %), cutting each head down to a job of 25 axes and keeping the other
80,640 blocks of rank 48 as they are.

**Status:** Theorem 1G needs every paid block to be smaller than every job. Here the blocks of rank 48 are larger
than the job of 25, and no theorem covers that case. It is not known to be false either: the agents' proof uses
separation to bound the jobs that are left over, and a different construction might not need it.

**First experiment:** simulate the lockstep tree of calls for small non-separated rules in the job model (the
simulations of the proof of Theorem 1G run in seconds) and look for one where the left-over jobs grow faster than
the theorem allows. A counterexample would close the lead; none would suggest a proof.

## 2. Heads with hosted helpers in front of them

**Price on #3945:** 1.7600092e-3 (+1.70 %), the half cut plus 27,648 partly hosted blocks cut fully (23,328 of
rank 3 in the same two children as the heads, and 4,320 of rank 5 in two other children).

**Status:** argued only, read by no skeptic. These blocks are the first operation on their arrays except for hosted
helper storage placed there earlier. The head lemma needs nothing before the head; a version that allows helper
roles that do not touch the head's axes would cover them. The Lean proof already has a commuting route (a block on
one array moves to the front past every step that leaves that array alone, `Work/Mixed/HeadComm`).

**First experiment:** for each such block, list the hosted roles in front of it and check whether their cells are
disjoint from the block's axes (a few seconds of Python on `results/g0/g0-3945.json` and the owner fill).

## 3. The other 80,640 blocks of rank 48

**Status:** they are neither heads nor tails: the community's controller works on their arrays before them and reads
them after. Turning them into tails would need the bridges that read them to be reordered, which changes what
the next bridge reads: a new construction, not a property of the current one.

**First experiment:** search over orders of the bridges involved (on T74 they were three), with the community's
exact operator replay as the test of correctness.

## 4. The next circuit of the community

The lever is not tied to #3945. Every newer circuit of the community that keeps T74's schedule has heads in the
same places, and the half cut gave about +1.1 % on each of them (the README's table of lever prices). The draft #4076, built on
#3945's words with three more staged helper entries, would give 1.7504268e-3.

**First experiment:** when a new circuit appears, run `tools/lever/heads_fill.py` on its owner data and
`tools/lever/lever_price.py` with its profile; the whole Lean route then needs only new data and one build.

## 5. Integer multiplication

The lever raises the complex side above the community's selected bit side, 1.7334270e-3 (#589), and above the higher
bit sides claimed there since (for example #3647, ai-nanachi, 1.7366460e-3, archived on 20 October and not selected).
The community's headline
figure, kappa, is capped both by the complex side and by a bit side, so it would rise only as far as its bit side. Raising it further needs work on the bit side,
which this repository does not touch.
