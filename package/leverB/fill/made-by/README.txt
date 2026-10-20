fill/made-by/: how the every-head-kept fill of #3945 was made

The idea. In the community's #3815 / #3945 construction, every owner copy has rank-48 bins in two child pairs: children
8 and 10 and children 12 and 13. The construction's LP prices an empty rank-48 bin the same in either pair, so its
solver splits the empty bins between them arbitrarily; as published, 57,312 empty bins sit in children 8/10. An empty
child-8/10 bin is a literal head in my lowering (the first operation on its own array, acting on that array alone),
and only heads can take the half cut. So I gave the LP a bonus for exactly that pattern, equal to what the half cut
saves on it, and let the community's own programs do the rest. The paid profile, and so the community's own figure,
does not change; only which bins carry the hosted roles.

The steps (on a Colab, CPython 3.12 with #3815's recorded numpy 2.5.4, scipy 1.18.1, highspy 1.15.1; hub files read
from fresh worktrees of PR heads 77ffbf28 (#3815) and 91f81d13 (#3945), never modified):

  1. probe/mkpatch3815.py hall_lp.py(of #3815) research/clock-heads/construction/hall_lp.py 1
       a copy of #3815's hall_lp.py with one change: the column cost of "class 9r+4, empty rank-48 pattern" is lowered
       by the half cut's value delta(a)/a = (48^t - 24^t - 73^t + 49^t)/a, t = 1 - a (12.853 at #3815's root); the
       same bonus enters the reduced costs. research/clock-heads/construction/hall_lp.py is that copy.
  2. probe/hall_mip.py <construction dir with the patched hall_lp.py and #3815's build_fill.py> <#3815 ORDER.json> out 900
       #3815's build_fill.build() with the patched LP: column generation and Hall-row separation to convergence (LP
       root 0.0017305784372663 = #3815's, 80,640 heads, 11 fractional classes), then the same restricted model as a
       MIP (HiGHS: Optimal in 195 s, the LP's objective), re-separation on the MIP point (0 new Hall rows), then
       #3815's placement unchanged (--distinct 40). Output: the #3815-stage fill (../old-FILL.json.gz, 0e61fce2...),
       with #3815's idle profile {1: 84,240, 3: 23,328, 5: 4,320, 48: 137,952}. PR3734's owner_audit.py: PASS.
  3. probe/offsets3945.py <#3945 owner dir> out/inputs/FILL.json.gz final3945
       #3945's offset step re-implemented (keep every (t, h, o, c, p); recompute compact offsets with #3945's final
       widths, 4 -> 5 for its 26 widened helpers; no overflow): ../FILL.json.gz (9c06677d...) and ../PROFILE.json
       (the content of #3945's PROFILE.json). #3945's own owner/build_fill.py does the same step; as a control it
       regenerates #3945's published fill from #3945's own #3815-stage fill (PASS_REGENERATED_LITERAL_FILL).
  probe/run_probe3945b.sh is the script that ran steps 1-3 (o-model's COMMANDS.txt section 5, probe P3b);
  p3b-mip1.out is its log.

The checks of the fill (all on Colabs, hub code unchanged, fresh copies; ../checks/ and the round's gate):
  - #3945's owner/verify.py, plain and -O, on a copy of owner/ with only old-FILL, FILL, PROFILE and their pins
    replaced: PASS_LITERAL_ENLARGED_OWNER_FILL, its four rejection controls rejected; the published files pass the same
    harness (verify3945-allheads.out).
  - #3945's complete admission replay (fresh natives on both providers, owners, parent lift, rebind, discovery, both
    outer ledgers, the release with all 21 verify sub-stages): PASS_TERMINAL_WIDTH5_COMPLETE_CONDITIONAL_RELEASE,
    119 of 119 strict controls rejected, kappa exactly #3945's claim.
  - Heads counted five ways (two independent counters, lever-reprice's recount.py and heads_fill.py, o-eval, my
    lowering's G0): 80,640 = 40,320 + 40,320 each, 0 in children 12/13; the paid profile equals #3945's PROFILE.json
    digit for digit; the half cut priced by two routes to 24 digits: 0.001758152016764631555031.
  - An exact LP-dual certificate (90-digit arithmetic, all 841 patterns): no fill of #3945's words, with any
    chronology, prices the half cut above 1.7581520e-3 by more than 4e-34. So this fill is the best re-allocation of
    #3945's hosted roles for the lever.
