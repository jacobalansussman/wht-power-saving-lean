fill/: the fill of this package (who hosts each of the 271,728 roles of #2979's circuit), how it was made, its checks

The files
  FILL.json.gz       the fill in #2979's format (sha256 7fc766c6c3e1561cee60b258138a918020b96d85ea782a2b01ad0ee035fd082e)
  PROFILE.json       its paid profile per gauge copy: 1: 1,630,639; 2: 470,196; 3: 154,512; 4: 71,928; 5: 11,988;
                     6: 5; 7: 3; 8: 5,328; 17: 8; 19: 1; 23: 25,905; 48: 16,576 (sha256 5df79411...)
  ORDER.json, DISPATCH.json   the call order and dispatch as the run that made the fill wrote them; equal as JSON to
                     #2979's own files (the Lean data and the community's checks used #2979's own bytes)

What changes and what does not
  #2979's words, its 74-copy clock, the call order, the dispatch and the providers are unchanged; only the fill is new.
  Every bin of children 8 and 10 is empty, so each of the 16,576 rank-48 blocks there is a literal head in my lowering
  (8,288 + 8,288; the first operation on its own array, acting on that array alone). Unlike #2979's own fill and the
  every-head-kept fill of package/leverA, no other bin is idle at rank 48: the idle remainders are 1 (2,269 bins),
  6 (5), 7 (3), 17 (8), 19 (1), 23 (25,905) and 48 (16,576, the heads). The fill was chosen for the price after the
  full cut, so its uniform figure without the lever (1.6897030e-3) is lower than #2979's (1.7290340e-3).

How it was made (made-by/, README.txt there)
  With #3815's allocation model (Erwan Secousses): its per-tick Hall LP with the full cut's objective, a MIP on the
  generated columns (made-by/SOL-MIP0.json is the class solution), then #3815's reverse-time placement program, made
  tolerant: a role with no planned slot left at its call goes to a free place that obeys the same window, same-call
  and capacity rules. 9 roles needed that (8 of width 6 at call 553, 1 of width 4 at call 0); 8 planned width-6 slots
  stayed idle (made-by/u1-2979sep.fallback.json). The same fill came out again, byte for byte in its content, from
  the class solution on a fresh copy of the community repository (checks/hub/hx-admit-NOTE.md, step 1).

Checks
  eval/    an independent evaluator of the fill (evalfill.json: valid, 16,576 literal heads, the recount agrees),
           the lever pricer's head count (heads_fill.json: 16,576 literal; heads_fill-control.json, where the bridge
           is moved ahead of the sweep: 0) and the price (price.json)
  hub/     the community's own checks of #2979's line run on this fill, unmodified (hx-admit-NOTE.md and its logs):
           #2979's PR2839 owner engine with its controls, verify.py, the complete replay's fresh natives and the
           admission's owner, lift, rebind, discovery and outer-ledger steps, and the release's own verify.py all
           pass. The fill is rejected only where the hub asks for a higher uniform figure than the one it improves
           on (its strict-gain steps and the near-optimality bound of check_cap.py), because its uniform figure is
           lower; the full cut is what makes it worth more.
  price/   the full-cut price re-derived by three independent routes (hx-price-NOTE.md: exact rational enclosures on
           the 1e-24 grid, a 60-digit root, and 70-digit interval arithmetic; all give 1.78870456430e-3), the full-cut
           profile (profile2979Fc.json) and the Lean rate generator's scan (rate-2979Fc-generator-scan.json: the
           largest grid numerator is 17887045; 17887046 is refused)
