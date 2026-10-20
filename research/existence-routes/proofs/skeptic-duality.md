<!-- header: existence-routes -->
# The skeptic's record: duality

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** This is an AI skeptic's record of its check, as it wrote it. Its verdict is the status this revision uses for the track it checked, unless a later verdict is named in the section README.

---

# v-ex2-duality: skeptic's verdict on ex2-duality

14 October 2026, two sessions (20:35 to 21:03, 23:01 to 23:14 EDT). Nothing posted, no Lean.

**Verdict: not refuted, and the calibration holds: the route is blind both ways, as the track says.**

## What I rebuilt with my own code
- The identity behind "one direction costs at most two kernel calls" (A5), exactly in Q(i), blocks 1 to 6.
- Certificate A's prices, exactly at widths 3, 4 and 5: a rank-r block moves 2^w - 2^(w-r) of the frame's labels.
- An independent fluid LP. It reproduces every number in the amplification table, in both the "=" and ">=" forms and
  with caps from n to 6n. At n = 243, the hypothetical pair R1 + Rs gives the fixed relation 53.16 x [1] >= [243]. Its
  powers have exponent 0.723, below the 0.834 of R1 alone, so amplification shows up as one finite relation.
- Every Theta value the track quotes (0.8811, 0.8214, 0.8844, 0.7813, 0.834; R2's two pieces; 0.1926 for R1 + R2).
  The grouped test rule and its rank-1 form have identical balances, which is why their LPs agree.
- A literal run of the pipelining schedule. It is deadlock-free; the boundary rounds hold L(L-1) calls.

## What I checked by reading
- The Fritz citations are right. Theorem 7.15 counts only tropical points with phi(u) > 1, and there are none when
  u is a number. Theorem 8.1 is the same statement for this Archimedean case, and Lemma 2 supplies its hypotheses.
- The one step the track left [argued] (phase letters after regrouping) is not a gap in the job semiring. A regrouped
  engine phase is again an engine phase.
- Certificate A works at every block for the engine's letters, and Theorem R allows targets of any rank. So the
  spreading family is excluded at every block.
- The Closure Lemma is correct, though close to a tautology. The amplification corollary is correct mathematics. Its
  correction of ex-note 2.3 is fair: that note says "no mixture beats the best single unit" with no qualification,
  while ex-laser's Theorem 3 needs one common exponent. Every amplifier tested so far is excluded by certificate A.

## Small slips (none changes a conclusion)
- The pipelining overhead is of order L^2/N, not L/N.
- Corollary 6(c) should cite the concavity of r^theta, not the convexity of D.
- The PROOF.md header says the track ran to about 21:30; it closed at 20:38.
- out_calib_A.txt shows 0.8843 where the root is 0.884386.

Not checked: that the job model's exponent equals the challenge's (ex-laser 1.3, argued there), and Fritz part II for S
(the track leaves it open). Files: STATUS.txt (predictions and log), RUNLOG.txt, the scripts and their out_*.txt.
