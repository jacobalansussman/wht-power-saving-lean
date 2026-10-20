/-
Local packaging by fs-build (13 Oct 2026) for an offline re-creation of the Big O Speedrun verifier check;
not a submission. Same packaging as nmamano's submissions/dft/nmamano-shea-c588 (m-build's script): the proof
closure of Work.FourierU8XMHH.Main is copied under Submission/, changing local import paths only.
Same-run hosting: DaysSky's #442 hosting argument (x-hosting; h-L3's hosted network Work/XHost); mutable sources
with 5 centre-source copies kept (x-mutable); exit frames (x-exit-lean, k-lower, k-combine).
Credits: the hub's U8 word (Jeremy Kothe #869/#870, Chafik Boukhalfa #890, Capalbo #942, on Kothe's #799/#588
lineage), FanPUYen/Hao Yan (owners), the hub authors (exit frames and tail #591/#593/#599); Ryan Shea's Lean
instantiation (generators, bank and seam pipeline); Nil Mamano's gauged Lean extension and record files;
Dugongue (#401) and jacklightChen (#377), Zhihao Chen (#425); OpenAI's model and Fourier reduction; Jacob Sussman's
engine, checkers and generators, with the exit-frame (XExit, KLow, KMX) and mutable-source (XMut) engine sides.
-/
import Submission.Work.FourierU8XMHH.Main

open OAI.PowerSaving Filter Asymptotics

theorem Race.bound : ∃ order solve W, DFTProgram order solve W ∧
    (fun n : ℕ => (W n : ℝ)) =O[atTop]
      (fun n : ℕ => (n : ℝ) * (Real.log (n : ℝ)) ^
        (1 - (3239541 : ℝ) / 2000000000)) := by
  obtain ⟨order, solve, W, hprogram, htime⟩ := transform_mainU8XMHH
  refine ⟨order, solve, W, hprogram, ?_⟩
  have hdecimal : decimalTimeU8XMHH =
      (fun n : ℕ => (n : ℝ) * (Real.log (n : ℝ)) ^
        (1 - (3239541 : ℝ) / 2000000000)) := by
    funext n
    have hexponent : decimalExponentU8XMHH = 1 - (3239541 : ℝ) / 2000000000 := by
      norm_num [decimalExponentU8XMHH]
    simp only [decimalTimeU8XMHH, hexponent]
  exact hdecimal ▸ htime.2.1
