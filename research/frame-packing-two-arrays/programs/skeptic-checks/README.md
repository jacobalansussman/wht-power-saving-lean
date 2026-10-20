# skeptic-checks: the skeptics' own programs

Each skeptic wrote its own engine, frame game or solver from the definitions, without importing the track's code
(except where a README line says the track's code was used only to export test programs as data). These are the
programs behind "a skeptic repeated it with its own code" in the section README. Outputs are in
`../../data/skeptic-checks/<same folder>/`.

| folder | checks | main commands |
|---|---|---|
| [slack/](slack/) | frame facts from actual matrices at blocks 1 to 3; the order-3 control as an exact identity; the block-2 stand-in 2 in 3 is right in one representation and is no engine unit; every ceiling and parity figure; the spectral cut ceiling | `python3 vframes.py`; `python3 vctl.py`; `python3 vstandin.py`; `python3 vunits.py`; `python3 vspectral.py` |
| [heights/](heights/) | an independent true-model engine; every named unit rebuilt; Theorem C's chain with mixed height systems; Theorem A's steps; the false reading of the common-frame height statement | `python3 vunits.py`; `python3 vchain.py`; `python3 vthmA.py`; `python3 vcf.py`; `python3 vraw.py 3:6,7 4:4,5` (raw arrays, large f) |
| [packing-measurement/](packing-measurement/) | its own replay of the forward and backward calls and the positions (they read a word in gcert/1 format); its own engine for the calibration units; the level formula; OpenAI's network derived from the preprint | `python3 vcal.py`; `python3 vlevel.py`; `python3 voa.py`; `python3 voa_p.py`; `python3 vfm.py <word> F1` |
| [packing-capacity/](packing-capacity/) | its own engine; Theorem A, A′ and P* on rebuilt units; Hu's LP; the geometry; the counterexample to "image frames are the sharpest count" | `python3 vsharp.py`; `python3 vunits.py SEED N_RANDOM` (outputs for seeds 1, 2 and 4); `python3 vgeom.py 3 0`; `python3 vhu.py SEED N_RANDOM` (output for seed 3); `python3 vtables.py` |
| [frame-game/](frame-game/) | the faithful block at block 2 from actual 4 × 4 matrices; its own frame-game arithmetic; landings; the recount of the design tallies; the refutation of the "closest design"; how many repairs only rescale; its own linearised position test | `python3 vgroup.py`; `python3 vland.py`; `python3 vtally.py`; `python3 vhit.py`; `python3 vselfinj.py 4 30 12 61`; `python3 vpsilin.py 4 30 12 61` |
| [frame-game-completion/](frame-game-completion/) | its own engine and lift; exact identities at blocks 2 and 3; all 10,800 products at block 2 | `python3 vexh2.py`; `python3 export_progs.py 2` then `python3 vcheck.py 2 8`; `python3 myprogs.py 2 11` then `python3 vcheck.py 2 8 full,f1 progs_mine_m2.json`; `python3 vhom.py m n_words length` (see its header) |
| [closed-patterns/](closed-patterns/) | its own frame game and closure solver; the parity law on 40 cases | `python3 v_lemma1.py` |
| [single-term-game/](single-term-game/) | its own groups, game and solver; the game against explicit modules; the block-4 cube 2 in 7 re-run; Conjecture S; the sub-game lemma with symplectic matrices; no involutive lift of an apartment swap | `python3 vb4.py B4 7 4 0 6` (part 0 of 6); `python3 vinv.py D7 6`; `python3 vmod.py D7 PLAYS STEPS SEED`; `python3 vlem43.py`; `python3 vlift.py`; `python3 vcross.py B3 DEPTH`; `python3 vbeam.py H43 7 300 1.0 2` |

Some argument lists above are examples. Each file's header gives its usage. In `frame-game/`, `vtally.py` recounts
every design tally from the build track's raw outputs in `../../data/frame-game/`, and `vhit.py` re-decides the
"closest design" from the build track's pickle `../../data/frame-game/design-runs/hit_m4_s102_4.pkl`, which
`loadhit.py` reads with a restricted unpickler (numpy arrays only). `vselfinj.py` and `vpsilin.py` generate designs
with the build track's `../frame-game/design.py` and judge them with the skeptic's own arithmetic. On 19 October all
four reproduced their shipped outputs (`vselfinj.py` and `vpsilin.py` at seed 61); the shipped `out_vhit.txt` repeats
one line, an artefact of the original run, which the re-run prints once.
