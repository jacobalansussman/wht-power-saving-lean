# Data of the section "why the record savers save"

- `inputs/` holds the small shared inputs that the programs read:
  - `U8X-parsed.json.gz`: my agents' lowering of the community's U8 word (Jeremy Kothe's, with Chafik Boukhalfa's
    frame changes and icekylinx's gauge44 construction; the pull requests are in the Sources of
    [../README.md](../README.md)) with exit frames, in this repository's certificate format (`gcert/1`: 112 ports,
    window 9). It is data derived from the community's files
    (CrocSwap/integer-mult-bounds, Apache-2.0), the same word whose Lean data is in `community/u8x/`.
  - `words.json`: the profiles of my words, the same file as `community/u8xmhh/data/words.json`.
  - `hub-checkpoint-profiles.jsonl`: the paid-move profiles and figures of six of the community's checkpoints
    (from Kothe's atm3 to his k8 frontier; each row names its checkpoint), as numbers read from its published files;
    `lineage.py` adds the T74 line as a seventh row.
  - `identity.json`: the per-call identity of my words (paid rank, narrow mass, residual helper widths and D = 320 per replica for each word).
  - `U8X0-checker-histograms.json`: my checker's per-class histograms of the U8X0 word, re-run on 15 October.
- The other files are copies of key results and certificates that also sit beside their programs in `../programs/`:
  the centre theorem's results (`centre-floor-64-results.json`), the copy-floor calibration, the k8 census and
  rank-7-pair results with their skeptics' verdicts, the invariant-decoder and word-rewrite summaries, the
  stage-equivariant summary, the lattice-unit ceilings and what-ifs, the anatomy of my words and the two centre
  design tables.
- `hub-inputs/README.md` lists the community's files that some programs read and that are not copied here.
