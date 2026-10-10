# tools/emit/: where the certificate of the fourth result comes from

`gcert_emit.py` (Python 3, standard library only; Chafik Boukhalfa with Anthropic Claude assistance, Apache-2.0,
contributed in the pull request that added the fourth result) lowers the data that the tool chain of
CrocSwap/integer-mult-bounds writes for a source-assisted complex word to the format gcert/1 of this repository:

- `complex_frame_flow.py --witness` of their pull request #184 (icekylinx): the flow DAG — one node per
  (phase, frame), the value channels on every edge, births, retirements and reads — as `flow.json` and
  `flow.witness.json`;
- `exact_complex_flow_lift.py` of the same pull request: at every node an invertible completion as a list of
  elementary row operations, and the coordinates of every read, as `lift.certificate.json.gz`;
- `graph.json` and `selection.json` of the aligned word (their `source_aligned_local_v4.py`, pull request #194,
  ikeboy).

The emitter numbers one slot per physical flow coordinate and, node by node in the lift's topological order, at
the node's frame, emits the source copies, the reads (side and deferred responses ±1/2, as their `contract_v4.py`
checks; a retained total is formed in place in a retiring input slot), the completion as single-target
fan-outs (the lift stores the operations taking the completion to the identity, replayed in reverse; a swap only
relabels positions; a scaling changes the unit of a slot, as `tools/gx/gxconv.py` does), and names every
coordinate that passes through a node untouched as an extra register, so that the block histogram is the flow's
ledger. After the last node of phase 2 it applies the K blocks 1 - J/2 of their decoder (pull request #144) on
each parity class of a cube at its parity 3-space, reads each mixed source into its antipodal target at that
target's cap frame, and undoes the K blocks at the full frame.

Nothing of the outside repository is included, imported or executed here: the emitter reads the files above as
data, and `tools/certificate/gcert1-p11-pr233-flow.json.gz` is its output. Its inputs are regenerated and the
output compared byte for byte by `research/gcert-program-233/verify.py` of CrocSwap/integer-mult-bounds (pull
request #256), which also runs `tools/gx/gx.py` of this repository (vendored there unchanged) on the result.

Control: on the layer of their pull request #194 (their frames of #168, pairs of #193) the same emitter produces a
certificate whose block histogram equals, class by class and rank by rank, the one of
`tools/certificate/gcert1-p11-pr193.json.gz`, which was rebuilt independently (`tools/gx/ORIGIN.md`), with the
same N = 262,944; `gx.check1` accepts it and `gxdry.py` prices it at the 7474547/10^10 of the second result.

Usage (paths of their tree):

    python3 -B tools/emit/gcert_emit.py --cache <aligned>/cache --flow flow.json --witness flow.witness.json \
        --lift lift.certificate.json.gz --out gcert1-p11-pr233-flow.json.gz --name 'pr233 flow word'
