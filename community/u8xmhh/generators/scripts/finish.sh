#!/bin/zsh
# fs-price: from a same-run hosting result of U8XM (host_price --samerun 1 + host_repack2) to the generator inputs and the
# Lean texts of GXD.U8XMHH (Python only, seconds).  usage: finish.sh TAG [BANKS.json]  (TAG = hu8xm_sr; BANKS default out/TAG.repacked2.json)
set -e
F=${F:-$PWD}  # as run: the agents' fs-price folder
BC=${F:h}
cd $F
tag=$1; banks=${2:-out/$tag.repacked2.json}
# 1. host table of one run (slot order of the emitter; same table in all 8 runs is asserted)
python3 -I tools/emit_hosttab.py data/words.json out/$tag.json U8XM_5cc out/HU8XM-hosttab.lean out/HU8XM-hosttab.json
# 2. data modules, rate certificate, unit price, glue (mirror of hostSchedK, histogram and plan_block at A / A+1 inside)
rm -rf lean && mkdir -p lean
python3 -I tools/gen_hu8xm.py lean out/HU8XM-hosttab.json $banks
# 3. Fourier files with the seam on XHost.engine_of_ginv_host (k-combine's Lean-checked U8XM Fourier files renamed)
python3 -I tools/mk_seam_spec.py U8XM U8XMHH out/U8XMHH.json 13868857 10939 out/seam_U8XMHH.json
python3 -I tools/mk_fourier_hx.py $BC/k-combine/pkg/k-combine-u8xm/overlay lean U8XM U8XMHH out/seam_U8XMHH.json
grep -rn "13868857\|13868856\|10939" lean/Work/FourierU8XMHH lean/Work/FourierU8XMHHNeg && { echo "old numbers left"; exit 1; } || true
echo "finish ok $tag"
