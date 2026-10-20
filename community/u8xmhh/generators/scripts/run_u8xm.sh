#!/bin/zsh
# fs-price: same-run hosting of U8XM (509 residual widths: U8 + exit frames + mutable sources + 5 centre copies)
set -e
F=${F:-$PWD}  # as run: the agents' fs-price folder
cd $F
tag=${1:-hu8xm_sr}
echo "== U8XM_5cc -> $tag $(date +%T)"
python3 -I tools/host_price.py data/words.json U8XM_5cc --copies 0 --defer 0 --samerun 1 --tl 600 --out out/$tag.json
python3 -I tools/audit_slots.py data/words.json U8XM_5cc out/$tag.json
echo "== repack2 $(date +%T)"
python3 -I tools/host_repack2.py data/words.json out/$tag.json out/$tag.repacked2.json
python3 -I tools/audit_slots.py data/words.json U8XM_5cc out/$tag.repacked2.json
echo "== done $(date +%T)"
