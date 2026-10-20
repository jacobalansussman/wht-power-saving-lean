#!/bin/zsh
# fs-price: same-run hosting runs (controls + the U8XM stack). Run through heavy.sh.
set -e
F=${F:-$PWD}  # as run: the agents' fs-price folder
cd $F
for spec in "588atm3_shea_route_control h588_c0_samerun" "869/870/890_U8_exit_frames(E1) h869_c0_samerun" "U8XM_5cc hu8xm_c0_samerun"; do
  name=${spec% *}; tag=${spec##* }
  echo "== $name -> $tag $(date +%T)"
  python3 -I tools/host_price.py data/words.json "$name" --copies 0 --defer 0 --samerun 1 --tl 600 --out out/$tag.json
  python3 -I tools/host_repack2.py data/words.json out/$tag.json out/$tag.repacked2.json
  python3 -I tools/audit_slots.py data/words.json "$name" out/$tag.json
  python3 -I tools/audit_slots.py data/words.json "$name" out/$tag.repacked2.json
done
echo "== done $(date +%T)"
