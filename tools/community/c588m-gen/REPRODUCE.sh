#!/bin/bash
# Regenerates the certificate package community/c588m from its inputs (Python only, no Lean; about a minute) and
# compares every generated Lean file, MODULES.txt and the certificate with the published ones, byte for byte.
#
# usage: tools/community/c588m-gen/REPRODUCE.sh [work dir]      (default: tools/community/c588m-gen/out, git-ignored)
#
# Inputs, all in this repository:
#   * Ryan Shea's published candidate certificate with 112 source copies, and his bank schedule and price
#     (third-party/shea256-fourier-transform-below-nlogn/lean/complex588/artifacts/);
#   * the same certificate with the copies folded back, and its annealed bank schedule (inputs/; made by
#     tools/community/mutable/REPRODUCE.sh);
#   * Shea's generator and verification scripts (third-party/shea256-.../research/complex-transfer-588-596/) and the
#     record templates of the Fourier chain (third-party/shea256-.../lean/complex588/third-party/record-templates/);
#   * the generators and checkers of this repository (tools/gx, tools/gen, tools/Compare.lean) and OpenAI's
#     statement (third-party/openai-math/);
#   * the Lean engine for mutable sources as written (tools/community/mutable/XMut/).
# The steps are those that made the package on 13 October 2026; each writes a JSON or log file into <work>/logs.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd)
W=${1:-$HERE/out}; rm -rf "$W"; mkdir -p "$W"/{logs,data/shea-artifacts,toolkit/tools,toolkit/third-party,scripts}
S=$ROOT/third-party/shea256-fourier-transform-below-nlogn
# the layout the generators were written for
cp -R "$HERE" "$W/tools"; rm -rf "$W/tools/out" "$W/tools/REPRODUCE.sh"
cp "$S"/lean/complex588/artifacts/* "$W/data/shea-artifacts/"
cp "$HERE/inputs/c588-no-copies.gcert.json.gz" "$W/data/"
cp -R "$S/research/complex-transfer-588-596/experiments" "$S/research/complex-transfer-588-596/verification" "$W/toolkit/"
cp -R "$ROOT/tools/gx" "$ROOT/tools/gen" "$W/toolkit/tools/"; rm -rf "$W/toolkit/tools/gx/out" "$W/toolkit/tools/gen/out"
cp "$ROOT/tools/Compare.lean" "$W/toolkit/tools/"
cp -R "$S/lean/complex588/third-party/record-templates" "$W/toolkit/third-party/"
mkdir -p "$W/toolkit/third-party/openai-math/lean/ComparatorChallenges"
cp "$ROOT/third-party/openai-math/lean/ComparatorChallenges/UniformFourier.lean" "$W/toolkit/third-party/openai-math/lean/ComparatorChallenges/"
cp "$0" "$W/REPRODUCE.sh"; echo "(see the repository's VERIFY.md)" > "$W/BUILD.md"
cp "$ROOT/community/scripts/"*.sh "$W/scripts/"
cd "$W"
PY="python3 -I"
# 1. the input certificate: Shea's candidate with the 112 source copies folded back; must equal inputs/
$PY tools/prep_cert.py data/shea-artifacts/candidate.gcert.json.gz data/c588-no-copies.gcert.json.gz data/c588m.gcert.json.gz > logs/prep_cert.json
# 2. exact preflights of the scalar side (Shea's rule gives SPar 10,21,4,8,8), kernel-term mirror, mirror controls
$PY tools/scalar_audit_z.py data/c588m.gcert.json.gz --out logs/scalar_audit_z.json > /dev/null
$PY tools/zmirror.py data/c588m.gcert.json.gz --spar 10,21,4,8,8 --espar 10,8,4,8,8 --margin > logs/zmirror-c588m.json
$PY tools/mirror_controls.py data/shea-artifacts/candidate.gcert.json.gz data/c588m.gcert.json.gz > logs/mirror-controls.json
# 3. local Lean modules (gxgen through Shea's driver, mutable-source kinds; the label mirror runs inside)
mkdir -p out/local
$PY -B tools/generate_local_lean_z.py data/c588m.gcert.json.gz --toolkit toolkit/tools/gx --out out/local > logs/generate_local.log
# 4. parse the generated data back; kernel-term mirror on the parsed text; tries
$PY tools/roundtrip.py out/local/Work/GCert/Data/Gen/C588M/S0.lean data/c588m.gcert.json.gz out/c588m-parsed.json.gz > logs/roundtrip.json
# 5. bank tables (annealed = headline, cyclic = fallback, packed = optional), schedK mirror, exact rate
mkdir -p out/price
$PY tools/make_sched.py data/c588m.gcert.json.gz "$HERE/inputs/588atm3_mutable_sources_no_copies.schedule.json.gz" toolkit out/price > logs/make_sched.log
$PY tools/pack_leftover.py out/price/price-annealed.json out/price > logs/pack_leftover.log
for v in annealed cyclic packed; do
  $PY tools/trie_check.py out/local/Work/GCert/Data/Gen/C588M/S0.lean data/c588m.gcert.json.gz out/price/price-$v.schedule.json.gz toolkit/tools/gx > logs/trie_check-$v.json
  $PY tools/generate_global_lean_z.py out/price/price-$v.json --toolkit toolkit/tools/gen --record toolkit/third-party/record-templates --out out/global-$v > logs/generate_global-$v.log
  $PY tools/generate_checks_z.py --out out/global-$v --original-challenge toolkit/third-party/openai-math/lean/ComparatorChallenges/UniformFourier.lean --comparator toolkit/tools/Compare.lean
done
# 6. predicted values of the Lean negative controls
$PY tools/control_facts.py data/c588m.gcert.json.gz out/price/price-annealed.schedule.json.gz > logs/control_facts.json
# 7. the engine modules at Work/XMut (import prefix only), the summary, the package
$PY tools/relocate_xmut.py "$HERE/../mutable/XMut" out/xmut > logs/relocate_xmut.json
$PY tools/preflight_summary.py "$W" > logs/PREFLIGHT.json
$PY tools/package.py "$W" > logs/package.json
# 8. compare with the published package
P=$ROOT/community/c588m; G=$W/pkg/c588m; bad=0; n=0
cmp -s data/c588m.gcert.json.gz "$P/data/c588m.gcert.json.gz" || { echo "DIFF data/c588m.gcert.json.gz"; bad=1; }
cmp -s "$G/MODULES.txt" "$P/MODULES.txt" || { echo "DIFF MODULES.txt"; bad=1; }
for d in overlay variants; do
  # the statement files of the package carry a modification notice added after the build (Apache-2.0, 4(b));
  # package_check.py same compares them without it
  while read -r f; do n=$((n+1)); cmp -s "$G/$d/$f" "$P/$d/$f" || python3 -I "$ROOT/tools/community/package_check.py" same "$G/$d/$f" "$P/$d/$f" || { echo "DIFF $d/$f"; bad=1; }; done < <(cd "$G/$d" && find . -type f | sed 's|^\./||' | sort)
  m=$(cd "$P/$d" && find . -type f | wc -l); k=$(cd "$G/$d" && find . -type f | wc -l)
  [ "$m" = "$k" ] || { echo "COUNT $d: published $m, regenerated $k"; bad=1; }
done
[ $bad = 0 ] && echo "REGENERATION CHECK PASSED: $n Lean and comparator files, MODULES.txt and the certificate are byte-identical (the statement files once their modification notice is removed)" \
             || { echo "REGENERATION CHECK FAILED"; exit 1; }
