#!/bin/bash
# build-and-check.sh <package dir> <tree = .../lean/complex588, overlay already applied>
#
# (Adapted from the script of the same name that I wrote for my earlier packages: it also re-derives the two
# statement files from OpenAI's file in third-party/openai-math and compares them byte for byte, and it runs the
# Lean controls that result.json lists under "lean_controls"; each must compile.)
# (Changed on 19 October 2026, when it was first run on the lever packages from a fresh clone of this repository: a
# Lean control listed with "expect": "fails" must fail to compile, with a Lean error and not a missing file; the
# lakefile lines added by the packages named under "requires" are allowed as well; and the axioms of every name in
# "axiom_targets" are printed from a file written into the log folder with the imports of the axioms file, because
# the axioms file of a lever package prints only some of those names.)
#
# Builds the new modules of a certificate package and runs the Lean checks named in its result.json:
#   1. build: one `lake build <module>` per line of MODULES.txt (or the build_targets of result.json), in order,
#      waiting while 2 or more Lean processes run on the machine (a Colab-sized machine has room for two);
#   2. axioms: `lake env lean <axioms_file>`; every target must depend on propext, Classical.choice, Quot.sound only;
#   3. compare: tools/Compare.lean, challenge against solution: must print RESULT: PASS;
#   4. compare-negative: the challenge one unit stronger: must FAIL, on the exponent definition;
#   5. source scan of the overlay files (sorry, admit, axiom, native_decide, implemented_by, extern, unsafe, opaque),
#      and no tracked file of the tree changed, except lakefile.lean by the line that the package's result.json
#      names under "tree_change" (u8x, u8xm and u8xmhh need the XExit library line).
# Each step runs under /usr/bin/time -v; one log per step in $LOGDIR (default ./logs/<package name>).
# Prints one PASS or FAIL line per check and ends with ALL CHECKS PASSED (exit 0) or SOME CHECK FAILED (exit 1).
# The official comparator is a separate step: run-comparator.sh.
set -u
PKG=$(cd "${1:?package dir}" && pwd); T=$(cd "${2:?tree}" && pwd)
RREPO=$(cd "$(dirname "$0")/.." && pwd)   # this repository (resolved before the script changes directory)
NAME=$(basename "$PKG"); LOGDIR=${LOGDIR:-$PWD/logs/$NAME}; mkdir -p "$LOGDIR"
export LEAN_NUM_THREADS=${LEAN_NUM_THREADS:-2}
J() { python3 -I -c 'import json,sys; d=json.load(open(sys.argv[1])); v=d
for k in sys.argv[2].split("."): v=v[k]
print(" ".join(v) if isinstance(v,list) else v)' "$PKG/result.json" "$1"; }
TIME=/usr/bin/time; [ -x $TIME ] && $TIME -v true 2>/dev/null && TV="-v" || TV="-l"
step() { local n=$1; shift; local s=$(date +%s)
  echo "=== STEP $n START $(date -u) cmd: $*"
  $TIME $TV -o "$LOGDIR/$n.time" "$@" > "$LOGDIR/$n.log" 2>&1; local rc=$?
  echo "=== STEP $n END rc=$rc wall=$(( $(date +%s) - s ))s last: $(tail -1 "$LOGDIR/$n.log")"
  return $rc; }
# wait while MAX_LEAN (default 2) Lean processes run on the machine, or while less than MIN_FREE_MB (default 0: no
# check) of memory is available; a Compare run needs about 16 GB
avail() { local a; a=$(free -m 2>/dev/null | awk '/Mem:/ {print $7}'); echo "${a:-999999}"; }
nlean() { local n; n=$(pgrep -c -x lean 2>/dev/null); echo "${n:-0}"; }   # pgrep -c prints 0 and exits 1 when none run
waitlean() { while [ "$(nlean)" -ge "${MAX_LEAN:-2}" ] || [ "$(avail)" -lt "${MIN_FREE_MB:-0}" ]; do sleep 15; done; }
cd "$T" || exit 90
echo "START $(date -u) tree=$T HEAD=$(git rev-parse HEAD 2>/dev/null) package=$NAME"
if [ -f "$PKG/MODULES.txt" ]; then MODS=$(grep -v '^#' "$PKG/MODULES.txt" | awk 'NF {print $1}'); else MODS=$(J build_targets); fi
i=0
for m in $MODS; do
  i=$((i+1)); waitlean
  step "$(printf '%02d' $i)-$m" lake build "$m" || { echo "BUILD FAILED at $m"; echo "SOME CHECK FAILED"; exit 1; }
done
waitlean; step build-targets lake build $(J build_targets) || { echo "SOME CHECK FAILED"; exit 1; }
waitlean; step axioms lake env lean "$(J axioms_file)"
# every name of axiom_targets, from a file with the imports of the axioms file
{ grep -E '^import ' "$(J axioms_file)"; for t in $(J axiom_targets); do echo "#print axioms $t"; done; } > "$LOGDIR/AxiomTargets.lean"
waitlean; step axiom-targets lake env lean "$LOGDIR/AxiomTargets.lean"
# the statement files, re-derived from OpenAI's statement (only the exponent and seven names change)
REF=$RREPO/third-party/openai-math/lean/ComparatorChallenges/UniformFourier.lean
SUF=$(J suffix); CHF=$(J challenge_module | tr . /).lean; NGF=$(J negative_challenge_module | tr . /).lean
bash "$RREPO/scripts/derive_challenges.sh" "$REF" "$(J numerator)" "$SUF" > "$LOGDIR/derived-pos.lean"
bash "$RREPO/scripts/derive_challenges.sh" "$REF" "$(J negative_numerator)" "$SUF" > "$LOGDIR/derived-neg.lean"
{ cmp "$LOGDIR/derived-pos.lean" "$CHF" && echo "challenge: identical to OpenAI's statement with $(J numerator)"; } > "$LOGDIR/derive.log" 2>&1
{ cmp "$LOGDIR/derived-neg.lean" "$NGF" && echo "negative challenge: identical to OpenAI's statement with $(J negative_numerator)"; } >> "$LOGDIR/derive.log" 2>&1
# Lean controls (each file must compile: for example the proof that the rate fact fails one unit higher, or that
# the host-table checker rejects a table with one entry off by one)
python3 -I -c 'import json,sys; [print(c["name"], c["file"]) for c in json.load(open(sys.argv[1])).get("lean_controls", [])]' "$PKG/result.json" > "$LOGDIR/controls.txt"
rm -f "$LOGDIR/controls.rc"
while read -r cn cf; do waitlean; step "control-$cn" lake env lean "$PKG/$cf" < /dev/null; echo "control-$cn rc=$?" >> "$LOGDIR/controls.rc"; done < "$LOGDIR/controls.txt"
CH=$(J challenge_module); NG=$(J negative_challenge_module); SO=$(J solution_module); TH=$(J theorems)
waitlean; step compare lake env lean --run tools/Compare.lean "$CH" "$SO" $TH
waitlean; step compare-negative lake env lean --run tools/Compare.lean "$NG" "$SO" $TH; NEG=$?
{ echo "# overlay files: word matches of sorry admit axiom native_decide implemented_by extern unsafe opaque sorryAx"
  (cd "$PKG/overlay" && find . -name '*.lean' | sed 's|^\./||' | sort) | while read -r f; do
    grep -nwE 'sorry|admit|axiom|native_decide|implemented_by|extern|unsafe|opaque|sorryAx' "$T/$f" | sed "s|^|$f:|"; done
  echo "# tracked files of the tree changed:"; git status --short --untracked-files=no
  echo "# lakefile.lean diff:"; git diff -U0 -- lakefile.lean | grep -E '^[-+]' | grep -vE '^(\+\+\+|---) '
} > "$LOGDIR/source-scan.txt" 2>&1
python3 -I - "$LOGDIR" "$NEG" "$PKG/result.json" "$PKG/overlay" <<'PY'
import json, re, sys
from pathlib import Path
L = Path(sys.argv[1]); neg = int(sys.argv[2]); R = json.load(open(sys.argv[3])); ov = Path(sys.argv[4])
ok = []
ok.append(("compare: RESULT: PASS", "RESULT: PASS" in (L/"compare.log").read_text()))
bad = (L/"compare-negative.log").read_text()
ok.append(("compare-negative: fails (rc != 0)", neg != 0))
ed = R.get("exponent_definition", "").split(".")[-1]
ok.append(("compare-negative: the exponent definition differs (%s)" % ed, "differs" in bad and ed in bad))
ax = (L/"axioms.log").read_text() + "\n" + ((L/"axiom-targets.log").read_text() if (L/"axiom-targets.log").exists() else "")
ok.append(("axioms: no sorryAx", "sorryAx" not in ax))
dl = (L/"derive.log").read_text()
ok.append(("statement files: both identical to OpenAI's statement with only the exponent and names changed",
           dl.count("identical to OpenAI") == 2))
rcs = (L/"controls.rc").read_text().splitlines() if (L/"controls.rc").exists() else []
ENVFAIL = ("unknown module prefix", "does not exist", "No such file", "unknown package", "file not found", "not found in")
for c in R.get("lean_controls", []):
    ran = any(r.startswith("control-%s rc=" % c["name"]) for r in rcs)
    if c.get("expect", "compiles") == "fails":
        lg = L/("control-%s.log" % c["name"]); txt = lg.read_text() if lg.exists() else ""
        ok.append(("Lean control %s fails, as it must" % c["name"],
                   ran and ("control-%s rc=0" % c["name"]) not in rcs and "error" in txt and not any(s in txt for s in ENVFAIL)))
    else:
        ok.append(("Lean control %s compiles" % c["name"], ("control-%s rc=0" % c["name"]) in rcs))
want = R.get("permitted_axioms", ["propext", "Classical.choice", "Quot.sound"])
for t in R["axiom_targets"]:
    m = re.search(r"'%s' depends on axioms: \[([^\]]*)\]" % re.escape(t), ax)
    ok.append(("axioms: %s only on %s" % (t, ", ".join(want)), bool(m) and set(x.strip() for x in m.group(1).split(",")) <= set(want)))
text = (L/"source-scan.txt").read_text()
words, _, changed = text.partition("# tracked files of the tree changed:")
scan = [l for l in words.splitlines() if l and not l.startswith("#")]
changed, _, lfdiff = changed.partition("# lakefile.lean diff:")
tracked = [l for l in changed.splitlines() if l.strip()]
# a package whose result.json declares a tree_change (u8x, u8xm, u8xmhh: the XExit library line) may change
# lakefile.lean, and only by adding the line(s) that tree_change quotes
# (here: also the lines of "lakefile_add", which apply-overlay.sh appends, with the blank line before them)
tc = R.get("tree_change", ""); lf = [l for l in lfdiff.splitlines() if l.strip() and l.strip() != "+"]
add = [x.strip() for x in R.get("lakefile_add", [])]
for req in R.get("requires", []):   # apply-overlay.sh lays these first, with their own lakefile_add lines
    rq = Path(sys.argv[3]).resolve().parent.parent / req / "result.json"
    if rq.exists():
        add += [x.strip() for x in json.load(open(rq)).get("lakefile_add", [])]
if (tc or add) and tracked and all(l.split()[-1].endswith("lakefile.lean") for l in tracked) and lf \
        and all(l.startswith("+") and l[1:].strip() and (l[1:].strip() in add or l[1:].strip() in tc) for l in lf):
    ok.append(("source scan: lakefile.lean changed only by the lines of tree_change / lakefile_add (this package and those it requires)", True))
    tracked = []
challenge_files = {R["challenge_module"].replace(".", "/") + ".lean", R["negative_challenge_module"].replace(".", "/") + ".lean"}
cs = R.get("control", {}).get("suffix")  # a control instance has its own pair of statement files (as in package_check.py)
if cs:
    challenge_files |= {"Work/Fourier%s/UniformFourierChallenge.lean" % cs, "Work/Fourier%sNeg/UniformFourierChallenge.lean" % cs}
def allowed(l):
    f, rest = l.split(":", 1)
    if f in challenge_files and re.search(r"^\d+:\s*sorry\s*$", rest): return True
    # any copy of OpenAI's statement file (a package may carry the statement files of several instances): its two
    # sorry lines, 327 and 331, as in OpenAI's file
    if f.endswith("/UniformFourierChallenge.lean") and re.search(r"^(327|331):\s*sorry\s*$", rest): return True
    body = rest.split(":", 1)[1] if ":" in rest else rest
    return "No `sorry`" in body or "No `native_decide`" in body or body.lstrip().startswith(("--", "/-", "*", "GENERATED", "Generated"))
hits = [l for l in scan if ":" in l and not allowed(l)]
ok.append(("source scan: no forbidden word outside comments and the challenges' sorry (%d hits)" % len(hits), not hits))
if R.get("mode", "add") == "add":
    ok.append(("source scan: no tracked file of the base tree changed", not tracked))
for name, v in ok: print(("PASS " if v else "FAIL ") + name)
for h in hits: print("   hit: " + h)
print("ALL CHECKS PASSED" if all(v for _, v in ok) else "SOME CHECK FAILED")
sys.exit(0 if all(v for _, v in ok) else 1)
PY
RC=$?
echo "END $(date -u) rc=$RC"
exit $RC
