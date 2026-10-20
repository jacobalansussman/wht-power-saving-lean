#!/bin/bash
# swap-certificate.sh: install another certificate package under community/, make it the headline, re-run the
# list of checks and, only when every check of the new package passed with a log, update the figures and names
# in the text files.
#
# usage: tools/swap-certificate.sh [options] <package folder or .tar.gz>
#   --name <folder>   folder name under community/ (default: "name" of the package's result.json, else its folder)
#   --remove-old      remove the old package folder with git rm (default: keep it; the text may still cite its logs)
#   --dry-run         run the static checks and print the substitutions, change nothing
#   --tree <T>        also lay the package over T (lean/complex588 of a BUILT clone of Ryan Shea's repository at the
#                     commit of result.json) and build and check it there (community/scripts/build-and-check.sh);
#                     for a Colab-sized machine, run inside tmux; logs go to the package's logs/rerun-<time>/
#   --variant <v>     with --tree: a folder of <package>/variants/ to lay over the overlay
#   --comparator      with --tree: also run the official comparator, positive and negative configurations
#                     (community/scripts/run-comparator.sh; needs CC and TC, see that script)
#
# A package is a folder like community/c588m: overlay/ (the files laid over Shea's tree), MODULES.txt,
# result.json, SHA256SUMS, and optionally data/, preflight/, variants/, logs/. A package made by the generator
# (tools/community/c588m-gen/package.py) has no result.json: one is derived from its Lean files, and its title,
# its credits and its check statuses must then be written by hand (the script says so).
#
# What it does, in order:
#   1. static checks of the incoming package as it arrived (tools/community/package_check.py check); stops if one
#      fails;
#   2. computes the old -> new substitutions for the figures and names (saving 1.3552606e-3 and 1.3553e-3 style,
#      the numerators, the theorem suffix, the folder name, the gain over the compared record);
#   3. installs the package as community/<name> and keeps the old one (unless --remove-old);
#   4. re-runs the checks: with --tree, the build, the axioms, Compare.lean positive and negative, the source
#      scan; with --comparator, the official comparator positive and negative; the statuses and the log paths go
#      into result.json; checks that were not run keep the status they had ("not run" for a derived result.json);
#   5. rewrites the package's SHA256SUMS, then re-runs the static checks on the installed copy;
#   6. the text files (README.md VERIFY.md NOTICE ORIGIN.md RELATED-WORK.md notes/**/*.md, or $FIGURE_FILES) are
#      changed ONLY IF every check listed in the new result.json has status "pass" and names a log file that
#      exists in the package (or the package's preflight/PREFLIGHT.json). Otherwise nothing in them is changed and
#      the script lists the checks that block it. When they are changed, community/current.json names the new
#      package, and:
#        - names match whole words only: "community/c588m" does not match "community/c588m-gen" or
#          "tools/community/c588m"; the suffix "C588M" matches "FourierC588M." and "FourierC588MNeg" but not
#          "C588MControls";
#        - a line that names the Walsh-Hadamard corollary ("WHT" or "Walsh") is never changed: those figures
#          come from a separate build;
#        - a note that records the figures of its day carries the HTML comment "<!-- swap-certificate: keep
#          figures" and is left unchanged;
#   7. rewrites the repository's MANIFEST.sha256 and prints what must be read by hand: the title, the credits,
#      and the prose of README.md and VERIFY.md (the check tables cite log files by path).
# Exit status: 0 when everything above was done (or for --dry-run); 7 when the package was installed but the text
# files were left unchanged because a check is not "pass" with a log; other non-zero values when it stopped earlier.
# It never commits and never pushes.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PC="python3 -I $ROOT/tools/community/package_check.py"
NAME=""; REMOVE=0; DRY=0; TREE=""; VAR=""; CMP=0; SRC=""
while [ $# -gt 0 ]; do
  case "$1" in
    --name) NAME=$2; shift 2 ;;
    --remove-old) REMOVE=1; shift ;;
    --keep-old) REMOVE=0; shift ;;  # the default; accepted for old command lines
    --dry-run) DRY=1; shift ;;
    --tree) TREE=$(cd "$2" && pwd); shift 2 ;;
    --variant) VAR=$2; shift 2 ;;
    --comparator) CMP=1; shift ;;
    -h|--help) sed -n '2,48p' "$0"; exit 0 ;;
    *) SRC=$1; shift ;;
  esac
done
[ -n "$SRC" ] || { sed -n '2,48p' "$0"; exit 2; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
if [ -f "$SRC" ]; then
  mkdir -p "$TMP/in" && tar -xzf "$SRC" -C "$TMP/in" 2>/dev/null
  P=$(dirname "$(find "$TMP/in" -maxdepth 3 -type d -name overlay | head -1)")
else
  mkdir -p "$TMP/in" && cp -R "$SRC" "$TMP/in/" && P="$TMP/in/$(basename "$SRC")"
fi
[ -d "$P/overlay" ] || { echo "no overlay/ in $SRC"; exit 3; }
find "$P" -name '.DS_Store' -delete
OLD=$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["headline"])' "$ROOT/community/current.json")
OLDR="$ROOT/community/$OLD/result.json"
echo "== current headline: community/$OLD ($(python3 -I -c 'import json,sys; d=json.load(open(sys.argv[1])); print(d["saving"], d["suffix"])' "$OLDR"))"
DERIVED=0
if [ ! -f "$P/result.json" ]; then
  echo "== no result.json in the package: deriving one from its Lean files (title, credits, statuses to be written by hand)"
  $PC derive "$P" "$OLDR" > /dev/null; DERIVED=1
fi
# 1. the static checks run on the package exactly as it arrived (its SHA256SUMS covers its own result.json)
echo "== 1. static checks of the incoming package"
$PC check "$P" "$ROOT" || { echo "STOPPED: a static check failed; nothing was changed"; exit 4; }
if [ -z "$NAME" ]; then
  NAME=$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1])).get("name",""))' "$P/result.json")
  [ -n "$NAME" ] || NAME=$(basename "$P")
fi
case "$NAME" in *[!A-Za-z0-9._-]*|"") echo "bad package name: $NAME"; exit 6 ;; esac
python3 -I - "$P/result.json" "$NAME" <<'PY'
import json, sys
p, n = sys.argv[1], sys.argv[2]; d = json.load(open(p))
if d.get("name") != n:  # written only when it changes; SHA256SUMS is rewritten in step 5
    d["name"] = n; open(p, "w").write(json.dumps(d, indent=2, ensure_ascii=False) + "\n")
PY
echo "== 2. substitutions (old -> new)"
$PC tokens "$OLDR" "$P/result.json" | tee "$TMP/tokens.json"
FILES=${FIGURE_FILES:-"README.md VERIFY.md NOTICE ORIGIN.md RELATED-WORK.md $(cd "$ROOT" && find notes -name '*.md' 2>/dev/null | tr '\n' ' ')"}
subst() { # subst <write 0/1>
  (cd "$ROOT" && python3 -I - "$TMP/tokens.json" "$1" $FILES <<'PY'
import json, re, sys
from pathlib import Path
pairs = json.load(open(sys.argv[1])); write = sys.argv[2] == "1"
if not pairs: print("   (no substitution)"); sys.exit(0)
def pat(k):
    e = re.escape(k)
    if k[0].isdigit() or k[0] == "+":  # numbers: not part of a longer number
        return r"(?<![0-9.])" + e + r"(?![0-9])"
    if k.startswith("community/"):     # folder: whole path component only
        return r"(?<![A-Za-z0-9_./-])" + e + r"(?![A-Za-z0-9_-])"
    if re.fullmatch(r"[A-Z][A-Za-z0-9]*", k):  # theorem suffix: end of a name, or before Neg
        return e + r"(?=Neg\b|[^A-Za-z0-9_]|$)"
    return r"(?<![A-Za-z0-9_])" + e + r"(?![A-Za-z0-9_])"
rx = re.compile("|".join(pat(k) for k in sorted(pairs, key=len, reverse=True)))
WHT = re.compile(r"WHT|Walsh")
for f in sys.argv[3:]:
    p = Path(f)
    if not p.is_file(): continue
    t = p.read_text()
    if "<!-- swap-certificate: keep figures" in t:  # notes that record the figures of their day
        print(f"   {f}: kept (it carries the marker comment: keep figures)"); continue
    out, n, kept = [], 0, 0
    for line in t.splitlines(keepends=True):
        if WHT.search(line):
            kept += bool(rx.search(line)); out.append(line); continue
        u, k = rx.subn(lambda m: pairs[m.group(0)], line); n += k; out.append(u)
    if n or kept:
        print(f"   {f}: {n} substitutions" + (f"; {kept} Walsh-Hadamard lines left unchanged" if kept else ""))
        if write and n: p.write_text("".join(out))
PY
  ) ; }
if [ $DRY = 1 ]; then echo "== dry run: the substitutions that would be made if every check passes"; subst 0; exit 0; fi
echo "== 3. install community/$NAME"
cd "$ROOT"
if [ -e "community/$NAME" ]; then echo "community/$NAME exists; choose another --name (nothing was changed)"; exit 5; fi
cp -R "$P" "community/$NAME"
if [ $REMOVE = 1 ]; then
  echo "   removing community/$OLD (--remove-old; it stays in the history)"
  git rm -r -q --cached "community/$OLD" 2>/dev/null || true; rm -rf "community/$OLD"
fi
echo "== 4. checks"
NEW="$ROOT/community/$NAME"
setstatus() { # setstatus <log dir relative to the package> id=status ...
  python3 -I - "$NEW/result.json" "$@" <<'PY'
import json, sys
p = sys.argv[1]; logdir = sys.argv[2]; d = json.load(open(p)); st = dict(a.split("=", 1) for a in sys.argv[3:])
logs = {"build": "build.out", "axioms": "axioms.log", "compare": "compare.log", "compare-negative": "compare-negative.log",
        "source-scan": "source-scan.txt", "comparator": "cmp-pos.out", "comparator-negative": "cmp-neg.out"}
ids = {c["id"] for c in d.setdefault("checks", [])}
for c in d["checks"]:
    if c["id"] in st: c["status"] = st[c["id"]]; c["log"] = f"{logdir}/{logs.get(c['id'], '')}"
for k, v in st.items():
    if k not in ids: d["checks"].append({"id": k, "status": v, "what": "", "log": f"{logdir}/{logs.get(k, '')}"})
open(p, "w").write(json.dumps(d, indent=2, ensure_ascii=False) + "\n")
PY
}
if [ -n "$TREE" ]; then
  STAMP=$(date -u +%Y%m%dT%H%MZ); REL="logs/rerun-$STAMP"; export LOGDIR="$NEW/$REL"; mkdir -p "$LOGDIR"
  bash "$ROOT/community/scripts/apply-overlay.sh" "$NEW" "$TREE" $VAR
  if bash "$ROOT/community/scripts/build-and-check.sh" "$NEW" "$TREE" > "$LOGDIR/build.out" 2>&1; then
    setstatus "$REL" build=pass axioms=pass compare=pass compare-negative=pass source-scan=pass
  else
    setstatus "$REL" build="see log" axioms="see log" compare="see log" compare-negative="see log" source-scan="see log"
  fi
  if [ $CMP = 1 ]; then
    POS=$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["comparator"]["positive"])' "$NEW/result.json")
    NEG=$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["comparator"]["negative"])' "$NEW/result.json")
    if bash "$ROOT/community/scripts/run-comparator.sh" "$NAME-pos" "$TREE" "$POS" > "$LOGDIR/cmp-pos.out" 2>&1 \
       && grep -q "Your solution is okay" "$LOGDIR/cmp-pos.out"; then setstatus "$REL" comparator=pass; else setstatus "$REL" comparator=FAIL; fi
    if bash "$ROOT/community/scripts/run-comparator.sh" "$NAME-neg" "$TREE" "$NEG" > "$LOGDIR/cmp-neg.out" 2>&1; then
      setstatus "$REL" comparator-negative="FAIL (accepted)"; else setstatus "$REL" comparator-negative=pass; fi
  fi
else
  echo "   no --tree: no Lean check was run here; the statuses of result.json are kept as they are"
fi
echo "== 5. SHA256SUMS, then the static checks of the installed copy"
$PC sums "$NEW"
$PC check "$NEW" "$ROOT" || echo "WARNING: a static check of the installed copy failed (see above)"
echo "== checks list (community/$NAME/result.json)"
python3 -I - "$NEW" > "$TMP/blocking.txt" <<'PY'
import json, sys
from pathlib import Path
pkg = Path(sys.argv[1]); d = json.load(open(pkg / "result.json")); block = []
for c in d.get("checks", []):
    log = c.get("log", ""); has = bool(log) and (pkg / log).is_file()
    print("   %-22s %-12s %s" % (c["id"], c["status"], log if has else (log + " (MISSING)" if log else "(no log)")), file=sys.stderr)
    if c["status"] != "pass" or not has: block.append(c["id"])
print(" ".join(block))
PY
BLOCK=$(cat "$TMP/blocking.txt")
echo "== 6. figures and names in the text files"
if [ -n "$BLOCK" ]; then
  echo "   NOT CHANGED: these checks are not 'pass' with a log in the package: $BLOCK"
  echo "   the text files and community/current.json still describe community/$OLD; the substitutions they would get are:"
  subst 0
else
  subst 1
  python3 -I - community/current.json "$NAME" "${VAR:-$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1])).get("variant",""))' "community/$NAME/result.json")}" <<'PY'
import json, sys
p = sys.argv[1]; d = json.load(open(p)); d["headline"] = sys.argv[2]; d["variant"] = sys.argv[3]
open(p, "w").write(json.dumps(d, indent=2) + "\n")
PY
  echo "   community/current.json now names $NAME"
fi
echo "== 7. MANIFEST.sha256"
$PC manifest "$ROOT"
echo "== read by hand before committing"
python3 -I -c 'import json,sys
d=json.load(open(sys.argv[1]))
print("   title:   " + d.get("title",""))
print("   credits: " + ("\n            ".join(d.get("credits",[]))))
if d.get("credits_note"): print("   NOTE: " + d["credits_note"])' "$NEW/result.json"
echo "   README.md and VERIFY.md: the check tables cite log files by path, and the prose that describes how the"
echo "   certificate was made (for example the number of source copies removed) is not a figure and was not"
echo "   changed. git diff shows every change made."
if [ -n "$BLOCK" ]; then
  echo "EXIT 7: the package was installed, but the text files and community/current.json were NOT changed (step 6)"
  exit 7
fi
