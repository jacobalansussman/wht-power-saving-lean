import json, re, sys
from pathlib import Path
L = Path(sys.argv[1]); neg = int(sys.argv[2]); R = json.load(open(sys.argv[3])); ov = Path(sys.argv[4])
ok = []
ok.append(("compare: RESULT: PASS", "RESULT: PASS" in (L/"compare.log").read_text()))
bad = (L/"compare-negative.log").read_text()
ok.append(("compare-negative: fails (rc != 0)", neg != 0))
ed = R.get("exponent_definition", "").split(".")[-1]
ok.append(("compare-negative: the exponent definition differs (%s)" % ed, "differs" in bad and ed in bad))
ax = (L/"axioms.log").read_text()
ok.append(("axioms: no sorryAx", "sorryAx" not in ax))
dl = (L/"derive.log").read_text()
ok.append(("statement files: both identical to OpenAI's statement with only the exponent and names changed",
           dl.count("identical to OpenAI") == 2))
rcs = (L/"controls.rc").read_text().splitlines() if (L/"controls.rc").exists() else []
for c in R.get("lean_controls", []):
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
if (tc or add) and tracked and all(l.split()[-1].endswith("lakefile.lean") for l in tracked) and lf \
        and all(l.startswith("+") and l[1:].strip() and (l[1:].strip() in add or l[1:].strip() in tc) for l in lf):
    ok.append(("source scan: lakefile.lean changed only by the line that result.json's tree_change names", True))
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
