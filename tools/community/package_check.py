#!/usr/bin/env python3
"""Static checks and bookkeeping for the certificate packages under community/ (no Lean needed).

usage:
  package_check.py check  <package dir> <repository root>   static checks; exit 1 if one fails
  package_check.py derive <package dir> <old result.json>   write <package>/result.json if it has none
  package_check.py sums   <package dir>                     rewrite <package>/SHA256SUMS
  package_check.py tokens <old result.json> <new result.json>   print the figure substitutions as JSON
  package_check.py manifest <repository root>              rewrite MANIFEST.sha256 (every file git would commit)
  package_check.py same   <file a> <file b>                 exit 0 if the two files are equal once the
                                                            modification notice (below) is removed from each

A package is a folder with overlay/ (files laid over lean/complex588 of Ryan Shea's repository),
MODULES.txt (new modules in build order), result.json (figures, names, checks) and SHA256SUMS.
The static checks are:
  1. every file matches SHA256SUMS;
  2. the challenge is OpenAI's UniformFourier.lean (third-party/openai-math) with only the exponent line and the
     seven suffixed names changed (blank lines ignored);
  3. the negative challenge is the challenge with only the exponent numerator raised by one;
  4. the exponent in Goal.lean and in the challenge is the numerator of result.json;
  5. the comparator configurations name the challenge, the solution and the theorems of result.json;
  6. MODULES.txt lists every .lean file of overlay/ exactly once (mode "add");
  7. no sorry, admit, axiom, native_decide, implemented_by, extern, unsafe or opaque outside comments, except
     the two `sorry` of each challenge file;
  8. the figures of result.json agree with its numerators;
  9. every check listed in result.json is "pass" and names a log file inside the package, or is "not shown" or
     "not run".
Check 2 also requires the modification notice that each copy of OpenAI's statement file carries (Apache License
2.0, section 4(b)): a comment that starts on line 3 with NOTICE_MARK and ends with "-/", followed by one blank line.
It was added after the build; checks 2 and 3 compare the files without it, and the notice gives the sha256 of the
file as compiled.
A Walsh-Hadamard package (result.json with "results", one entry per theorem) gets checks 1, 5-7 and 9, and per
theorem: the negative challenge is the challenge with only the numerator raised by one (and its module name in the
header comment), and the figures agree with the numerator.
"""
import hashlib, json, re, sys
from pathlib import Path

NAMES = ["decimalExponent", "decimalTime", "TimeBounds", "DFTGoal", "ConvGoal", "transform_main", "convolution_main"]
NOTICE_MARK = "/- Modified copy (Apache License 2.0, section 4(b)) of lean/ComparatorChallenges/UniformFourier.lean"
STATUSES = ("pass", "not shown", "not run")


def strip_notice(t):
    """(text without the modification notice, whether it had one)."""
    ls = t.split("\n")
    if len(ls) > 3 and ls[2] == NOTICE_MARK:
        j = next(i for i in range(2, len(ls)) if ls[i].endswith("-/"))
        if ls[j + 1] == "":
            return "\n".join(ls[:2] + ls[j + 2:]), True
    return t, False


def check_statuses(pkg, checks, label=""):
    bad = []
    for c in checks:
        st, log = c.get("status", ""), c.get("log", "")
        if st == "pass":
            if not log or not (pkg / log).is_file():
                bad.append(f"{label}{c.get('id')}: pass without a log file in the package ({log or 'no log'})")
        elif st not in STATUSES:
            bad.append(f"{label}{c.get('id')}: status {st!r}")
    return bad


WORDS = re.compile(r"\b(sorry|admit|axiom|native_decide|implemented_by|extern|unsafe|opaque|sorryAx)\b")


def sci(n, digits):
    s = f"{n / 10**10:.{digits}e}"
    m, e = s.split("e")
    return f"{m}e{int(e)}"


def mod_path(m):
    return m.replace(".", "/") + ".lean"


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_sums(pkg):
    pkg = Path(pkg)
    files = sorted(p for p in pkg.rglob("*") if p.is_file() and p.name != "SHA256SUMS" and ".DS_Store" not in p.name)
    (pkg / "SHA256SUMS").write_text("".join(f"{sha256(p)}  {p.relative_to(pkg).as_posix()}\n" for p in files))
    return len(files)


def derive(pkg, old):
    """result.json for a package in m-gen's format (no result.json): names and numerators read from the Lean files."""
    pkg = Path(pkg); ov = pkg / "overlay"; old = json.loads(Path(old).read_text())
    ch = [p for p in ov.glob("Work/Fourier*/UniformFourierChallenge.lean") if not p.parent.name.endswith("Neg")]
    if len(ch) != 1:
        sys.exit(f"cannot derive: {len(ch)} challenge files under overlay/Work/Fourier*/")
    s = ch[0].parent.name[len("Fourier"):]
    num = int(re.search(r"def decimalExponent%s : ℝ := 1 - (\d+)/\(10\^\(10:ℕ\)\)" % re.escape(s), ch[0].read_text()).group(1))
    negp = ov / f"Work/Fourier{s}Neg/UniformFourierChallenge.lean"
    neg = int(re.search(r":= 1 - (\d+)/", negp.read_text()).group(1)) if negp.exists() else num + 1
    ax = ov / f"Work/Fourier{s}/Axioms.lean"
    targets = re.findall(r"#print axioms (\S+)", ax.read_text()) if ax.exists() else []
    ctrl = [f"Work.GCert.Data.{s}Controls"] if (ov / f"Work/GCert/Data/{s}Controls.lean").exists() else []
    cmp = sorted(p.name for p in ov.glob("comparator/*.json"))
    r = {
        "format": "community-certificate-package/1", "name": pkg.name, "suffix": s,
        "title": "(to be written)", "statement": old.get("statement", ""), "variant": "",
        "numerator": num, "kernel_numerator": num + 1, "negative_numerator": neg,
        "saving": sci(num, 7), "saving_short": sci(num, 4), "exponent": f"1 - {num}/10^10",
        "compared_with": old.get("compared_with", {}), "mode": "add", "base": old["base"],
        "theorems": [f"OAI.PowerSaving.transform_main{s}", f"OAI.PowerSaving.convolution_main{s}"],
        "axiom_targets": targets, "permitted_axioms": ["propext", "Classical.choice", "Quot.sound"],
        "challenge_module": f"Work.Fourier{s}.UniformFourierChallenge",
        "negative_challenge_module": f"Work.Fourier{s}Neg.UniformFourierChallenge",
        "solution_module": f"Work.Fourier{s}.Main", "axioms_file": f"Work/Fourier{s}/Axioms.lean",
        "goal_file": f"Work/Fourier{s}/Goal.lean", "exponent_definition": f"OAI.PowerSaving.decimalExponent{s}",
        "build_targets": [f"Work.Fourier{s}.Main", f"Work.Fourier{s}.UniformFourierChallenge",
                          f"Work.Fourier{s}Neg.UniformFourierChallenge"] + ctrl + [f"Work.Fourier{s}.Axioms"],
        "comparator": {"positive": next((f"comparator/{c}" for c in cmp if "neg" not in c.lower()), ""),
                       "negative": next((f"comparator/{c}" for c in cmp if "neg" in c.lower()), "")},
        "checks": [{"id": i, "status": "not run", "what": "", "log": ""} for i in
                   ["build", "axioms", "compare", "compare-negative", "comparator", "comparator-negative", "source-scan"]],
        "credits": old.get("credits", []),
        "credits_note": "copied from the package this one replaces: CHECK THEM BY HAND",
    }
    (pkg / "result.json").write_text(json.dumps(r, indent=2, ensure_ascii=False) + "\n")
    return r


def check(pkg, root):
    pkg = Path(pkg); root = Path(root); ov = pkg / "overlay"
    R = json.loads((pkg / "result.json").read_text())
    if "results" in R and "numerator" not in R:
        return check_wht(pkg, root, R)
    res = []
    # 1
    bad = []
    if (pkg / "SHA256SUMS").exists():
        for line in (pkg / "SHA256SUMS").read_text().splitlines():
            h, f = line.split("  ", 1)
            if not (pkg / f).exists() or sha256(pkg / f) != h:
                bad.append(f)
        res.append((f"1. SHA256SUMS ({len(bad)} files differ)", not bad))
    else:
        res.append(("1. SHA256SUMS present", False))
    s = R["suffix"]; n = R["numerator"]; neg = R["negative_numerator"]
    chf = ov / mod_path(R["challenge_module"]); ngf = ov / mod_path(R["negative_challenge_module"])
    oai = (root / "third-party/openai-math/lean/ComparatorChallenges/UniformFourier.lean").read_text()
    if chf.exists():
        t, has_n = strip_notice(chf.read_text())
        line = f"def decimalExponent{s} : ℝ := 1 - {n}/(10^(10:ℕ))"
        ok_line = t.count(line) == 1
        u = t.replace(line, "def decimalExponent : ℝ := 1 - 1/(10^(13:ℕ))")
        for nm in NAMES:
            u = re.sub(r"\b%s%s\b" % (nm, re.escape(s)), nm, u)
        strip = lambda x: [l for l in x.splitlines() if l.strip()]
        res.append(("2. challenge = OpenAI's statement with only the exponent and seven names changed", ok_line and strip(u) == strip(oai)))
        if ngf.exists():
            v, has_m = strip_notice(ngf.read_text())
            res.append(("2b. both statement files carry the modification notice", has_n and has_m))
            res.append(("3. negative challenge = challenge with numerator %d -> %d only" % (n, neg),
                        neg == n + 1 and v == t.replace(line, line.replace(f"1 - {n}/", f"1 - {neg}/"))))
        else:
            res.append(("3. negative challenge present", False))
        g = ov / R.get("goal_file", "")
        res.append(("4. Goal.lean states the same exponent", g.is_file() and line in g.read_text()))
    else:
        res.append(("2. challenge file present: " + str(chf.relative_to(pkg)), False))
    # 5
    okc = True
    for k, want_ch in (("positive", R["challenge_module"]), ("negative", R["negative_challenge_module"])):
        p = ov / R.get("comparator", {}).get(k, "-")
        if not p.is_file():
            okc = False; continue
        c = json.loads(p.read_text())
        okc &= c.get("challenge_module") == want_ch and c.get("solution_module") == R["solution_module"] \
            and c.get("theorem_names") == R["theorems"]
    res.append(("5. comparator configurations name the challenge, solution and theorems", okc))
    # 6
    if R.get("mode", "add") == "add":
        mods = [l.split()[0] for l in (pkg / "MODULES.txt").read_text().splitlines() if l.strip() and not l.startswith("#")]
        leans = sorted(p.relative_to(ov).as_posix()[:-5].replace("/", ".") for p in ov.rglob("*.lean"))
        res.append((f"6. MODULES.txt lists the {len(leans)} overlay modules once each", sorted(mods) == leans and len(set(mods)) == len(mods)))
    # 7
    chal = {mod_path(R["challenge_module"]), mod_path(R["negative_challenge_module"])}
    cs = R.get("control", {}).get("suffix")  # a control instance has its own pair of statement files
    if cs:
        chal |= {f"Work/Fourier{cs}/UniformFourierChallenge.lean", f"Work/Fourier{cs}Neg/UniformFourierChallenge.lean"}
    hits = scan_words(ov, chal)
    res.append((f"7. no forbidden word in Lean code ({len(hits)} hits)", not hits))
    # 8
    res.append(("8. figures agree with the numerators", R["saving"] == sci(n, 7) and R["saving_short"] == sci(n, 4)
                and R["exponent"] == f"1 - {n}/10^10"))
    # 9
    sbad = check_statuses(pkg, R.get("checks", []))
    res.append((f"9. every check is pass with a log in the package, or not shown / not run ({len(sbad)} not)", not sbad))
    return report(res, hits + sbad, bad)


def scan_words(ov, chal):
    hits = []
    for p in sorted(ov.rglob("*.lean")):
        rel = p.relative_to(ov).as_posix(); depth = 0
        for i, l in enumerate(p.read_text().splitlines(), 1):
            code = l
            if depth == 0 and "/-" not in l:
                code = l.split("--", 1)[0]
            opened, closed = l.count("/-"), l.count("-/")
            in_comment = depth > 0 or (opened and l.strip().startswith("/-"))
            depth = max(0, depth + opened - closed)
            if in_comment:
                continue
            if WORDS.search(code):
                if rel in chal and l.strip() == "sorry":
                    continue
                hits.append(f"{rel}:{i}: {l.strip()}")
    return hits


def check_wht(pkg, root, R):
    """Static checks of a Walsh-Hadamard package (one result per theorem under "results")."""
    ov = pkg / "overlay"; res = []
    bad = []
    if (pkg / "SHA256SUMS").exists():
        for line in (pkg / "SHA256SUMS").read_text().splitlines():
            h, f = line.split("  ", 1)
            if not (pkg / f).exists() or sha256(pkg / f) != h:
                bad.append(f)
        res.append((f"1. SHA256SUMS ({len(bad)} files differ)", not bad))
    else:
        res.append(("1. SHA256SUMS present", False))
    chal = set(); sbad = []
    for r in R["results"]:
        n = r["numerator"]; neg = r["negative_numerator"]; th = r["theorem"]
        ch, ng = r["modules"][0], r["modules"][1]
        chal |= {mod_path(ch), mod_path(ng)}
        chf, ngf = ov / mod_path(ch), ov / mod_path(ng)
        ok = chf.is_file() and ngf.is_file() and neg == n + 1
        if ok:
            t, v = chf.read_text(), ngf.read_text()
            ok = f"1 - {n}/(10:ℝ)^10" in t and \
                v == t.replace(str(n), str(neg)).replace(f"`{ch}`", f"`{ng}`")
        res.append((f"3. {th}: negative challenge = challenge with {n} -> {neg} only", ok))
        okc = True
        for k, want_ch in (("positive", ch), ("positive_nanoda", ch), ("negative", ng)):
            p = ov / r.get("comparator", {}).get(k, "-")
            if not p.is_file():
                okc = k == "positive_nanoda" and k not in r.get("comparator", {}); continue
            c = json.loads(p.read_text())
            okc &= c.get("challenge_module") == want_ch and c.get("theorem_names") == [th]
        res.append((f"5. {th}: comparator configurations name its challenge and theorem", okc))
        res.append((f"8. {th}: figures agree with the numerator", r["saving"] == sci(n, 7)))
        sbad += check_statuses(pkg, r.get("checks", []), th.split(".")[-1] + " ")
    if (pkg / "MODULES.txt").exists():
        mods = [l.split()[0] for l in (pkg / "MODULES.txt").read_text().splitlines() if l.strip() and not l.startswith("#")]
        leans = sorted(p.relative_to(ov).as_posix()[:-5].replace("/", ".") for p in ov.rglob("*.lean"))
        res.append((f"6. MODULES.txt lists the {len(leans)} overlay modules once each", sorted(mods) == leans and len(set(mods)) == len(mods)))
    hits = scan_words(ov, chal)
    res.append((f"7. no forbidden word in Lean code ({len(hits)} hits)", not hits))
    res.append((f"9. every check is pass with a log in the package, or not shown / not run ({len(sbad)} not)", not sbad))
    return report(res, hits + sbad, bad)


def report(res, hits, bad):
    for name, v in res:
        print(("PASS " if v else "FAIL ") + name)
    for h in hits[:20]:
        print("   hit: " + h)
    for f in bad[:20]:
        print("   sha256 differs: " + f)
    ok = all(v for _, v in res)
    print("STATIC CHECKS PASSED" if ok else "SOME STATIC CHECK FAILED")
    return ok


def tokens(old, new):
    """Old -> new strings for the figures and names that the text files quote."""
    o = json.loads(Path(old).read_text()); n = json.loads(Path(new).read_text())
    pairs = {}
    def add(a, b):
        if a and b and str(a) != str(b):
            pairs[str(a)] = str(b)
    add(o["saving"], n["saving"]); add(o["saving_short"], n["saving_short"])
    add(sci(o["kernel_numerator"], 7), sci(n["kernel_numerator"], 7))
    add(o["numerator"], n["numerator"]); add(o["kernel_numerator"], n["kernel_numerator"])
    add(o["negative_numerator"], n["negative_numerator"])
    add(o["suffix"], n["suffix"]); add("community/" + o["name"], "community/" + n["name"])
    c = o.get("compared_with", {}).get("numerator")
    if c:
        po = (o["numerator"] / c - 1) * 100; pn = (n["numerator"] / c - 1) * 100
        add(f"+{po:.1f} %", f"+{pn:.1f} %"); add(f"+{po:.1f}%", f"+{pn:.1f}%")
        add(f"{o['numerator'] / c:.3f} times", f"{n['numerator'] / c:.3f} times")
    return pairs


def manifest(root):
    import subprocess
    root = Path(root)
    out = subprocess.run(["git", "-C", str(root), "ls-files", "-co", "--exclude-standard", "-z"],
                         check=True, capture_output=True).stdout.decode().split("\0")
    files = sorted({f for f in out if f and f != "MANIFEST.sha256" and (root / f).is_file()}, key=lambda s: s.encode())
    (root / "MANIFEST.sha256").write_text("".join(f"{sha256(root / f)}  {f}\n" for f in files))
    return len(files)


if __name__ == "__main__":
    if len(sys.argv) < 3 and not (len(sys.argv) == 2 and sys.argv[1] == "manifest"):
        sys.exit(__doc__)
    cmd = sys.argv[1]
    if cmd == "check":
        sys.exit(0 if check(sys.argv[2], sys.argv[3]) else 1)
    elif cmd == "derive":
        print(json.dumps(derive(sys.argv[2], sys.argv[3]), indent=2, ensure_ascii=False))
    elif cmd == "sums":
        print(f"{write_sums(sys.argv[2])} files in SHA256SUMS")
    elif cmd == "manifest":
        print(f"{manifest(sys.argv[2])} files in MANIFEST.sha256")
    elif cmd == "same":
        a = strip_notice(Path(sys.argv[2]).read_text())[0]; b = strip_notice(Path(sys.argv[3]).read_text())[0]
        sys.exit(0 if a == b else 1)
    elif cmd == "tokens":
        print(json.dumps(tokens(sys.argv[2], sys.argv[3]), indent=2, ensure_ascii=False))
    else:
        sys.exit(__doc__)
