#!/usr/bin/env python3
"""run.py -- from the published search data to the certificate tools/certificate/gcert1-e8-r783.json.gz.

usage: python3 -B tools/e8/run.py [--work DIR] [--pool] [--solve] [--no-input-check]

Default (standard library only, about 2 minutes):
   data/hierarchy-l0.55-m0.85.json.gz        the recorded output of the solver stage
     -> prune.py     compile, drop the target stops that do not pay
     -> reuse2.py    copies written on dead helpers (slot reuse)
     -> xsub.py      single-source helpers replaced by the source array itself
     -> xexp.py 6    helpers that live only before the scatter replaced by direct reads
     -> sched.py 1 60 1.5     the same additions at other frames (seed 1, 60 annealing sweeps, start temperature 1.5)
   then: sha256 of the result against the published certificate, the repository's two checkers
   (tools/gx/refcheck.py, tools/gx/gxdry.py) and the stand-alone replay (replay.py) as a second opinion.

--pool    first rebuild the pool of 1,470 sums from the four hierarchies in data/pool-sources/ (pool1.py) and
          compare it with data/pool-a.json.
--solve   re-run the solver stage (solve.py data/pool-a.json 0.55 0.85; needs numpy and scipy) and go on from ITS
          hierarchy.  Says whether the solver returned the recorded hierarchy.  With another scipy the solver may
          return another optimum of equal value; the checkers then still decide, but the checksum will differ.
--no-input-check   do not compare the data files with data/SHA256 first.

Exit code 0 only if the JSON written is byte for byte the JSON of the published certificate AND all three checkers
accept it at the figure 8762479.  Anything else ends with a line starting FAILED and exit code 1."""
import sys
sys.dont_write_bytecode = True
import gzip, hashlib, json, os, re, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.dirname(HERE)
TAG = "oh-l0.55-m0.85-s0"          # the working name of the recorded run; the last pass writes its input file's name into the certificate
EXPECT_JSON = "3594f19c4d01cf11fb930c5c61baeed399620acbc5b94096b16bcacb23997dd3"      # sha256 of the uncompressed JSON
EXPECT_GZ = "4b92f00fc7454b9b15e6d71a3b05062792eacca0a0ca4efee46f87827d57c8aa"        # sha256 of the .gz file
EXPECT_FIGURE = 8762479
EXPECT_R = 783
PUBLISHED = TOOLS + "/certificate/gcert1-e8-r783.json.gz"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fail(msg):
    print("FAILED: " + msg, flush=True)
    sys.exit(1)


def stage(name, args, work, keep=4):
    """run one program in the working directory; its whole output goes to logs/<name>.log, the last lines are shown"""
    t0 = time.time()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", GX_OUT=work + "/gxout")
    p = subprocess.run([sys.executable, "-B"] + args, cwd=work, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open("%s/logs/%s.log" % (work, name), "w") as fh:
        fh.write(p.stdout)
    lines = [ln for ln in p.stdout.rstrip("\n").split("\n") if ln.strip()]
    print("== %s  (%.0f s, exit code %d)" % (name, time.time() - t0, p.returncode))
    for ln in (lines[-keep:] if keep else []):
        print("   " + ln[:230])
    sys.stdout.flush()
    if p.returncode != 0:
        fail("%s ended with exit code %d; full output in %s/logs/%s.log" % (name, p.returncode, work, name))
    return p.stdout


def main():
    a = sys.argv[1:]
    work = None
    if "--work" in a:
        work = os.path.abspath(a[a.index("--work") + 1])
    for w in a:
        if w.startswith("--") and w not in ("--work", "--pool", "--solve", "--no-input-check"):
            print(__doc__); sys.exit(2)
    work = work or os.path.abspath("e8-run")
    os.makedirs(work + "/logs", exist_ok=True)
    T0 = time.time()
    print("Python %s; working directory %s" % (sys.version.split()[0], work))

    # ---- 0. the data files
    if "--no-input-check" not in a:
        n = 0
        with open(HERE + "/data/SHA256") as fh:
            for ln in fh:
                if not ln.strip() or ln.startswith("#"): continue
                want, rel = ln.split(None, 1)
                rel = rel.strip()
                try:
                    with open(HERE + "/" + rel, "rb") as f2: have = sha(f2.read())
                except OSError as e:
                    fail("data file %s cannot be read (%s)" % (rel, e))
                if have != want:
                    fail("data file %s has sha256 %s, data/SHA256 says %s" % (rel, have, want))
                n += 1
        print("== data: %d files agree with data/SHA256" % n)

    sys.path.insert(0, HERE)
    import hio
    recorded = HERE + "/data/hierarchy-l0.55-m0.85.json.gz"
    pool = HERE + "/data/pool-a.json"

    # ---- optional: the pool from the four hierarchies of the earlier rounds
    if "--pool" in a:
        src = [HERE + "/data/pool-sources/%s.json.gz" % nm for nm in ("best-grammar", "g-b3-F8uu", "g-b4", "g-b5")]
        stage("pool1", [HERE + "/pool1.py", work + "/pool.json"] + src, work, keep=3)
        same = hio.load_pool(work + "/pool.json") == hio.load_pool(pool)
        print("   the rebuilt pool %s data/pool-a.json" % ("EQUALS" if same else "DIFFERS FROM"))
        if not same: fail("pool1.py did not rebuild data/pool-a.json")
        pool = work + "/pool.json"

    # ---- 1. the hierarchy: recorded, or from the solver
    if "--solve" in a:
        stage("solve", [HERE + "/solve.py", pool, "0.55", "0.85", work + "/hierarchy.json", "30"], work, keep=3)
        same = hio.load_h(work + "/hierarchy.json") == hio.load_h(recorded)
        print("   the solver's hierarchy %s the recorded one" % ("EQUALS" if same else "DIFFERS FROM"))
    else:
        try:
            hio.dump_h(hio.load_h(recorded), work + "/hierarchy.json")
        except Exception as e:
            fail("the recorded hierarchy cannot be read: %r" % (e,))
        print("== hierarchy: the recorded output of the solver stage (data/hierarchy-l0.55-m0.85.json.gz)")

    # ---- 2. the chain of passes
    t, xs, xe, out = "t-" + TAG, "c-%s-xsub.json.gz" % TAG, "c-%s-xexp.json.gz" % TAG, "c-%s-sched.json.gz" % TAG
    for f in (t + ".json.gz", xs, xe, out, "drop.json"):
        if os.path.exists(work + "/" + f): os.remove(work + "/" + f)
    o0 = stage("prune", [HERE + "/prune.py", "hierarchy.json", "drop.json"], work, keep=1)
    m0 = re.search(r"figure (\d+)", o0)
    if not m0 or int(m0.group(1)) == 0:
        fail("the circuit compiled from the hierarchy is refused by the checkers (prune.py reports the figure 0): the hierarchy is not a valid one")
    stage("reuse2", [HERE + "/reuse2.py", t, "maxk=8", "gen", "join", "passes=2", "seeds=10"], work, keep=4)
    stage("xsub", [HERE + "/xsub.py", t + ".json.gz", xs], work, keep=3)
    stage("xexp", [HERE + "/xexp.py", xs, xe, "6"], work, keep=4)
    stage("sched", [HERE + "/sched.py", work + "/" + xe, work + "/" + out, "1", "60", "1.5"], work, keep=4)
    cert = work + "/" + out
    if not os.path.exists(cert):
        fail("the chain wrote no certificate (a pass refused its own result; see %s/logs/)" % work)

    # ---- 3. checksums
    with open(cert, "rb") as fh: gz = fh.read()
    raw = gzip.decompress(gz)
    print("== result %s" % cert)
    print("   sha256 of the .gz file           %s" % sha(gz))
    print("   sha256 of the uncompressed JSON  %s" % sha(raw))
    print("   published certificate (JSON)     %s" % EXPECT_JSON)
    same_json, same_gz = sha(raw) == EXPECT_JSON, sha(gz) == EXPECT_GZ
    print("   the JSON %s the published certificate; the .gz bytes %s" % ("IS" if same_json else "IS NOT", "are the same too" if same_gz else
          "differ (another zlib compresses differently; the JSON is what counts)" if same_json else "differ"))
    if os.path.exists(PUBLISHED):
        with open(PUBLISHED, "rb") as fh: pub = fh.read()
        print("   tools/certificate/gcert1-e8-r783.json.gz: %s" % ("SAME BYTES" if pub == gz else "same JSON inside, other gzip bytes" if gzip.decompress(pub) == raw else "DIFFERENT"))
    else:
        print("   (tools/certificate/gcert1-e8-r783.json.gz is not in this tree; compared with the recorded checksums only)")

    # ---- 4. the repository's checkers, and the stand-alone replay
    o1 = stage("refcheck", [TOOLS + "/gx/refcheck.py", cert], work, keep=0)
    m1 = re.search(r"^ACCEPTED by gx\.check1: h=(\d+) v=(\d+) R=(\d+) N=(\d+)[^;]*", o1, re.M)
    print("   refcheck.py: %s" % (m1.group(0) if m1 else "NOT ACCEPTED: " + o1.strip()[-200:]))
    o2 = stage("gxdry", [TOOLS + "/gx/gxdry.py", cert], work, keep=0)
    m2 = re.search(r"MIRROR ACCEPTS", o2)
    m3 = re.search(r"10 decimals, fill 1 - 2\^-40: whole-block (\d+)[^;]*", o2)
    print("   gxdry.py:    %s; %s" % ("MIRROR ACCEPTS" if m2 else "NO ACCEPTANCE", m3.group(0) if m3 else "no figure"))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    p = subprocess.run([sys.executable, "-I", "-B", HERE + "/replay.py", cert, "e8", "json", "seed=1"], cwd=work, env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open(work + "/logs/replay.log", "w") as fh: fh.write(p.stdout)
    m4 = re.search(r"^REPLAY ACCEPTED: R=(\d+) .* figure=(\d+)  E8 family: True", p.stdout, re.M)
    print("== replay (exit code %d)" % p.returncode)
    print("   replay.py:   %s" % (m4.group(0) if m4 else "NOT ACCEPTED: " + p.stdout.strip()[-300:]))

    bad = []
    if not m1: bad.append("refcheck.py does not accept")
    elif int(m1.group(3)) != EXPECT_R: bad.append("refcheck.py sees %s helpers, not %d" % (m1.group(3), EXPECT_R))
    if not (m2 and m3): bad.append("gxdry.py does not accept")
    elif int(m3.group(1)) != EXPECT_FIGURE: bad.append("gxdry.py gives the figure %s, not %d" % (m3.group(1), EXPECT_FIGURE))
    if not m4 or p.returncode != 0: bad.append("replay.py does not accept")
    elif int(m4.group(2)) != EXPECT_FIGURE: bad.append("replay.py gives the figure %s, not %d" % (m4.group(2), EXPECT_FIGURE))
    if not same_json: bad.append("the result is not the published certificate (other sha256 of the JSON)")
    print("total %.0f s" % (time.time() - T0))
    if bad:
        fail("; ".join(bad))
    print("RESULT: REPRODUCED.  The chain wrote the JSON of the published certificate byte for byte (sha256 %s);" % EXPECT_JSON)
    print("        refcheck.py ACCEPTED, gxdry.py and replay.py give the figure %d (a saving of 8.7625e-4), %d helpers." % (EXPECT_FIGURE, EXPECT_R))
    print("        This is a price computed by the Python tools, not a proof.  For what Lean checks about this file, see VERIFY.md.")


if __name__ == "__main__":
    main()
