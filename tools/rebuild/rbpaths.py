"""rbpaths.py: the directories and the outside input files of the rebuild scripts in tools/rebuild/.

  INPUTS  the directory that holds the outside input files, laid out as fetch_inputs.py writes them:
          INPUTS/pr168/<path in the outside repository>     the eight files at the commit of pull request #168
          INPUTS/pr193/<path in the outside repository>     the one file at the commit of pull request #193
          INPUTS/pr233/<path in the outside repository>     optional: one file at the commit of pull request #233
  WORK    the work directory: WORK/cache/ (intermediate results, Python pickle files written and read by
          these scripts only) and WORK/logs/
The step scripts (build.py, aligned.py, rewrite.py, phys.py, emit.py) take the two directories from the
environment variables REBUILD_INPUTS and REBUILD_WORK.  rebuild.py sets both.

"The outside repository" is github.com/CrocSwap/integer-mult-bounds.  FILES lists every file of it that the
rebuild opens: group, path, size in bytes, sha256, and what the file is used for.  All nine are JSON and are
loaded with json.load.  No program of the outside repository is executed or imported.
OPTIONAL lists one more file in the same way.  The rebuild of the published certificate does not read it;
rebuild.py reads it only when it is named with pairs=<file>.
"""
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE)) + "/"
HUB = "CrocSwap/integer-mult-bounds"
COMMIT = {"pr168": "4a3c769e5c5430e7114c4d3e099ff34664677f17",      # head of pull request #168 when the data were taken
          "pr193": "187e1010ac8b259af8e9b5166f68b64bc27b4b47",      # head of pull request #193 when the data were taken
          "pr233": "109a857a329d18ed5552d5573f17ddfa886ae57f"}      # head of pull request #233 when the file was taken

FILES = [
    ("pr168", "references/paired-cube/sources/local_L1.json", 185,
     "1861724ee3b58ea8acb3334416fbfd9172a9ae02ad11a6f6885a4909f8692f05", "local configuration of a cube"),
    ("pr168", "references/paired-cube/sources/tmod_TE_TD_TB3_1_1_4_full_6.0617964e-4.json", 30417,
     "9dcdb6f4e33410418c01f55f092862edfa7c05001ceada3974f8d1196bcd894b", "triple module"),
    ("pr168", "references/paired-cube/sources/pmod_J0_full_6.0666810e-4.json", 4424,
     "7d1ae23430e110d32aa1cab84415abdbc918dfab6b317cb24d30756fc90c5925", "pair module"),
    ("pr168", "references/paired-cube/sources/qmod_climb3u_best.json", 412,
     "3f85f6b2134af8781f16229186d0d7c1de71ab4d5235244259e1b79f2e3d14ba", "all-but-one module"),
    ("pr168", "references/paired-cube/selected-module/matching-arcs.json", 146101,
     "46c77eb75ebcb2435595a3b6522a54e1d12d3c61eecb3b51c489edd1f78447e1", "carrier arcs"),
    ("pr168", "references/paired-cube/physical/frames.json", 603959,
     "3681fc6ef151cdfa5453512f56f7165607ece64ccf6b71478af8c327350f5800", "moved frames"),
    ("pr168", "references/paired-cube/selected-module/SOURCE.json", 340,
     "7689968c429109119a232c3f7631cdfa5b7a7895c61af36ee82e3229988e45ab",
     "comparison only: the graph hash recorded by the outside repository"),
    ("pr168", "certificates/paired-cube-complex-input.json", 1235,
     "b1898950f95a7ce618f376d14e628061b345a519e7bcac76b67d277deeec9af5",
     "comparison only: the profile of the #168 word recorded by the outside repository"),
    ("pr193", "research/source-assisted-v4/data/physical-pairs.json", 42931,
     "7095501f054e9ec18aa212bfed1b01e497ca2ae997f4c9ece2779909b3fd72e2", "hand-over pairs"),
]

OPTIONAL = [
    ("pr233", "research/source-assisted-v4-layer/witness/physical-pairs.json", 44022,
     "e6409fffbf5f3b1ac6434acf2e426503de3a8a18537b8fee761fc4fc931fb896",
     "optional: another choice of the hand-over pairs for the same word"),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def slash(d):
    return os.path.abspath(d).rstrip("/") + "/"


def verify(inputs, quiet=False, files=None):
    """compare every input file under `inputs` with its sha256 in FILES (or in the list `files`); returns the
    number of files that are missing or different (0 = all nine are the expected files)"""
    inputs, bad = slash(inputs), 0
    for group, path, size, want, _ in (FILES if files is None else files):
        p = inputs + group + "/" + path
        if not os.path.isfile(p):
            state = "MISSING"
        else:
            state = "ok" if sha256(p) == want else "DIFFERENT"
        bad += state != "ok"
        if not quiet:
            print("  %-9s %s  %s/%s" % (state, want[:16], group, path))
    return bad


def _env(var):
    d = os.environ.get(var)
    if not d:
        sys.exit("%s is not set.  Run tools/rebuild/rebuild.py, which sets it; its first lines say how." % var)
    return slash(d)


def inputs():
    """INPUTS, with a trailing slash"""
    return _env("REBUILD_INPUTS")


def work():
    """WORK, with a trailing slash; WORK/cache and WORK/logs are created"""
    w = _env("REBUILD_WORK")
    for d in ("cache", "logs"):
        os.makedirs(w + d, exist_ok=True)
    return w
