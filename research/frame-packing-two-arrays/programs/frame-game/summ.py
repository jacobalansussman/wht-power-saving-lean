"""summ.py: aggregate design.py outputs: designs by final status pair, and the trend of orphan counts."""
import sys, re, glob, ast, collections
files = sys.argv[1:]
for f in files:
    txt = open(f).read().splitlines()
    head = [l for l in txt if l.startswith("design.py") or l.startswith("m=")]
    cnt = collections.Counter(); n = 0
    for l in txt:
        l = l.strip()
        if l.startswith("(") and l.rstrip()[-1].isdigit():
            try:
                k, v = l.rsplit(" ", 1); k = ast.literal_eval(k); v = int(v)
            except Exception: continue
            st = k[-1] if isinstance(k, tuple) else k
            cls = []
            pl = [x for x in st if str(x).startswith("psilin")]
            st = [x for x in st if not str(x).startswith("psilin")]
            for s in st:
                s = str(s)
                cls.append("cap" if s.startswith("cap") else ("unproduced" if "target" in s else ("inconsistent" if "linear" in s else s)))
            cnt[tuple(cls) + tuple(pl)] += v; n += v
        elif l.startswith("'skel-fail'") or "skel-fail" in l:
            pass
    cand = sum(1 for l in txt if "CANDIDATE" in l)
    print(f"{f}: {' | '.join(head)}")
    print(f"   designs {n}; unit candidates {cand}; by final status: {dict(cnt.most_common())}")
