#!/usr/bin/env python3
"""m-gen: read the frame table `tab` and the label states `st0`, `st2` (SSC.Trie literals) back from the GENERATED
Lean data module (text parse only) and evaluate, with Lean's semantics (key bit 0 first: left = even keys;
GLab.cdim = code & 255), the kernel facts that the glue states about them:
  res_ok : forall q < R, 9 - cdim (tab.get (st0.get (2 v + q))) = schedRes[q]   (schedRes from the price schedule)
  cen_ok : every retained register r: cdim (tab.get (st2.get r)) = 8
Also: st0 = the certificate's start frames, tab codes = gxbind.code of the frames (gx-labels' code).
usage: trie_check.py <S0.lean> <cert.json.gz> <price.schedule.json.gz> <toolkit/tools/gx>"""
import sys, re, json, gzip

def parse_tries(path):
    env = {}
    for line in open(path, encoding="utf-8"):
        m = re.match(r"noncomputable def ([A-Za-z0-9_]+) : SSC\.Trie := (.*)$", line.rstrip("\n"))
        if m: env[m.group(1)] = m.group(2)
    tok = re.compile(r"\s*(\.node|\.leaf|\(|\)|\d+|[A-Za-z_][A-Za-z0-9_]*)")
    memo = {}
    def ev(name):
        if name in memo: return memo[name]
        s = env[name]; toks = []; pos = 0
        while pos < len(s):
            m = tok.match(s, pos)
            if not m: raise ValueError(s[pos:pos + 30])
            toks.append(m.group(1)); pos = m.end()
        i = [0]
        def term():
            t = toks[i[0]]; i[0] += 1
            if t == "(":
                x = term(); assert toks[i[0]] == ")"; i[0] += 1; return x
            if t == ".leaf":
                x = int(toks[i[0]]); i[0] += 1; return [x]
            if t == ".node":
                l = term(); r = term(); out = [None] * (len(l) + len(r))
                out[0::2] = l; out[1::2] = r; return out
            return ev(t)
        v = term(); assert i[0] == len(toks); memo[name] = v; return v
    return env, ev

def main(s0, cert, sched, gx):
    sys.path.insert(0, gx)
    import gxbind as gb
    env, ev = parse_tries(s0)
    c = json.loads(gzip.open(cert, "rt").read())
    s = json.loads(gzip.decompress(open(sched, "rb").read()))
    v, R, h = c["v"], c["R"], c["h"]
    tab, st0, st2 = ev("tab"), ev("st0"), ev("st2")
    cdim = lambda code: code & 255
    out = {}
    out["tab_codes_equal_frames"] = all(tab[i] == gb.code(h, f) for i, f in enumerate(c["frames"])) and \
        all(x == 0 for x in tab[len(c["frames"]):])
    out["st0_equals_start"] = st0[:2 * v + R] == c["start"] and all(x == 0 for x in st0[2 * v + R:])
    res = [h - cdim(tab[st0[2 * v + q]]) for q in range(R)]
    out["res_ok"] = res == s["widths"]
    out["cen_ok"] = all(cdim(tab[st2[r]]) == 8 for _, r, _ in c["ret"])
    print(json.dumps(out))
    assert all(out.values())

if __name__ == "__main__":
    main(*sys.argv[1:5])
