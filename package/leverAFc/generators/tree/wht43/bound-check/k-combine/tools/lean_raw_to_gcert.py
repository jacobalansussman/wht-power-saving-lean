#!/usr/bin/env python3
"""x-mutable: rebuild the gcert/1 JSON of Ryan Shea's complex588 certificate from his GENERATED Lean data
(Work/GCert/Data/Gen/C588/S0.lean, record `raw : GXD.Raw`, produced by Sussman's gxgen.py from
pr588-atm3-source-copies-no-cleanup.gcert.json.gz).  Pure text parsing; nothing of the input is executed.
usage: lean_raw_to_gcert.py <S0.lean> <out.json.gz>
"""
import sys, re, json, gzip

TOK = re.compile(r"\s*(\+\+|\.out|\.inn|\.table|\.star|\.node|\.leaf|[\[\](),⟨⟩]|\d+|[A-Za-z_][A-Za-z0-9_.]*)")

def tokens(s):
    pos, out = 0, []
    s = s.strip()
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m: raise ValueError("cannot tokenise at %r" % s[pos:pos+40])
        out.append(m.group(1)); pos = m.end()
    return out

class P:
    def __init__(self, toks, env): self.t, self.i, self.env = toks, 0, env
    def peek(self): return self.t[self.i] if self.i < len(self.t) else None
    def eat(self, x=None):
        tok = self.t[self.i]
        if x is not None and tok != x: raise ValueError("expected %s got %s" % (x, tok))
        self.i += 1; return tok
    def expr(self):
        a = self.atom()
        while self.peek() == "++":
            self.eat(); b = self.atom(); a = a + b
        return a
    def atom(self):
        tok = self.peek()
        if tok == "[":
            self.eat(); xs = []
            while self.peek() != "]":
                xs.append(self.expr())
                if self.peek() == ",": self.eat()
            self.eat("]"); return xs
        if tok == "(":
            self.eat(); xs = [self.expr()]
            while self.peek() == ",": self.eat(); xs.append(self.expr())
            self.eat(")"); return xs[0] if len(xs) == 1 else tuple(xs)
        if tok == "⟨":
            self.eat(); xs = [self.expr()]
            while self.peek() == ",": self.eat(); xs.append(self.expr())
            self.eat("⟩"); return ("anon",) + tuple(xs)
        if tok == ".out":
            self.eat(); f = int(self.eat()); s = int(self.eat()); tg = self.atom(); ex = self.atom()
            return ("out", f, s, tg, ex)
        if tok == ".inn":
            self.eat(); f = int(self.eat()); t = int(self.eat()); sr = self.atom()
            return ("inn", f, t, sr)
        if tok == ".table":
            self.eat(); return ("table", self.atom())
        if tok.isdigit(): self.eat(); return int(tok)
        if tok in ("true", "false"): self.eat(); return tok == "true"
        self.eat()
        if tok not in self.env: raise KeyError(tok)
        return self.env[tok]

def co(c):
    _, neg, num, den = c
    return (-num if neg else num, den)

def main(src, dst):
    env = {}
    for line in open(src, encoding="utf-8"):
        m = re.match(r"noncomputable def ([A-Za-z0-9_]+) : (.*?) := (.*)$", line.rstrip("\n"))
        if not m: continue
        name, typ, body = m.groups()
        if "Trie" in typ: continue          # label tables / replay states are not needed here
        env[name] = P(tokens(body), env).expr()
    flat = lambda xs: [y for x in xs for y in x]
    frames = flat(env["frames"])
    def gate(g):
        if g[0] == "out":
            _, f, s, tg, ex = g
            return ["out", f, s, [[t, *co(c)] for t, c in tg], list(ex)]
        _, f, t, sr = g
        return ["in", f, t, [[s, *co(c)] for s, c in sr]]
    A = [gate(g) for g in flat(env["gA"])]
    B = [gate(g) for g in flat(env["gB"])]
    rows = env["scat"][1]
    cert = dict(format="gcert/1", h=9, v=len(env["ports"]), R=len(env["start"]) - 2 * len(env["ports"]),
                ports=env["ports"], frames=frames, start=env["start"], final=env["final"],
                ext=[[r, f] for r, f in env["ext"]], ret=[[k, r, f] for k, (r, f) in enumerate(env["ret"])],
                scat=dict(table=[[[k, *co(c)] for k, c in row] for row in rows]), A=A, B=B)
    gzip.open(dst, "wt").write(json.dumps(cert))
    print("v", cert["v"], "R", cert["R"], "frames", len(frames), "A", len(A), "B", len(B), "ret", len(cert["ret"]),
          "ext", len(cert["ext"]))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
