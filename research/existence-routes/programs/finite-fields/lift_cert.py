"""lift_cert.py -- ex2-modp, refutation side (PROOF.md Proposition 6(c)): turn "GB = {1} modulo two primes" into a PROOF
of emptiness over C by an exact Nullstellensatz certificate 1 = sum g_i f_i over Q.

  python3 lift_cert.py make PREFIX OUT.sing     # write a Singular script (system rebuilt over Q from the two primes)
  python3 lift_cert.py check PREFIX OUT.txt     # parse Singular's cofactors and check sum g_i f_i == 1 EXACTLY (Python)

Singular (on a Colab) only proposes the cofactors; the proof is the exact identity checked here in pure Python.
The system checked is the one rebuilt over Q from nag-rerun's mod-p exports (rational reconstruction mod p1, checked
mod p2); that rebuild is itself only as good as the export (see NOTE.md)."""
import sys, os, re
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zsys import *

def sing_poly(P, names):
    terms = []
    for k, c in P.items():
        mon = '*'.join(f"{names[i]}^{e}" if e > 1 else names[i] for i, e in enumerate(k) if e)
        cs = f"({c.numerator}/{c.denominator})" if Fraction(c).denominator != 1 else f"({int(c)})"
        terms.append(cs + ('*' + mon if mon else ''))
    return '+'.join(terms) if terms else '0'

def parse_sing(s, names):
    s = s.replace(' ', '').replace('\n', '')
    if s in ('0', ''):
        return {}
    s = s.replace('-', '+-')
    idx = {v: i for i, v in enumerate(names)}; P = {}
    for term in s.split('+'):
        if not term: continue
        c = Fraction(1); e = [0] * len(names)
        for f in term.split('*'):
            if f == '-': c = -c; continue
            neg = f.startswith('-')
            if neg: f = f[1:]; c = -c
            if re.fullmatch(r'\d+(/\d+)?', f):
                c *= Fraction(f)
            elif f == '':
                continue
            else:
                if '^' in f:
                    v, k = f.split('^'); e[idx[v]] += int(k)
                else:
                    e[idx[f]] += 1
        k = tuple(e); P[k] = P.get(k, 0) + c
    return {k: v for k, v in P.items() if v}

def mul(P, Q):
    out = {}
    for k1, c1 in P.items():
        for k2, c2 in Q.items():
            k = tuple(a + b for a, b in zip(k1, k2)); out[k] = out.get(k, 0) + c1 * c2
    return {k: v for k, v in out.items() if v}

if __name__ == '__main__':
    mode, pre, out = sys.argv[1:4]
    names, polys = load_Q(pre)
    # drop the appended 'extras' (alpha-only elements of the RELAXATION's GB mod p, appended last by fullsys.py): they are
    # sound for nag-rerun's screening but would make an exact certificate depend on the relaxation; we certify without them
    m = re.match(r'(.*)\.a(\d+)$', pre)
    if m and os.path.exists(m.group(1) + '.log') and not os.environ.get('KEEP_EXTRAS'):
        mm = re.search(rf"p1 array {m.group(2)}:.*extras (\d+)", open(m.group(1) + '.log').read())
        if mm and int(mm.group(1)):
            k = int(mm.group(1)); polys = polys[:-k]
            used = sorted({i for P in polys for kk in P for i, e in enumerate(kk) if e})
            names = [names[i] for i in used]; polys = [{tuple(kk[i] for i in used): c for kk, c in P.items()} for P in polys]
            print(f"dropped {k} appended extra row(s); {len(polys)} rows, {len(names)} unknowns remain")
    if mode == 'makem':
        # beta-free certificate (a2-nag's pattern shape): find y(alpha) with sum_r y_r A_r = 0 and sum_r y_r c_r = 1, i.e.
        # lift e = (0,...,0,1) through the rows (A_r1..A_rc, c_r) of the module over Q[alpha] (alphas only: small ring).
        # Computed mod p1, reconstructed and checked exactly like the others.  Incomplete in principle (a certificate
        # may need beta), complete enough in practice; 'make'/'makep' remain the fallback.
        ai, bi = split_linear(names, polys)
        an = [names[i] for i in ai]
        def apoly(P):
            return {tuple(k[i] for i in ai): c for k, c in P.items()}
        gens = []
        for P in polys:
            comps = []
            for b in bi:
                comps.append({k: c for k, c in P.items() if k[b] == 1})
            comps.append({k: c for k, c in P.items() if all(k[b] == 0 for b in bi)})
            gens.append('[' + ','.join(sing_poly(apoly(Q), an) if Q else '0' for Q in comps) + ']')
        e = '[' + ','.join(['0'] * len(bi) + ['1']) + ']'
        s = ["short = 0;", f"ring r = 1073741789, ({','.join(an) if an else 'zz'}), dp;", "option(redSB);",
             "module M = " + ",\n  ".join(gens) + ";", f"module E = {e};",
             "module G = std(M);", "if (size(reduce(E, G)) == 0) {", "  matrix T = lift(M, E);",
             f"  int i; for (i = 1; i <= {len(polys)}; i++) {{ print(\"COF \" + string(i)); print(string(T[i,1])); }}",
             "} else { print(\"NOT IN MODULE\"); }", "quit;"]
        open(out, 'w').write('\n'.join(s) + '\n')
        print('wrote', out)
    elif mode in ('make', 'makep'):
        ch = '0' if mode == 'make' else '1073741789'      # makep: cofactors computed mod p1, reconstructed + checked over Q
        s = ["short = 0;", f"ring r = {ch}, ({','.join(names)}), dp;", "option(redSB);",
             "ideal I = " + ",\n  ".join(sing_poly(P, names) for P in polys) + ";",
             "ideal G = std(I);", 'print("GB:"); print(G);',
             "if (size(G) == 1 && deg(G[1]) == 0) {", "  matrix T = lift(I, ideal(1));",
             f"  int i; for (i = 1; i <= {len(polys)}; i++) {{ print(\"COF \" + string(i)); print(string(T[i,1])); }}",
             "}", "quit;"]
        open(out, 'w').write('\n'.join(s) + '\n')
        print('wrote', out)
    else:
        txt = open(out).read()
        cofs = re.findall(r'COF (\d+)\n(.*?)(?=\nCOF |\Z)', txt, re.S)
        if not cofs:
            print('NO CERTIFICATE (Singular found no lift: GB over Q is not {1}?)'); sys.exit(1)
        tot = {}
        P1 = 1073741789
        modp = 'field characteristic' in txt or bool(re.search(r'\d{9,}', ''.join(c for _, c in cofs)))
        def recon(g):
            out = {}
            for k, v in g.items():
                r = ratrec(int(v) % P1, P1) if v.denominator == 1 else v
                if r is None: raise SystemExit('CERTIFICATE CHECK FAILED [rational reconstruction failed]')
                out[k] = r
            return out
        for i, s in cofs:
            g = parse_sing(s.strip(), names)
            if modp: g = recon(g)
            for k, v in mul(g, polys[int(i) - 1]).items():
                tot[k] = tot.get(k, 0) + v
        tot = {k: v for k, v in tot.items() if v}
        one = {tuple([0] * len(names)): Fraction(1)}
        maxdeg = max((sum(k) for _, s in cofs for k in parse_sing(s.strip(), names)), default=0)
        print(('PROVED EMPTY over C [exact certificate 1 = sum g_i f_i over Q verified; ' if tot == one else
               'CERTIFICATE CHECK FAILED [') + f"{len(cofs)} cofactors, max cofactor degree {maxdeg}]")
