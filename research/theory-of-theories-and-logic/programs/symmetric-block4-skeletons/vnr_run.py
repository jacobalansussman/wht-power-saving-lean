"""vnr_run.py -- v-nag-rerun driver (runs on a Colab with Singular).
  relax NAME           : relaxation with every coefficient free -> Singular std over Q -> HIT (GB != 1) / EMPTY (GB = 1)
  fz NAME              : forced zeros: for every alpha a_i, std(I_relax + <t a_i - 1>) = 1 ?
  full NAME TMO        : full closure system of array 0 with the forced zeros (from out/NAME.fz) dropped, every other
                         alpha FREE (no charts), constant-pivot elimination, Singular std over Q.
Outputs in out/."""
import sys, os, json, time, subprocess
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import vnr
REC = os.path.join(HERE, 'recipes')
OUT = os.path.join(HERE, 'out'); os.makedirs(OUT, exist_ok=True)

def sing(text, path, tmo):
    open(path, 'w').write(text)
    t0 = time.time()
    try:
        r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=tmo)
        out = r.stdout
    except subprocess.TimeoutExpired:
        return 'TIMEOUT', time.time() - t0
    for l in out.splitlines():
        if l.startswith('RESULT'):
            return ('GB=1' if 'GB=1' in l and 'GB!=1' not in l else 'GB!=1'), time.time() - t0
    return 'ERR ' + out[-300:], time.time() - t0

def main():
    mode, name = sys.argv[1], sys.argv[2]
    rec = json.load(open(os.path.join(REC, name + '.json')))
    if mode == 'relax':
        sk = vnr.Skel(rec)
        txt, av, bv = vnr.relax_singular(sk, name)
        res, dt = sing(txt, os.path.join(OUT, name + '.relax.sing'), 600)
        print('RELAX', name, 'W', rec['W'], 'cost', sk.cost, 'gens', len(sk.gens), 'alphas', len(av), res, round(dt, 1), flush=True)
    elif mode == 'fz':
        sk = vnr.Skel(rec)
        fz = []
        for a in range(sk.na):
            txt, av, bv = vnr.relax_singular(sk, name)
            if a not in av:
                print('FZ', name, a, 'absent from relaxation: kept free', flush=True); continue
            txt = txt.replace('ring r = 0,(', 'ring r = 0,(tt,').replace(';\nideal G', f',\ntt*a{a}-1;\nideal G', 1)
            txt = txt.replace('ideal I =\n', 'ideal I =\n', 1)
            # insert the extra generator before the closing ';' of ideal I
            res, dt = sing(txt, os.path.join(OUT, f'{name}.fz{a}.sing'), 600)
            print('FZ', name, a, res, round(dt, 1), flush=True)
            if res == 'GB=1': fz.append(a)
            elif res != 'GB!=1': print('FZ-UNDECIDED', name, a, res, flush=True)
        json.dump(fz, open(os.path.join(OUT, name + '.fz'), 'w'))
        print('FORCED', name, fz, flush=True)
    elif mode == 'full':
        tmo = int(sys.argv[3])
        fz = json.load(open(os.path.join(OUT, name + '.fz')))
        sk = vnr.Skel(rec, zero=set(fz))
        t0 = time.time()
        cl = sk.closure(0, 60000)
        if cl is None:
            print('FULL', name, 'closure > cap', flush=True); return
        A, nc = cl
        core = vnr.const_eliminate(A, nc)
        contr = any(set(R) == {'r'} and set(R['r']) == {()} for R in core)
        nb = len({c for R in core for c in R if c != 'r'})
        print('FULL', name, 'forced', fz, 'closure', len(A), 'x', nc, 'core', len(core), 'x', nb, 'contradiction' if contr else '',
              round(time.time() - t0, 1), flush=True)
        if contr:
            print('FULLRES', name, 'GB=1 (constant contradiction)', flush=True); return
        txt, av, bv = vnr.core_singular(core, name)
        res, dt = sing(txt, os.path.join(OUT, name + '.full.sing'), tmo)
        print('FULLRES', name, res, 'alphas', len(av), 'betas', len(bv), round(dt, 1), flush=True)

if __name__ == "__main__":
    main()
