"""vall.py -- run vcert.py on every certificate output of ex2-modp, mapping each to its RAW export.
Also checks that rows dropped as 'extras' are alpha-only (so the dropped rows really are the appended extras)."""
import os, re, glob, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import read_ms
M = os.environ.get('MODP_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
HERE = os.path.dirname(os.path.abspath(__file__))
jobs = []
for f in sorted(glob.glob(f'{M}/cert*/*.txt')):
    d, b = os.path.basename(os.path.dirname(f)), os.path.basename(f)[:-4]
    keep = False
    if d == 'cert':
        if b.startswith('test'): continue
        m = re.match(r't27F5_(.*)', b) or None
        if b.startswith('t27F5_'): pre = f'{M}/data/dec/t27_Fmid3_m3_s5/full_c_{b[6:]}'
        elif b.startswith('t27s1_scr_'): pre = f'{M}/data/dec/t27_mid3_m2_s1/{b[10:]}'
        elif b.startswith('t27s1_'): pre = f'{M}/data/dec/t27_mid3_m2_s1/full_c_{b[6:]}'
        else: pre = None
    elif d in ('cert3', 'certm', 'cert2'):
        sk, sy = b.split('__', 1); pre = f'{M}/data/dec/{sk}/{sy}'
    elif d == 'certh':
        sk, sy = b.split('__', 1); pre = f'{M}/data/hard/dec/{sk}/{sy}'
    elif d == 'certh2':
        if '__SUB__' in b:
            sk, sy = b.split('__SUB__'); pre = f'{M}/data/hard/sub_{sk}/{sy}'
        else:
            sk, sy = b.split('__', 1); pre = f'{M}/data/hard/dec/{sk}/{sy}'
    elif d == 'certe':
        if b.startswith('ext_'): pre = f'{M}/data/hard/ext_s17/{b}'
        elif b.startswith('withext_'): pre = f'{M}/data/hard/dec/g27_Fany4_s17/full_c_{b[8:]}'; keep = True
        else: pre = None
    elif d == 'certp':
        pre = f'{M}/data/hard/dec/g27_Fany4_s17/full_c_{b[8:]}'; keep = True
    elif d == 'certc':
        pre = f'{M}/data/hard/comb/{b}'
    elif d == 'cert5':
        if b.startswith('SUB__'): pre = f'{M}/data/sub_f27_s39/{b[5:]}'
        else:
            sk, sy = b.split('__', 1); pre = f'{M}/data/dec/{sk}/{sy}'
    else:
        pre = None
    jobs.append((d, b, pre, keep, f))
for d, b, pre, keep, f in jobs:
    if pre is None or not os.path.exists(pre + '.p1.ms'):
        print(f'UNMAPPED {d}/{b} -> {pre}'); continue
    if os.path.getsize(f) == 0:
        print(f'EMPTY-OUTPUT {d}/{b}'); continue
    r = subprocess.run(['python3', f'{HERE}/vcert.py', pre, f] + (['KEEP_EXTRAS'] if keep else []),
                       capture_output=True, text=True, timeout=600)
    line = (r.stdout.strip().split('\n') or [''])[-1] if r.stdout.strip() else 'ERROR ' + r.stderr.strip().split('\n')[-1]
    # extras sanity: dropped rows alpha-only?
    note = ''
    m = re.search(r'drop=(\d+)', line)
    if m and int(m.group(1)):
        k = int(m.group(1)); names, p, P = read_ms(pre + '.p1.ms')
        bad = [i for i in range(len(P) - k, len(P)) if any(any(e and not names[j].startswith('a') for j, e in enumerate(mon)) for mon in P[i])]
        note = ' dropped-rows-alpha-only=' + ('YES' if not bad else f'NO {bad}')
    print(f'{d}/{b}: {line.split("  /")[0]}{note}', flush=True)
