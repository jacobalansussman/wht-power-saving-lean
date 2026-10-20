"""Copy the Fourier files of C588M (community/c588m) (Shea's renaming of Mamano's record files) to a new suffix NEW, with
THE SEAM calling XHost.engine_of_ginv_host on a host schedule.  Certificate-level names (GXD.C588M.*, BridgeRate.C588M,
BlockAccounting.*C588Mx) are kept unless overridden.
usage: python3 -I mk_fourier_h.py SRC_OVERLAY OUT NEW SEAMSPEC.json"""
import sys, os, re, json, hashlib
src, out, NEW, spec = sys.argv[1], sys.argv[2], sys.argv[3], json.load(open(sys.argv[4]))
OLD = 'C588M'
prov = {}
for sub in ('Work/FourierC588M', 'Work/FourierC588MNeg'):
    d = os.path.join(src, sub)
    dsub = sub.replace('FourierC588M', 'Fourier' + NEW)
    os.makedirs(os.path.join(out, dsub), exist_ok=True)
    for fn in sorted(os.listdir(d)):
        t = open(os.path.join(d, fn)).read()
        prov[sub + '/' + fn] = hashlib.sha256(t.encode()).hexdigest()
        t = t.replace(OLD, NEW)
        # certificate-level names stay those of the banked certificate (unless the spec says otherwise)
        for a, b in spec.get('keep', {}).items():
            t = t.replace(a, b)
        for a, b in spec.get('num', {}).items():
            t = t.replace(a, b)
        if fn == 'Seam.lean' and sub.endswith('C588M'):
            # imports
            t = re.sub(r'import Work\.GCert\.Data\.\S+\n', '', t)
            t = t.replace('import Work.Bridge401.End\n', 'import %s\n' % spec['import'])
            m = re.search(r'  GX\.engine_of_ginv_401 [\s\S]*?cl m k v\n', t)
            assert m, 'engine call not found'
            t = t[:m.start()] + spec['call'] + t[m.end():]
            t = t.replace('that calls `GX.engine_of_ginv_401`', 'that calls `XHost.engine_of_ginv_host` (h-L3-build: the HOSTED network)')
        open(os.path.join(out, dsub, fn), 'w').write(t)
for fn in ('fourier-c588m.json', 'fourier-c588m-neg.json'):
    t = open(os.path.join(src, 'comparator', fn)).read().replace(OLD, NEW)
    os.makedirs(os.path.join(out, 'comparator'), exist_ok=True)
    open(os.path.join(out, 'comparator', fn.replace('c588m', NEW.lower())), 'w').write(t)
json.dump(prov, open(os.path.join(out, 'fourier-%s-provenance.json' % NEW), 'w'), indent=1)
print('ok', NEW)
