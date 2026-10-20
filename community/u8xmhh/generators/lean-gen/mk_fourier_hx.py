"""fs-price (from h-L3-build's mk_fourier_h.py, OLD made a parameter): copy the Fourier files of OLD (C588M: m-gen; U8XM: k-combine) (Shea's renaming of Mamano's record files) to a new suffix NEW, with
THE SEAM calling XHost.engine_of_ginv_host on a host schedule.  Certificate-level names (GXD.C588M.*, BridgeRate.C588M,
BlockAccounting.*C588Mx) are kept unless overridden.
usage: python3 -I mk_fourier_hx.py SRC_OVERLAY OUT OLD NEW SEAMSPEC.json"""
import sys, os, re, json, hashlib
src, out, OLD, NEW, spec = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], json.load(open(sys.argv[5]))
prov = {}
for sub in ('Work/Fourier' + OLD, 'Work/Fourier' + OLD + 'Neg'):
    d = os.path.join(src, sub)
    dsub = sub.replace('Fourier' + OLD, 'Fourier' + NEW)
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
        if fn == 'Seam.lean' and sub.endswith(OLD):
            # imports
            t = re.sub(r'import Work\.GCert\.Data\.\S+\n', '', t)
            t = t.replace('import Work.Bridge401.End\n', 'import %s\n' % spec['import'])
            m = re.search(r'  GX\.engine_of_ginv_401 [\s\S]*?cl m k v\n', t)
            assert m, 'engine call not found'
            t = t[:m.start()] + spec['call'] + t[m.end():]
            t = t.replace('that calls `GX.engine_of_ginv_401`', 'that calls `XHost.engine_of_ginv_host` (h-L3-build: the HOSTED network)')
        open(os.path.join(out, dsub, fn), 'w').write(t)
for fn in ('fourier-%s.json' % OLD.lower(), 'fourier-%s-neg.json' % OLD.lower()):
    t = open(os.path.join(src, 'comparator', fn)).read().replace(OLD, NEW)
    os.makedirs(os.path.join(out, 'comparator'), exist_ok=True)
    open(os.path.join(out, 'comparator', fn.replace(OLD.lower(), NEW.lower())), 'w').write(t)
json.dump(prov, open(os.path.join(out, 'fourier-%s-provenance.json' % NEW), 'w'), indent=1)
print('ok', NEW)
