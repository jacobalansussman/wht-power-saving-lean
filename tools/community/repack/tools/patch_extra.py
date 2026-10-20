"""The three files of Shea's complex588 tree that depend on the bank count / exponent but are NOT written by
generate_global_lean.py: the Lean negative control (bank count) and the two challenge statements (positive: DFT
numerator; Neg = one-unit-stronger control: kernel numerator).  Usage:
patch_extra.py <shea lean/complex588> <out root> <banks> <kernel A> <dft numerator>"""
import sys, os
C, out, nb, A, d = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
jobs = [('Work/GCert/Data/C588Controls.lean', [('GX.schedK 5040 ', 'GX.schedK %s ' % nb)]),
        ('Work/FourierC588/UniformFourierChallenge.lean', [('1 - 12331685/(10^(10:ℕ))', '1 - %s/(10^(10:ℕ))' % d)]),
        ('Work/FourierC588Neg/UniformFourierChallenge.lean', [('1 - 12331686/(10^(10:ℕ))', '1 - %s/(10^(10:ℕ))' % A)])]
for rel, subs in jobs:
    t = open(os.path.join(C, rel)).read()
    for a, b in subs:
        assert t.count(a) == 1, (rel, a)
        t = t.replace(a, b)
    os.makedirs(os.path.dirname(os.path.join(out, rel)), exist_ok=True)
    open(os.path.join(out, rel), 'w').write(t)
    print('patched', rel)
