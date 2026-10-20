# Data of this section

| file | size | sha256 | what it is |
|---|---|---|---|
| `finite-fields-systems-and-certificates.tar.gz` | 4,488,787 bytes | `d252f97824329ef5e18dd9354cd07e03c63a6608aff53d5c224187c7c7c558e5` | 2,293 text files: the polynomial systems of the block-4 skeletons and relaxations exported modulo two primes (`data/dec/`, `data/hard/`, `data/ms/`, `data/msctrl/`, `data/sub_f27_s39/`, `data/toy/`), and the impossibility certificates (`cert*/`: Singular's proposed cofactors `.txt`, its scripts `.sing`, logs) |

The other small data of this section sit beside the programs that read them: the exact certificates of Theorems 17,
18 and E (`../programs/families-commutative-slice/rc_cert_b4.json`, `rc_cert_b5.json`, `rc_vars_b5.json`;
`../programs/star-lemma-and-slice/theorem-E-certificate/`), the stand-in schemes
(`../programs/certified-numerics/sx/out/scheme_2in3.json`, `scheme_3in5.json`), the exported relaxations used by the
calibration (`../programs/certified-numerics/sys/`) and the saved numerical solutions (`*.npz`, plain arrays, no
pickled objects).

**How to use it.** Unpack it inside `../programs/finite-fields/`, so that `data/` and `cert*/` sit beside `lift_cert.py`:

```
cd ../programs/finite-fields && tar xzf ../../data/finite-fields-systems-and-certificates.tar.gz
cd skeptic && python3 vall.py      # every certificate: the identity 1 = sum g_i f_i, exactly over Q and modulo p2
python3 vcover.py                  # does every exclusion step of the 8 skeletons have an accepted certificate?
python3 vpoints.py                 # every claimed existence point is an exact zero of its system
```

`vall.py` takes under a minute to about 3 minutes (pure Python, under 200 MB). Expected: 240 lines
`OK   Q-identity=True p2-raw-identity=True`, 13 lines `FAIL no cofactors` and 3 lines `EMPTY-OUTPUT`, which are the
documented failed attempts (timeouts and "not in module" answers of Singular, and runs that wrote nothing, replaced by
later runs); `vcover.py` ends with `TOTAL 121/122 steps have a
certificate accepted by my independent checker` (the last step is a system that reduces to the constant −1). The
programs find the folders through `MODP_DIR` (default: the parent of `skeptic/`).

**What is left out, and why.** msolve's own outputs (`*.gb`), timing files (`*.tm`) and the archives written before the
msolve 0.6.5 notice (`dec_neg.tgz`, `dec_more.tgz`, `dec_more2.tgz`, `dec_hard.tgz`, `msctrl.tgz`, `nagrr.tgz`). msolve
0.6.5 reads a product written `x*x` as `x`; 72 of the 112 full-chart exports contain such products. No result here
rests on an msolve verdict: every certificate and point is checked on polynomials read with `x*x` = `x²`, and Singular
re-read 15 raw exports with the same result. Never feed these exports to msolve 0.6.5 without rewriting powers as
`x^k`.

**What the certificates assume.** That the exported systems are the right ones (the closure model of the symmetric
block-4 skeletons, the relaxation lemma and the split into charts of the rounds told in
`../../theory-of-theories-and-logic/`). The skeptic rebuilt the systems exactly from their recipes and found 160 of 162
full-chart cores identical; the two others were certified separately. That rebuild needs the skeleton builders of
those rounds, which are not part of this section.
