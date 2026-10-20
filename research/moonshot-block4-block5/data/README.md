# Data of this section

The data files sit next to the programs that read them, so that every command in the program READMEs runs as given.
This page lists them, with their sha256, and names the large files that are not published and how to rebuild them.

## Published data

| file | what it is | sha256 |
|---|---|---|
| `../programs/two-in-seven-prefix-end-tests/prefixes/all_prefixes.jsonl` | the 1,261 saved five-move prefixes of the 2-in-7 searches (one JSON line each) | `cd1ad4ede41706dd87b5f9bc898cd843924b070242c47413c3f4bbea9f084cd5` |
| `../programs/seven-in-23-gl14-rank-bounds/dec_s21_gl14.npz` | the decomposition of the slot module in gl14 (lift seed 21) used by `pair.py` and `one65.py` | `06a7c73c8ef22a736e89bad82af663a7046a95b9a193956330930e99ed706539` |
| `../programs/seven-in-23-gl14-rank-bounds/dec_s21_V.npz` | the same in the natural 16-dimensional representation V | `461b3f93d456df5b2d62f42b916dc444e05eaac1491a0bfe903453b449301bb7` |
| `../programs/one-shared-move-exact-decisions/sing*/` | the Singular inputs of the one-shared-move exact decisions (with `out_singular*.txt`, the verdicts) | per file; regenerate and compare with `exact_small.py` |
| `../programs/block5-ball-test/systems/m801_k2_d0.sing` | s801, radius-2 ball, modulo p = 1,073,741,789 | `8e7e907b887297a6f3443249533f34640507c5dafc0a4b6ea5c45dcb567e6299` |
| `../programs/block5-ball-test/systems/m801_k2_d1.sing` | the same with the layer-4 rows | `d97dd246cfec1bf4fc122fabc72fdcce51b3abeb2876635da98f43cdf75935cd` |
| `../programs/block5-ball-test/systems/m1401_k2_d0.sing` | s1401, radius-2 ball | `56b547888b6efd124f3fbfc3a1abfc204d1c830c3b6a07700031205c95ad2074` |
| `../programs/block5-ball-test/systems/m1401_k2_d1.sing` | the same with the layer-4 rows | `855285a15faf14d78e1c76b899adbe64a3298ca1987164413d410e3d9bdff45b` |
| `../programs/block5-ball-test/systems/m1402_k2_d0.sing` | s1402, radius-2 ball | `33697da74dbb408481cc4f004bfed4bbd083081a08b3eab226202bbc86d49c77` |
| `../programs/block5-ball-test/systems/m1402_k2_d1.sing` | the same with the layer-4 rows | `65f67bb102884c6895850cdfcdf526b71c59499fe862a66c39511d05cc499bed` |
| `../programs/block5-ball-test/systems/q*_k2_*.sing` | the same six systems over the rationals | see `shasum -a 256` |
| `../programs/block5-ball-test/systems/m*_k2_*.ms` | the msolve inputs of the six modular systems | see `shasum -a 256` |
| `../programs/block5-ball-test/systems/m1501_k2_d1.sing` | the open design s1501, radius 2, with the layer-4 rows | `926d147d4cc5ddc94c92d44c900594d1a41935546efc463a7be19609c9874272` |
| `../programs/block5-ball-test/systems/m1501_k3_d1.sing` | s1501, radius 3, with the layer-4 rows (1,220 unknowns, 2,154 equations) | `af00e1ba2fbc9ae00c5042c805a7cf28288a0f8241a05e8dcb6db28ba81dbeee` |
| `../programs/block5-ball-test/systems/q1501_k2_d1.sing`, `q1501_k3_d1.sing` | the same two over the rationals | `ccd394612be5492ccf29b755b98cdbc0d0c49aa887ee761e7ea2527990ea30e2`, `c00b52732d91b2dbf33ad20da15528b7c32665297e2c55ada99f0c4439203658` |

## Not published, with how to rebuild them

- **The layered block-5 systems** (Python pickles, 1 KB to 4.4 MB each). They are not published because loading a
  pickle file can run code. They are rebuilt by `../programs/block5-sparse-designs/dumpsys.py`; the commands that made
  them are in `t34.sh` (s1501), `t35.sh`, `t42_1401.sh` and `t42_1402.sh` (s1401, s1402) and `t36.sh` (s801) in that
  folder. A rebuild need not be byte-identical (pickles of dictionaries); what matters is that `vball.py` then writes the
  same systems as above. The files that were used had these sha256:

  | file | sha256 |
  |---|---|
  | `sys_d3_801.pkl` | `a27c10ab71b3190b5f5a3680bf4f9914eb91a85d8f5cf6f7c06c3338201750dc` |
  | `sys_d4_801.pkl` | `7cb4984e2ccbac8c61dedfb38fdbd3c601c2d716f945da3401cedecc77bf3db2` |
  | `sys_d3_1401.pkl` | `a18484e6c10b6da442bb2e70607c33ac587ef4b259fac7ea4b53ce759dcf6fdb` |
  | `sys_d4_1401.pkl` | `59dedaf635158ecf501e3485257bad8ad106599f8ed95889e71027bf26192d22` |
  | `sys_d3_1402.pkl` | `d274e064bae371d98d5c14f8d27aa419d4933d8fe1418ba99ec81805280e061a` |
  | `sys_d4_1402.pkl` | `47cf47eeeb5ec35d7ab24597bec41ada54825504eb96fb96a90416eabe7b50cc` |
  | `sys_d3_1501.pkl` | `790030a225e6a5e492ba218d84e2c8f6b412542cdef5ac1c0738e5e80245e582` |
  | `sys_d4_1501.pkl` | `ee2d401231950d4e45445ef385c4eb9eb74afd7055aacc0a690d7f13b8be6059` |

- **Search outputs and solver logs** of the order-3 round (79 MB in all), of the gl14 track's run folders (63 MB) and of
  the last block-5 round and its skeptic (17 MB and 58 MB). They record searches that found nothing and solver runs
  whose verdicts are summarised in the program folders. Each program README gives the command that regenerates its
  part.
