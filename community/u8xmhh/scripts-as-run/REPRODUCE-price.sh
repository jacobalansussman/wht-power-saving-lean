#!/bin/zsh
# fs-price (wht43 bound-check, 13 Oct 2026): price U8XM + same-run hosting and produce the generator inputs and Lean texts
# of GXD.U8XMHH.  Python (python3 -I) except step 6 (two light Lean modules on the Mac through locked.sh, about 2.8 GB, 1-2 min).
set -e
F=<verify>/checks/wht43/bound-check/fs-price
H=~/misc/fft/verify/checks/wht42/common/heavy.sh
cd $F
# 1. controls of the pricing pipeline: Shea 12331685, C588M 13552606, U8XM 13868856 (seconds)
python3 -I tools/controls.py data/words.json out/controls-price.json
# 2. same-run hosting: C588MHH control (15699305, tables byte-identical to h-L3-design/h-L3-build), U8 no-copies control
#    (16254943), and the target U8XM_5cc (509 widths of k-combine's U8XM)  [~6 min each for the re-pack]
$H ./run_samerun.sh            # (its third item repeats run_u8xm.sh: a determinism check, identical output)
$H ./run_u8xm.sh hu8xm_sr
# 3. pad consolidation (concave pad price; 310 banks = width lower bound)
python3 -I tools/pad_sa.py out/hu8xm_sr.repacked2.json out/hu8xm_sr.padsa.json 4000000 11
python3 -I tools/audit_slots.py data/words.json U8XM_5cc out/hu8xm_sr.padsa.json
# 4. host table, Lean data/rate/instance/glue, Fourier files + seam on XHost.engine_of_ginv_host, comparator configs
./finish.sh hu8xm_sr out/hu8xm_sr.padsa.json
K=$F:h/k-combine/pkg/k-combine-u8xm/overlay/Work/GCert/Data
python3 -I tools/indep_audit.py lean/Work/XHost/Data/HU8XMSched.lean $K/Gen/U8XMSched.lean $K/InstU8XM.lean
# 5. cross-run reference (not expressible in tonight's Lean format)  [~10-25 min each]
$H ./run_cross.sh              # defer 1 (600 s, time limit) ; its defer 0 leg needs 1500 s:
$H ./run_cross_d0.sh
python3 -I tools/summary.py
# 6. Lean on the Mac: hostSchedK on the generated tables (true) and on a one-off copy (false)
python3 -I tools/mk_mac_test.py lean/Work/XHost/Data/HU8XMSched.lean $K/Gen/U8XMSched.lean mac/FsPrice
mac/mk.sh FsPrice/HsU8XMHH.lean 900 ; mac/mk.sh FsPrice/HsU8XMHHNeg.lean 900
