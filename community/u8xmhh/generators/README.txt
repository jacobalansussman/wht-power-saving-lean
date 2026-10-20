generators/ of community/u8xmhh: the programs that chose the host table and wrote the Lean files
=================================================================================================

lean-gen/
  gen_hu8xm.py      writes the data modules (Work/XHost/Data/HU8XMSched.lean: host rows, bank rows, pads and
                    remainders, checked in Lean by hostSchedK; RateU8XMHHx.lean: the rate certificate;
                    InstU8XMHH.lean: the unit price) and the glue Work/XHost/U8XMHH.lean. It is adapted from
                    gen_h588m.py of community/c588mhh (generators/lean-gen/), for the certificate of
                    community/u8xm in the place of that of community/c588m.
  mk_seam_spec.py, mk_fourier_hx.py
                    write the Fourier chain Work/FourierU8XMHH/ and the two statement files: a renamed copy of the
                    Fourier files of community/u8xm, with the seam pointed at XHost.engine_of_ginv_host.
  mk_mac_test.py    writes the two Lean test files with which hostSchedK was first run on the Mac on the generated
                    tables (logs/generation/mac-FsPrice_HsU8XMHH*.log); the test files themselves are not included.
  The twelve network modules of Work/XHost (HostK to HPriceK in MODULES.txt) were not written for this package:
  they are byte-identical to those of community/c588mhh. Six of them (HRoles, HDef, HStage, HCert, HGeom, HEnd)
  are written by gen_hnet.py and gen_hgeom.py in community/c588mhh/generators/lean-gen/; the other six (HostK,
  Frames, HSched, HBank, HSchedK, HPriceK) have no generator in this repository (lean-gen/hsched_text.lean there
  is a copy of HSched.lean).

design-tools/
  host_price.py     the integer programme that places the leftover helper pieces in idle data registers of the
                    same run (same-run hosting) and puts the rest in banks
  host_repack.py, host_repack2.py   re-pack the pieces that are not hosted into fewer banks
  pad_sa.py         merges the padding of the banks (310 banks, pads 1 and 61)
  emit_hosttab.py   writes the host table of one run (the same table serves all 8 runs)
  audit_slots.py, indep_audit.py    independent checks of the placement and of the generated Lean tables
  controls.py, summary.py           the pricing controls (the figures of Shea's result, community/c588m and
                    community/u8xm are reproduced) and the table of all prices (data/PRICES.json)
  exact_moment.py, foldlib.py, gen_rates.py, gen_sharp.py, route_price.py
                    unchanged copies of files that are already elsewhere in this repository (see below)

scripts/
  run_u8xm.sh, run_samerun.sh, finish.sh   the order in which the programs above ran (scripts-as-run/
                    REPRODUCE-price.sh is the whole sequence as it ran)

Running them. The programs were run from the agents' working folder, whose sub-folders are named in the code:
  m-gen/toolkit/tools/gen    is this repository's tools/gen (foldlib.py, foldgen.py, gen_rates.py, gen_sharp.py
                             are byte-identical);
  x-pack/tools               holds exact_moment.py, foldlib.py and rate_gap.py, byte-identical to the files of the
                             same names in this repository's tools/community/repack/tools/;
  k-combine/out/global-U8XM/Work/GCert/Data/Gen/U8XMSched.lean
                             is community/u8xm/overlay/Work/GCert/Data/Gen/U8XMSched.lean (byte-identical);
  k-combine/out/price/U8XM.json  is data/U8XM-price.json of this package;
  k-combine/pkg/k-combine-u8xm/overlay  is community/u8xm/overlay with the notice of its two statement files
                             removed (tools/community/package_check.py same compares them that way).
Set the variable BC to a folder laid out like that. host_price.py also needs SciPy (1.13.1 was used) and
NumPy (2.0.2). Regenerating the Lean files from these copies was not re-run for this package.

Changes. Nine lines differ from the copies that ran: the absolute paths of the agents' working folder are read
from the variables BC and F, and one sys.path line of host_price.py is a comment. CHANGES.txt lists each of them
with the line as it ran. Nothing else was changed; the agent labels in the docstrings (fs-price, h-L3-design,
h-L3-build, s-stack, x-hosting) are the names the agents went by.

Copies. exact_moment.py is an unchanged copy of Ryan Shea's research/complex-transfer-588-596/verification/
exact_moment.py (NOTICE, section 9, item e); route_price.py copies functions of his and says so. foldlib.py,
gen_rates.py and gen_sharp.py are unchanged copies of this repository's tools/gen/. audit_slots.py, exact_moment.py,
foldlib.py, gen_rates.py, gen_sharp.py and route_price.py are byte-identical to the files of the same names in
community/c588mhh/generators/design-tools/.
