import Work.XHostT.KWhole
import Work.XHostT.Data.H2979.Order
import Work.XHostT.Data.H2979.P0
import Work.XHostT.Data.H2979.P1
import Work.XHostT.Data.H2979.P2
import Work.XHostT.Data.H2979.P3
import Work.XHostT.Data.H2979.P4
import Work.XHostT.Data.H2979.P5
import Work.XHostT.Data.H2979.P6
import Work.XHostT.Data.H2979.P7
import Work.XHostT.Data.H2979.P8
import Work.XHostT.Data.H2979.P9
import Work.XHostT.Data.H2979.P10
import Work.XHostT.Data.H2979.P11
import Work.XHostT.Data.H2979.P12
import Work.XHostT.Data.H2979.P13
import Work.XHostT.Data.H2979.P14
import Work.XHostT.Data.H2979.P15
import Work.XHostT.Data.H2979.P16
import Work.XHostT.Data.H2979.P17
import Work.XHostT.Data.H2979.P18
import Work.XHostT.Data.H2979.P19
import Work.XHostT.Data.H2979.P20
import Work.XHostT.Data.H2979.P21
import Work.XHostT.Data.H2979.P22
import Work.XHostT.Data.H2979.P23
import Work.XHostT.Data.H2979.P24
import Work.XHostT.Data.H2979.P25
import Work.XHostT.Data.H2979.P26
/-! p-hostck mk_parts.py: the whole host check of H2979 (H2979-hosttab.json.gz, sha256 4c6720e49ac262a2...): 27 parts, 666 host kinds, 666 calls x 408 slots = 271728 records. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979
open OAI.PowerSaving.XHostT.K
/-- every host kind of the table passes the kernel check (assembled from the 27 parts). -/
theorem all_kinds : KOK ct rt resT NC R cK 0 NK :=
  (KOK.app (a := 0) (b := 325) (c := 666) (KOK.app (a := 0) (b := 150) (c := 325) (KOK.app (a := 0) (b := 75) (c := 150) (KOK.app (a := 0) (b := 25) (c := 75) rng0 (KOK.app (a := 25) (b := 50) (c := 75) rng1 rng2)) (KOK.app (a := 75) (b := 100) (c := 150) rng3 (KOK.app (a := 100) (b := 125) (c := 150) rng4 rng5))) (KOK.app (a := 150) (b := 225) (c := 325) (KOK.app (a := 150) (b := 175) (c := 225) rng6 (KOK.app (a := 175) (b := 200) (c := 225) rng7 rng8)) (KOK.app (a := 225) (b := 275) (c := 325) (KOK.app (a := 225) (b := 250) (c := 275) rng9 rng10) (KOK.app (a := 275) (b := 300) (c := 325) rng11 rng12)))) (KOK.app (a := 325) (b := 500) (c := 666) (KOK.app (a := 325) (b := 400) (c := 500) (KOK.app (a := 325) (b := 350) (c := 400) rng13 (KOK.app (a := 350) (b := 375) (c := 400) rng14 rng15)) (KOK.app (a := 400) (b := 450) (c := 500) (KOK.app (a := 400) (b := 425) (c := 450) rng16 rng17) (KOK.app (a := 450) (b := 475) (c := 500) rng18 rng19))) (KOK.app (a := 500) (b := 575) (c := 666) (KOK.app (a := 500) (b := 525) (c := 575) rng20 (KOK.app (a := 525) (b := 550) (c := 575) rng21 rng22)) (KOK.app (a := 575) (b := 625) (c := 666) (KOK.app (a := 575) (b := 600) (c := 625) rng23 rng24) (KOK.app (a := 625) (b := 650) (c := 666) rng25 rng26)))))
/-- global fact: the record counts sum to `NC * R`. -/
theorem sum_ok : ((List.range NK).map cK).sum = NC * R := by decide +kernel
/-- **the whole table**: a valid cross-replica host schedule (TableOK) with a valid global order. -/
theorem whole : TableOK ct rt resT NC R NK ∧ chkOrder ct NC NC 0 0 = true :=
  ⟨tableOK_of_KOK cK all_kinds sum_ok, order_ok⟩
end OAI.PowerSaving.XHostT.K.DH2979
#print axioms OAI.PowerSaving.XHostT.K.DH2979.whole
