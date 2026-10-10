import Work.FourierE8.Main

/-! Axiom report for the every-length Fourier theorems at our exponent (nothing is proved here). -/

#print axioms OAI.PowerSaving.transform_mainE8
#print axioms OAI.PowerSaving.convolution_mainE8
#print axioms OAI.PowerSaving.transform_main_orderE8
#print axioms OAI.PowerSaving.convolution_main_orderE8
#print axioms OAI.PowerSaving.RAM.hillsE8_program
#print OAI.PowerSaving.decimalExponentE8
#print OAI.PowerSaving.alphaE8
#check (OAI.PowerSaving.transform_mainE8 : OAI.PowerSaving.DFTGoalE8)
#check (OAI.PowerSaving.convolution_mainE8 : OAI.PowerSaving.ConvGoalE8)
