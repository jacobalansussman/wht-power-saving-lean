import Work.FourierC588.Main
#print axioms OAI.PowerSaving.transform_mainC588
#print axioms OAI.PowerSaving.convolution_mainC588
#print axioms OAI.PowerSaving.RAM.hillsC588_program
#print OAI.PowerSaving.decimalExponentC588
#print OAI.PowerSaving.TimeBoundsC588
#print OAI.PowerSaving.DFTGoalC588
#print OAI.PowerSaving.ConvGoalC588
#check @OAI.PowerSaving.transform_mainC588
#check @OAI.PowerSaving.convolution_mainC588
#check @OAI.PowerSaving.RAM.hillsC588_program
example : OAI.PowerSaving.decimalExponentC588 = 1 - 12383291/(10^(10:ℕ)) := rfl
