LoadPackage("ctbllib");
t := CharacterTable("S8(2)");
Print("S8(2) classes ", NrConjugacyClasses(t), "\n");
Print("degrees ", List(Irr(t), x -> x[1]), "\n");
Print("maxes ", Maxes(t), "\n");
Print("classnames ", ClassNames(t), "\n");
ia := Position(ClassNames(t), "2a");
Print("chi(2a) ", List(Irr(t), x -> x[ia]), "\n");
for nm in Maxes(t) do
  s := CharacterTable(nm);
  if Size(t)/Size(s) in [120, 136, 255, 2295] then
    pc := InducedClassFunction(TrivialCharacter(s), t);
    Print("PERM ", nm, " index ", Size(t)/Size(s), " : ", MatScalarProducts(t, Irr(t), [pc])[1], "\n");
  fi;
  if Size(t)/Size(s) = 2295 then
    Print("P0 ", nm, " degrees ", List(Irr(s), x -> x[1]), "\n");
    dec := MatScalarProducts(s, Irr(s), List(Irr(t), chi -> RestrictedClassFunction(chi, s)));
    for i in [1..Length(dec)] do
      if Irr(t)[i][1] <= 2000 then Print("RES ", Irr(t)[i][1], " ", Irr(t)[i][ia], " : ", dec[i], "\n"); fi;
    od;
  fi;
od;
QUIT;
