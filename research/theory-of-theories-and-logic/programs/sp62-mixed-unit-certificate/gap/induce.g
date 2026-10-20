LoadPackage("ctbllib");
t := CharacterTable("S6(2)");
for nm in Maxes(t) do
  s := CharacterTable(nm);
  lin := Filtered(Irr(s), x -> x[1] = 1 and ForAll(x, v -> v = 1 or v = -1));
  for k in [1..Length(lin)] do
    ind := InducedClassFunction(lin[k], t);
    dec := MatScalarProducts(t, Irr(t), [ind])[1];
    Print("IND ", nm, " index ", Size(t)/Size(s), " lin#", k, " : ", dec, "\n");
  od;
od;
QUIT;
