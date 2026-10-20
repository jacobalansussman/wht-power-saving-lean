LoadPackage("ctbllib");
d := CharacterTable("2.S6(2)");
Print("maxes ", Maxes(d), "\n");
cen := ClassPositionsOfCentre(d);
faith := Filtered(Irr(d), chi -> chi[cen[2]] = -chi[1]);
Print("faithful degrees ", List(faith, x -> x[1]), "\n");
mx := Maxes(d);
nm := First(mx, x -> PositionSublist(x, "L3(2)") <> fail);
Print("preimage of P0: ", nm, "\n");
s := CharacterTable(nm);
Print("its order ", Size(s), "\n");
for chi in faith do
  r := RestrictedClassFunction(chi, s);
  dec := MatScalarProducts(s, Irr(s), [r])[1];
  Print("FAITH ", chi[1], " : max mult ", Maximum(dec), " pieces ", Sum(dec), " dims ", Concatenation(List([1..Length(dec)], i -> ListWithIdenticalEntries(dec[i], Irr(s)[i][1]))), "\n");
od;
QUIT;
