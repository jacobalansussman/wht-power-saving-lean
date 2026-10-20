# ctbl.g -- character-table facts for Sp(6,2) = S6(2) and the double cover, restrictions to P_0 = 2^6:L3(2).
LoadPackage("ctbllib");
t := CharacterTable("S6(2)");
Print("S6(2) order ", Size(t), " classes ", NrConjugacyClasses(t), "\n");
Print("class names ", ClassNames(t), "\n");
Print("degrees ", List(Irr(t), x -> x[1]), "\n");
Print("maxes ", Maxes(t), "\n");
nm := First(Maxes(t), x -> x = "2^6:L3(2)");
Print("P0 table name ", nm, "\n");
s := CharacterTable(nm);
Print("P0 order ", Size(s), " irr degrees ", List(Irr(s), x -> x[1]), "\n");
Print("P0 class names ", ClassNames(s), "\n");
Print("P0 class orders ", OrdersClassRepresentatives(s), "\n");
fus := GetFusionMap(s, t);
Print("fusion P0->G ", fus, "\n");
rest := List(Irr(t), chi -> RestrictedClassFunction(chi, s));
dec := MatScalarProducts(s, Irr(s), rest);
Print("RESTRICTION MATRIX (rows: G-irreps by degree order, cols: P0-irreps)\n");
for i in [1..Length(dec)] do Print(Irr(t)[i][1], " : ", dec[i], "  sumsq=", Sum(dec[i], x->x^2), "\n"); od;
# character values on classes for transvection etc.
Print("Irr values table (rows irreps):\n");
for chi in Irr(t) do Print(List(chi), "\n"); od;
Print("centralizer orders ", SizesCentralizers(t), "\n");
Print("element orders ", OrdersClassRepresentatives(t), "\n");
# P0-irreps values
Print("P0 Irr values:\n");
for chi in Irr(s) do Print(List(chi), "\n"); od;
Print("P0 centralizers ", SizesCentralizers(s), "\n");
# double cover
d := CharacterTable("2.S6(2)");
if d <> fail then
  Print("2.S6(2) degrees ", List(Irr(d), x -> x[1]), "\n");
  Print("2.S6(2) maxes ", Maxes(d), "\n");
fi;
QUIT;
