Read("sp6_perms.g");
G := Group(gens);; P0 := Stabilizer(G, L0);;
T := CharacterTable(G);; irr := Irr(T);;
cl := ConjugacyClasses(T);;
pos := g -> PositionProperty(cl, c -> g in c);
ia := pos(hs[1]);; ib := pos(hs[1]*hs[2]^0);;
lab := List(irr, chi -> [chi[1], chi[ia]]);
Print("labels (deg, chi(h)) ", lab, "\n");
Print("class of h, h2, k: ", List(hs, pos), "  class sizes ", List(hs, g -> Size(cl[pos(g)])), "\n");
B := SylowSubgroup(P0, 2);; Print("B order ", Size(B), " abelianization ", AbelianInvariants(B), "\n");
for lam in LinearCharacters(B) do
  ind := InducedClassFunction(lam, G);
  dec := List(irr, chi -> ScalarProduct(T, ind, chi));
  Print("IND_B ker-index ", Index(B, KernelOfCharacter(lam)), " : ", List(Filtered([1..Length(irr)], i -> dec[i] > 0), i -> [lab[i], dec[i]]), "\n");
od;
QUIT;
