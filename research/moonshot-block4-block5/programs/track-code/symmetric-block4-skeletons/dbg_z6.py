from symsolve import *
R = Run(4, 6, 'Z6rot', 1); G = R.G; rng = random.Random(9)
fd = fixed_dirs(R, rng, 4000)
print("|H|", len(R.Hc), "Z6rot fixed dirs", len(fd), "progressing on input:", sum(1 for d in fd if G.progress_letters(G.ID, rng, want=2, Dmat=d)))
