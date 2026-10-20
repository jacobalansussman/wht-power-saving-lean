# The hook values behind the 640 near miss

**What it shows.** The k = 8 hook certificate at block 4 (the module V (x) Lambda^8 V with f = 9 digit columns) values
one paid move at 27,456 and a depth-3 piece at 68,320. Five moves then allow 5 x 27,456 = 137,280, two depth-3 pieces
need 2 x 68,320 = 136,640, and 640 is left. `chirank.py` computes a hook value as a sum of ranks modulo the prime
33,554,393 over the 16 translation characters, ranking one block per orbit of the address maps that fix the element.

**Commands.** Build the rank program: `cc -O2 -o rankp rankp.c` (add `-fopenmp` where OpenMP is available). Then
- `HOOKK=8 CHIS=0 python3 chirank.py d1,d3` ranks the character-0 blocks of one move (d1) and of a depth-3 element (d3);
- `HOOKK=8 CHIS=1 python3 chirank.py d1` and `HOOKK=8 CHIS=7 python3 chirank.py d3` rank one more character each;
- without `CHIS` it ranks every orbit and prints the full value; `HOOKK=k` selects the hook (1 to 15).

Needs numpy. Each block is a 6,435 x 6,435 matrix; one block took 54 to 75 seconds on an Apple M4 Mac mini (16 GB) and stayed under 3 GB.

**Expected output.** `rerun_k8_2026-10-19.txt`: rank 1,716 for d1 and 4,270 for d3 at character 0, and the same ranks at
characters 1 and 7. At even k every character block has the same rank (argued in the track's notes, and checked here at
two characters), so the values are 16 x 1,716 = 27,456 and 16 x 4,270 = 68,320.

**Re-run for this revision:** yes, on 19 October 2026 (the four blocks above). The odd-k values (k = 7: 2.4591; k = 9:
2.4591; k = 5: 2.2552) were not re-run; at odd k the 16 blocks differ, and a shortcut that multiplied one block by 16
gave wrong values there in an earlier round.
