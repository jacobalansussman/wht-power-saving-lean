/* bpm.c (r3-classes, 15 Oct 2026; written fresh, shares no code with starrec2/3/4): block-plus-matching test.
   Input (stdin, text): nf, then nf lines "mask weight slack" (mask: bit z-1 <-> direction z; slack = |D_Q| - r_Q),
   then 31 hex masks Z[x] (bit y-1 set <-> E_xy is a derived zero).
   For every form F with slack == weight (r_F = 0, so E[D_F,D_F] is nonsingular) of weight <= MAXD and every form G
   with S = D_F n D_G and |S| > s_G:  |D_F| <= s_G + termrank(D_F x D_F minus S x S, avoiding zeros)  must hold.
   Prints the first violation (an exclusion) or "pass".                                                          */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static int nf; static uint32_t *msk; static int *wt, *sl; static uint32_t Z[33];
static int mt[33]; static uint32_t seen, curS, curX;
static int aug(int x) {
  uint32_t c = curX & ~Z[x]; if (curS >> (x - 1) & 1) c &= ~curS; c &= ~seen;
  while (c) { int y = __builtin_ctz(c) + 1; c &= c - 1; seen |= 1u << (y - 1);
    if (!mt[y] || aug(mt[y])) { mt[y] = x; return 1; } }
  return 0;
}
static int tr(uint32_t X, uint32_t S) {
  memset(mt, 0, sizeof mt); curX = X; curS = S; int r = 0;
  for (uint32_t T = X; T; T &= T - 1) { seen = 0; r += aug(__builtin_ctz(T) + 1); }
  return r;
}
int main(int argc, char **argv) {
  int maxd = argc > 1 ? atoi(argv[1]) : 20;
  if (scanf("%d", &nf) != 1) return 1;
  msk = malloc(4 * nf); wt = malloc(4 * nf); sl = malloc(4 * nf);
  for (int i = 0; i < nf; i++) if (scanf("%u %d %d", &msk[i], &wt[i], &sl[i]) != 3) return 1;
  for (int x = 1; x <= 31; x++) if (scanf("%x", &Z[x]) != 1) return 1;
  long npairs = 0;
  for (int F = 0; F < nf; F++) {
    if (wt[F] == 0 || sl[F] != wt[F] || wt[F] > maxd) continue;
    uint32_t D = msk[F]; int n = wt[F];
    for (int G = 0; G < nf; G++) {
      if (G == F || wt[G] == 0) continue;
      uint32_t S = D & msk[G]; int k = __builtin_popcount(S);
      if (k <= sl[G]) continue;
      if (2 * (n - k) + sl[G] < n) { printf("EXCLUDED F=%d G=%d n=%d sG=%d tr<=%d (no zeros needed) S=%u\n", F, G, n, sl[G], 2 * (n - k), S); return 0; }
      npairs++;
      int t = tr(D, S);
      if (n > sl[G] + t) { printf("EXCLUDED F=%d G=%d n=%d sG=%d tr=%d S=%u D=%u\n", F, G, n, sl[G], t, S, D); return 0; }
    }
  }
  printf("pass (%ld pairs)\n", npairs);
  return 0;
}
