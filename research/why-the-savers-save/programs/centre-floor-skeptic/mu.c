/* v-centre / mu.c: own exact maximum over every subspace W of F2^9 of H(W)/dim W, H(W) = sum of integer weights
   c_S over pool sets S (8 lines, bases of R^8) meeting L(W).  Floor = 8 * sum c / max(H/d).  Own code.
   Input: labT.txt (labels of the 112 lines, from mappool.py) and pool_mine.txt (c_S and 8 line indices). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXP 11000
#define PW ((MAXP + 63) / 64)
static int lab[112], portof[512], P; static int64_t c[MAXP]; static uint64_t inc[112][PW]; static int pw;
static int64_t Hb = 0; static int db = 1; static long nsub = 0;
int main(void) {
  FILE *f = fopen("labT.txt", "r"); for (int i = 0; i < 512; i++) portof[i] = -1;
  for (int t = 0; t < 112; t++) { fscanf(f, "%d", &lab[t]); portof[lab[t]] = t; } fclose(f);
  f = fopen("pool_mine.txt", "r"); fscanf(f, "%d", &P); pw = (P + 63) / 64;
  int64_t C = 0;
  for (int s = 0; s < P; s++) { long long ci; fscanf(f, "%lld", &ci); c[s] = ci; C += ci;
    for (int j = 0; j < 8; j++) { int t; fscanf(f, "%d", &t); inc[t][s / 64] |= 1ULL << (s % 64); } }
  fclose(f);
  for (int Pm = 1; Pm < 512; Pm++) {
    int d = __builtin_popcount(Pm); int pv[9], k = 0; for (int p = 0; p < 9; p++) if (Pm >> p & 1) pv[k++] = p;
    int fpos[9][9], nf[9], tot = 0;
    for (int i = 0; i < d; i++) { nf[i] = 0; for (int q = 0; q < pv[i]; q++) if (!(Pm >> q & 1)) fpos[i][nf[i]++] = q; tot += nf[i]; }
    for (long a = 0; a < (1L << tot); a++) { int bas[9]; long b = a;
      for (int i = 0; i < d; i++) { int v = 1 << pv[i]; for (int j = 0; j < nf[i]; j++) { if (b & 1) v |= 1 << fpos[i][j]; b >>= 1; } bas[i] = v; }
      nsub++;
      uint64_t m[PW]; memset(m, 0, sizeof(uint64_t) * pw); int any = 0;
      for (int cc = 1; cc < (1 << d); cc++) { int x = 0; for (int i = 0; i < d; i++) if (cc >> i & 1) x ^= bas[i];
        int t = portof[x]; if (t >= 0) { any = 1; for (int w = 0; w < pw; w++) m[w] |= inc[t][w]; } }
      if (!any) continue;
      int64_t H = 0;
      for (int w = 0; w < pw; w++) { uint64_t x = m[w]; while (x) { int bit = __builtin_ctzll(x); H += c[w * 64 + bit]; x &= x - 1; } }
      if ((__int128)H * db > (__int128)Hb * d) { Hb = H; db = d; }
    }
  }
  printf("subspaces %ld  C=%lld  best d=%d H=%lld  floor = 8*C*d/H = %.6f\n", nsub, (long long)C, db, (long long)Hb, 8.0 * C * db / (double)Hb);
  return 0;
}
