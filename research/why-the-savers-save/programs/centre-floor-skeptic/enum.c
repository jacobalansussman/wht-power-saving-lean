/* v-centre / enum.c: own enumeration of every subspace W of F2^9 (RREF, highest bit = pivot).
   For each W: L(W) = ports with label in W.  Outputs per dim d:
     n(d)    = max |L(W)|
     Tmax(d) = max #frames meeting L(W)
     cov[c][d] = max |L(W) & coflat c|   (c over the 120 E7-type co-flats)
   and dumps the distinct L(W) masks with label-closure exactly d for d = 6, 7 (for exact-cover checks).
   Reads data.txt from build.py.  Own code. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef struct { uint64_t lo, hi; } m128;
static int lab[112], portof[512];
static int nfr; static m128 fr[1300];
static int ncf; static m128 cf[200];
static long nsub[10]; static int nmax[10]; static int Tmax[10]; static int cov[200][10];
static FILE *dump6, *dump7;
static int rank9(int *v, int k) { int b[9] = {0}; int r = 0;
  for (int i = 0; i < k; i++) { int x = v[i];
    for (int p = 8; p >= 0; p--) if (x >> p & 1) { if (b[p]) x ^= b[p]; else { b[p] = x; r++; break; } } }
  return r; }
static void visit(int d, int *bas) {
  m128 L = {0, 0}; int cnt = 0; int labsin[112]; int nl = 0;
  int span = 1 << d;
  for (int c = 0; c < span; c++) { int x = 0; for (int i = 0; i < d; i++) if (c >> i & 1) x ^= bas[i];
    int p = portof[x]; if (p >= 0) { if (p < 64) L.lo |= 1ULL << p; else L.hi |= 1ULL << (p - 64); cnt++; labsin[nl++] = x; } }
  nsub[d]++;
  if (cnt > nmax[d]) nmax[d] = cnt;
  int T = 0; for (int f = 0; f < nfr; f++) if ((fr[f].lo & L.lo) | (fr[f].hi & L.hi)) T++;
  if (T > Tmax[d]) Tmax[d] = T;
  for (int c = 0; c < ncf; c++) { int k = __builtin_popcountll(cf[c].lo & L.lo) + __builtin_popcountll(cf[c].hi & L.hi);
    if (k > cov[c][d]) cov[c][d] = k; }
  if ((d == 6 || d == 7) && rank9(labsin, nl) == d) fprintf(d == 6 ? dump6 : dump7, "%016llx%016llx\n",
      (unsigned long long)L.hi, (unsigned long long)L.lo);
}
int main(void) {
  FILE *f = fopen("data.txt", "r"); if (!f) return 1;
  for (int i = 0; i < 512; i++) portof[i] = -1;
  for (int t = 0; t < 112; t++) { if (fscanf(f, "%d", &lab[t]) != 1) return 2; portof[lab[t]] = t; }
  if (fscanf(f, "%d", &nfr) != 1) return 3;
  for (int i = 0; i < nfr; i++) { fr[i].lo = fr[i].hi = 0; for (int j = 0; j < 8; j++) { int p; if (fscanf(f, "%d", &p) != 1) return 4;
      if (p < 64) fr[i].lo |= 1ULL << p; else fr[i].hi |= 1ULL << (p - 64); } }
  if (fscanf(f, "%d", &ncf) != 1) return 5;
  { char line[4096]; fgets(line, sizeof line, f);
    for (int i = 0; i < ncf; i++) { cf[i].lo = cf[i].hi = 0; if (!fgets(line, sizeof line, f)) return 6; char *s = line; int p, k;
      while (sscanf(s, "%d%n", &p, &k) == 1) { s += k; if (p < 64) cf[i].lo |= 1ULL << p; else cf[i].hi |= 1ULL << (p - 64); } } }
  fclose(f);
  dump6 = fopen("L6.txt", "w"); dump7 = fopen("L7.txt", "w");
  /* enumerate RREF: pivot set P (bitmask over 9 positions), each basis vector i: bit p_i set, zeros at other pivots,
     free bits at non-pivot positions below p_i */
  for (int P = 0; P < 512; P++) {
    int d = __builtin_popcount(P); int pv[9], k = 0; for (int p = 0; p < 9; p++) if (P >> p & 1) pv[k++] = p;
    int fpos[9][9], nf[9], tot = 0;
    for (int i = 0; i < d; i++) { nf[i] = 0; for (int q = 0; q < pv[i]; q++) if (!(P >> q & 1)) fpos[i][nf[i]++] = q; tot += nf[i]; }
    for (long a = 0; a < (1L << tot); a++) { int bas[9]; long b = a;
      for (int i = 0; i < d; i++) { int v = 1 << pv[i]; for (int j = 0; j < nf[i]; j++) { if (b & 1) v |= 1 << fpos[i][j]; b >>= 1; } bas[i] = v; }
      visit(d, bas); }
  }
  fclose(dump6); fclose(dump7);
  long all = 0; for (int d = 0; d <= 9; d++) all += nsub[d];
  printf("subspaces total %ld\n", all);
  for (int d = 0; d <= 9; d++) printf("d=%d  #W=%ld  n(d)=%d  Tmax(d)=%d\n", d, nsub[d], nmax[d], Tmax[d]);
  double best = 1e9; int bd = 0; for (int d = 1; d <= 9; d++) { double r = 8.0 * nfr * d / Tmax[d]; if (r < best) { best = r; bd = d; } }
  printf("clique floor = %.4f at d = %d\n", best, bd);
  /* coverage maxima: report min/max over Gamma (0..111) and N (112..119) */
  for (int d = 1; d <= 9; d++) { int gmn = 999, gmx = 0, nmn = 999, nmx = 0;
    for (int c = 0; c < 112; c++) { if (cov[c][d] < gmn) gmn = cov[c][d]; if (cov[c][d] > gmx) gmx = cov[c][d]; }
    for (int c = 112; c < ncf; c++) { if (cov[c][d] < nmn) nmn = cov[c][d]; if (cov[c][d] > nmx) nmx = cov[c][d]; }
    printf("cov d=%d  Gamma %d..%d   N %d..%d\n", d, gmn, gmx, nmn, nmx); }
  return 0;
}
