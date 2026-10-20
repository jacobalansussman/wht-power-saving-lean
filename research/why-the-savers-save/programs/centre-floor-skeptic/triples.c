/* v-centre / triples.c: every triple of distinct dim-6 closed label sets S1,S2,S3 (from enum.c, L6u.txt) with
   |C \ (S1 u S2 u S3)| <= EPS for some E7-type co-flat C (Gamma_a, N_i; all 120, no symmetry used).
   Writes the triples (global ids, sorted) to stdout as "a b c"; counts to stderr.  Own code. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef struct { uint64_t lo, hi; } m128;
static m128 *S; static int nS;
static m128 cf[200]; static int ncf;
static inline int pc(m128 a) { return __builtin_popcountll(a.lo) + __builtin_popcountll(a.hi); }
static inline m128 and_(m128 a, m128 b) { m128 r = {a.lo & b.lo, a.hi & b.hi}; return r; }
static inline m128 or_(m128 a, m128 b) { m128 r = {a.lo | b.lo, a.hi | b.hi}; return r; }
static int cmpi(const void *a, const void *b) { int x = *(int *)a, y = *(int *)b; return x - y; }
int main(int argc, char **argv) {
  int EPS = atoi(argv[1]);
  FILE *f = fopen("L6u.txt", "r"); char buf[64]; int cap = 1 << 20; S = malloc(sizeof(m128) * cap); nS = 0;
  while (fscanf(f, "%32s", buf) == 1) { char h[17], l[17]; memcpy(h, buf, 16); h[16] = 0; memcpy(l, buf + 16, 16); l[16] = 0;
    S[nS].hi = strtoull(h, 0, 16); S[nS].lo = strtoull(l, 0, 16); nS++; }
  fclose(f);
  f = fopen("data.txt", "r"); int x; for (int i = 0; i < 112; i++) fscanf(f, "%d", &x);
  int nfr; fscanf(f, "%d", &nfr); for (int i = 0; i < nfr * 8; i++) fscanf(f, "%d", &x);
  fscanf(f, "%d", &ncf); char line[4096]; fgets(line, sizeof line, f);
  for (int i = 0; i < ncf; i++) { cf[i].lo = cf[i].hi = 0; fgets(line, sizeof line, f); char *s = line; int p, k;
    while (sscanf(s, "%d%n", &p, &k) == 1) { s += k; if (p < 64) cf[i].lo |= 1ULL << p; else cf[i].hi |= 1ULL << (p - 64); } }
  fclose(f);
  fprintf(stderr, "dim-6 sets %d, coflats %d\n", nS, ncf);
  int *cand = malloc(sizeof(int) * nS); int *cv = malloc(sizeof(int) * nS);
  long total = 0;
  for (int c = 0; c < ncf; c++) {
    m128 C = cf[c]; int sz = pc(C), nc = 0, mx = 0;
    for (int i = 0; i < nS; i++) { int k = pc(and_(S[i], C)); if (k > mx) mx = k; }
    for (int i = 0; i < nS; i++) { int k = pc(and_(S[i], C)); if (k >= sz - EPS - 2 * mx) { cand[nc] = i; cv[nc] = k; nc++; } }
    long cnt = 0;
    for (int a = 0; a < nc; a++) for (int b = a + 1; b < nc; b++) {
      if (cv[a] + cv[b] < sz - EPS - mx) continue;
      m128 U = and_(or_(S[cand[a]], S[cand[b]]), C); int u = pc(U);
      if (u < sz - EPS - mx) continue;
      for (int e = b + 1; e < nc; e++) {
        if (u + cv[e] < sz - EPS) continue;
        int w = pc(and_(or_(U, S[cand[e]]), C));
        if (w >= sz - EPS) { int t[3] = {cand[a], cand[b], cand[e]}; qsort(t, 3, sizeof(int), cmpi);
          printf("%d %d %d\n", t[0], t[1], t[2]); cnt++; } } }
    total += cnt;
    fprintf(stderr, "coflat %d size %d max-cover %d candidates %d triples %ld\n", c, sz, mx, nc, cnt);
  }
  fprintf(stderr, "TOTAL triple-incidences %ld\n", total);
  return 0;
}
