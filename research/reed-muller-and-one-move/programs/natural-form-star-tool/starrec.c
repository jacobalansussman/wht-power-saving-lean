/* starrec.c (rm2-exact, 15 Oct 2026): the star-lemma recursion (e3-slice6 section 7.2) on the block-5 multi-target
system, one move per nonzero direction, canonical Singer targets T_k = mult(alpha^(2k)) (tables5.bin).

For an even form F (r_F = 0), E[D_F, D_F] is nonsingular (Theorem 10).  For every affine flat A of F_2^5 put
S_A = D_F n A and let u(A) be an upper bound on rank E[S_A, S_A]:
  start: u(A) = min(|S_A|, s_G for every form G != F with S_A inside D_G);
  star lemma: if A has dimension d <= 3, the g = 2^(5-d) - 1 flats B of dimension d+1 through A split the rest of
  the space, so with k = rank E[S_A, S_A] (and a nonsingular k x k block inside, Schur complement 0 on the rest):
      (g - 2) k <= sum_B u(B) - 2 |S_A|.
  A negative right side is a contradiction (no program); otherwise u(A) <= floor((sum_B u(B) - 2|S_A|)/(g-2)),
  and u(A) <= u(B) for B > A.  Dimensions are processed 3, 2, 1, 0 (each uses the finished level above).
usage: starrec <tables5.bin> list < sets   |   starrec <tables5.bin> enum <W> <p> <P>
Output per set: "S ks Z F d flat tot 2|S_A|" (contradiction) or "S ks P minval" (passes; minval = least value of
sum_B u(B) - 2|S_A| met), or "S ks D" (deletion fails).
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int nf; static uint32_t *mask; static int *wt, *vid; static uint32_t pk[32];
#define FULL31 0x7fffffffu
static int nfl[5]; static uint32_t *fl[5]; static int *up[4]; /* up[d][i*g + j] = index of (d+1)-flat */
static int G_[5] = {31, 15, 7, 3, 1};
static int nT; static int *Tf; static uint8_t *mwf; static uint8_t *szf; static int NFL;
static int off[5];

static void load(const char *fn) {
  FILE *f = fopen(fn, "rb"); if (!f) { perror(fn); exit(1); }
  if (fread(&nf, 4, 1, f) != 1) exit(1);
  mask = malloc(4 * nf); wt = malloc(4 * nf); vid = malloc(4 * nf);
  if (fread(mask, 4, nf, f) != (size_t)nf || fread(wt, 4, nf, f) != (size_t)nf || fread(vid, 4, nf, f) != (size_t)nf) exit(1);
  int nv; if (fread(&nv, 4, 1, f) != 1 || nv != 32) exit(1);
  if (fread(pk, 4, 32, f) != 32) exit(1);
  fclose(f);
}

/* flats as 32-bit masks over points 0..31 */
static void mkflats(void) {
  for (int d = 0; d <= 4; d++) {
    int cap = 4096; fl[d] = malloc(4 * cap); nfl[d] = 0;
    /* linear subspaces of dim d: closure of generator sets; collect unique */
    uint32_t *subs = malloc(4 * 4096); int ns = 0;
    int gens[5];
    /* enumerate d-tuples of increasing nonzero vectors, span, keep if dim d */
    for (long code = 0; code < (1L << (5 * d)); code++) {
      int ok = 1; long c = code;
      for (int i = 0; i < d; i++) { gens[i] = c & 31; c >>= 5; if (gens[i] == 0 || (i && gens[i] <= gens[i - 1])) { ok = 0; break; } }
      if (!ok) continue;
      uint32_t S = 1u; /* point 0 */
      for (int i = 0; i < d; i++) { uint32_t T = S; for (int p = 0; p < 32; p++) if (S >> p & 1) T |= 1u << (p ^ gens[i]); S = T; }
      if (__builtin_popcount(S) != (1 << d)) continue;
      int dup = 0; for (int i = 0; i < ns; i++) if (subs[i] == S) { dup = 1; break; }
      if (!dup) subs[ns++] = S;
    }
    for (int i = 0; i < ns; i++) for (int a = 0; a < 32; a++) {
      uint32_t T = 0; for (int p = 0; p < 32; p++) if (subs[i] >> p & 1) T |= 1u << (p ^ a);
      int dup = 0; for (int j = 0; j < nfl[d]; j++) if (fl[d][j] == T) { dup = 1; break; }
      if (!dup) fl[d][nfl[d]++] = T;
    }
    free(subs);
  }
  for (int d = 0; d < 4; d++) {
    int g = G_[d]; up[d] = malloc(4 * nfl[d] * g);
    for (int i = 0; i < nfl[d]; i++) {
      int c = 0;
      for (int j = 0; j < nfl[d + 1]; j++) if ((fl[d][i] & ~fl[d + 1][j]) == 0) { if (c < g) up[d][i * g + c] = j; c++; }
      if (c != g) { fprintf(stderr, "incidence error d=%d c=%d\n", d, c); exit(3); }
    }
  }
  NFL = 0; for (int d = 0; d <= 4; d++) { off[d] = NFL; NFL += nfl[d]; }
  fprintf(stderr, "flats %d %d %d %d %d\n", nfl[0], nfl[1], nfl[2], nfl[3], nfl[4]);
}

static void precompute(int maxd) {
  Tf = malloc(4 * nf); nT = 0;
  for (int F = 1; F < nf; F++) if (vid[F] == 0 && wt[F] <= maxd) Tf[nT++] = F;
  mwf = malloc((long)nT * NFL * 32); szf = malloc((long)nT * NFL);
  int dirs[32];
  for (int t = 0; t < nT; t++) {
    int F = Tf[t]; uint32_t A = mask[F]; int d = 0;
    for (int z = 1; z < 32; z++) if (A >> (z - 1) & 1) dirs[d++] = z;
    long N = 1L << d; uint8_t *T = malloc(N * 32); memset(T, 255, N * 32);
    for (int G = 1; G < nf; G++) {
      if (G == F) continue;
      uint32_t I = mask[G] & A; uint32_t u = 0;
      for (int k = 0; k < d; k++) if (I >> (dirs[k] - 1) & 1) u |= 1u << k;
      uint8_t *e = T + (long)u * 32 + vid[G]; if (wt[G] < *e) *e = wt[G];
    }
    for (int b = 0; b < d; b++) for (long u = 0; u < N; u++) if (!(u >> b & 1)) {
      uint8_t *x = T + u * 32, *y = T + (u | (1L << b)) * 32;
      for (int v = 0; v < 32; v++) x[v] = x[v] < y[v] ? x[v] : y[v];
    }
    for (int dd = 0; dd <= 4; dd++) for (int i = 0; i < nfl[dd]; i++) {
      uint32_t P = fl[dd][i] >> 1;   /* points 1..31 -> bits 0..30 (direction z <-> bit z-1) */
      uint32_t S = P & A; uint32_t u = 0;
      for (int k = 0; k < d; k++) if (S >> (dirs[k] - 1) & 1) u |= 1u << k;
      memcpy(mwf + ((long)t * NFL + off[dd] + i) * 32, T + (long)u * 32, 32);
      szf[(long)t * NFL + off[dd] + i] = __builtin_popcount(S);
    }
    free(T);
  }
}

static int dump = 0;
static int test(const int *cnt, char *buf, int *minval_out) {
  int w[32];
  for (int v = 0; v < 32; v++) { int s = 0; for (int k = 0; k < 31; k++) if (pk[v] >> k & 1) s += cnt[k]; w[v] = s; }
  for (int v = 1; v < 32; v++) if (w[v] > 8) return 1;
  static int u[4096];
  int minval = 1000;
  for (int t = 0; t < nT; t++) {
    uint8_t *mw = mwf + (long)t * NFL * 32; uint8_t *sz = szf + (long)t * NFL;
    for (int i = 0; i < NFL; i++) {
      int b = sz[i];
      for (int v = 0; v < 32; v++) if (mw[i * 32 + v] != 255) { int s = mw[i * 32 + v] - w[v]; if (s < b) b = s; }
      u[i] = b;
    }
    for (int d = 3; d >= 0; d--) {
      int g = G_[d];
      for (int i = 0; i < nfl[d]; i++) {
        int tot = 0, mn = 1000;
        for (int j = 0; j < g; j++) { int x = u[off[d + 1] + up[d][i * g + j]]; tot += x; if (x < mn) mn = x; }
        int n = sz[off[d] + i]; if (n == 0) continue;
        int val = tot - 2 * n;
        if (val < minval) minval = val;
        if (val < 0) { sprintf(buf, "Z %d %d %d %d %d", Tf[t], d, i, tot, 2 * n); *minval_out = val; return 2; }
        int k = val / (g - 2);
        int *ui = &u[off[d] + i];
        if (k < *ui) *ui = k;
        if (mn < *ui) *ui = mn;
        if (dump && *ui == 0) printf("ZB %d %d %d %u %d\n", Tf[t], d, i, (fl[d][i] >> 1) & mask[Tf[t]], val);
      }
    }
  }
  *minval_out = minval;
  return 0;
}

static uint32_t ptab[5][4][256];
static void mktab(void) {
  for (int j = 0; j < 5; j++) for (int b = 0; b < 4; b++) for (int x = 0; x < 256; x++) {
    uint32_t r = 0;
    for (int t = 0; t < 8; t++) if (x >> t & 1) { int k = 8 * b + t; if (k < 31) r |= 1u << (((1 << j) * k) % 31); }
    ptab[j][b][x] = r;
  }
}
static inline uint32_t mulj(int j, uint32_t S) { return ptab[j][0][S & 255] | ptab[j][1][S >> 8 & 255] | ptab[j][2][S >> 16 & 255] | ptab[j][3][S >> 24 & 255]; }
static inline uint32_t rotsub(uint32_t S, int e) { uint64_t x = S; return (uint32_t)(((x >> e) | (x << (31 - e))) & FULL31); }
static int canonical(uint32_t S) {
  for (int j = 0; j < 5; j++) {
    uint32_t M = mulj(j, S), T = M;
    while (T) { int e = __builtin_ctz(T); T &= T - 1; if (rotsub(M, e) < S) return 0; }
  }
  return 1;
}
static void report(const int *cnt, int res, const char *buf, int mv) {
  printf("S "); int first = 1;
  for (int k = 0; k < 31; k++) for (int r = 0; r < cnt[k]; r++) { printf(first ? "%d" : ",%d", k); first = 0; }
  if (res == 0) printf(" P %d\n", mv); else if (res == 1) printf(" D\n"); else printf(" %s\n", buf);
}

int main(int argc, char **argv) {
  load(argv[1]); mktab(); mkflats();
  const char *md = getenv("MAXD"); precompute(md ? atoi(md) : 20);
  const char *du = getenv("DUMP"); dump = du && atoi(du);
  fprintf(stderr, "test forms %d\n", nT);
  char buf[256]; int cnt[31]; long ntot = 0, nP = 0, nD = 0, nZ = 0; int mv;
  if (!strcmp(argv[2], "list")) {
    char line[1024];
    while (fgets(line, sizeof line, stdin)) {
      memset(cnt, 0, sizeof cnt); char *p = line; int any = 0;
      while (*p) { if (*p >= '0' && *p <= '9') { int k = strtol(p, &p, 10); cnt[k % 31]++; any = 1; } else p++; }
      if (!any) continue;
      int res = test(cnt, buf, &mv); report(cnt, res, buf, mv); ntot++; if (res == 0) nP++; else if (res == 1) nD++; else nZ++;
      fflush(stdout);
    }
  } else {
    int W = atoi(argv[3]), pp = atoi(argv[4]), PP = atoi(argv[5]);
    long idx = 0; int r = W - 1; uint32_t c = (1u << r) - 1;
    while (c < (1u << 30)) {
      uint32_t S = 1u | (c << 1);
      if (canonical(S)) {
        if (idx % PP == pp) {
          memset(cnt, 0, sizeof cnt); for (int k = 0; k < 31; k++) if (S >> k & 1) cnt[k] = 1;
          int res = test(cnt, buf, &mv); report(cnt, res, buf, mv); ntot++; if (res == 0) nP++; else if (res == 1) nD++; else nZ++;
        }
        idx++;
      }
      if (r == 0) break;
      uint32_t u = c & -c, v = c + u; c = v + (((v ^ c) / u) >> 2);
    }
  }
  printf("SUMMARY tested %ld pass %ld deletion %ld star-recursion %ld\n", ntot, nP, nD, nZ);
  return 0;
}
