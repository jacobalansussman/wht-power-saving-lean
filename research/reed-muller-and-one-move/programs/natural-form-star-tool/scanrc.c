/* scanrc.c (rm2-exact, 15 Oct 2026): fast ten-piece rank-cover filter for the block-5 multi-target system with one
move per nonzero direction and canonical Singer targets T_k = mult(alpha^(2k)) (tables from gen_tables.py).

For a test form F (support D_F, s_F = |D_F| - r_F) and two linear functionals phi != psi, D_F x D_F is cut into
  H x H, H' x H'  (H = D_F n {phi = 0}, H' = D_F n {phi = 1})  and the eight cross blocks
  K x K', K' x K (union D_F n {psi = 0}),  K x L', L' x K (union {phi+psi = 0}),
  L x K', K' x L (union {phi+psi = 1}),  L x L', L' x L (union {psi = 1}),
with K = H n {psi=0}, L = H n {psi=1}, K' = H' n {psi=0}, L' = H' n {psi=1}.  Each piece R x C has
rank E[R,C] <= min(|R|, |C|, s_G) for every form G != F with R u C inside D_G (Theorem 10 with a target per array).
If the ten bounds sum to less than s_F, no program exists (exact; every hit is re-verified independently in Python).

usage: scanrc <tables5.bin> <mode> ...
  mode "list": read target multisets from stdin, one per line "k1,k2,..." (repeats allowed)
  mode "enum <W> <p> <P>": enumerate W-subsets of Z_31 up to k -> 2^j k + c (one representative per orbit), and
        test those with (index % P == p)
  test forms: all 1023 even forms (r_F = 0 for every target), plus with flag ALLF=1 in the environment all forms.
Output: one line per target set: "S k1,..,kW D" (deletion fails), "S ... X F phi psi bound sF G1..G10" (excluded),
        "S ... P" (passes the filter).
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int nf; static uint32_t *mask; static int *wt, *vid; static uint32_t pk[32];
static uint32_t lin[32];
static int nT; static int *Tf; static uint32_t *hs; static uint8_t *minwt; static uint16_t *argG;
#define FULL31 0x7fffffffu

static void load(const char *fn) {
  FILE *f = fopen(fn, "rb"); if (!f) { perror(fn); exit(1); }
  if (fread(&nf, 4, 1, f) != 1) exit(1);
  mask = malloc(4 * nf); wt = malloc(4 * nf); vid = malloc(4 * nf);
  if (fread(mask, 4, nf, f) != (size_t)nf || fread(wt, 4, nf, f) != (size_t)nf || fread(vid, 4, nf, f) != (size_t)nf) exit(1);
  int nv; if (fread(&nv, 4, 1, f) != 1 || nv != 32) exit(1);
  if (fread(pk, 4, 32, f) != 32) exit(1);
  fclose(f);
  for (int c = 1; c < 32; c++) { uint32_t m = 0; for (int z = 1; z < 32; z++) if (__builtin_popcount(c & z) & 1) m |= 1u << (z - 1); lin[c] = m; }
}

static void precompute(int allf) {
  Tf = malloc(4 * nf); nT = 0;
  for (int F = 1; F < nf; F++) if (allf || vid[F] == 0) Tf[nT++] = F;
  hs = malloc(4L * nT * 62); minwt = malloc((long)nT * 62 * 32); argG = malloc(2L * nT * 62 * 32);
  memset(minwt, 255, (long)nT * 62 * 32);
  for (int t = 0; t < nT; t++) {
    int F = Tf[t]; uint32_t A = mask[F];
    for (int c = 1; c < 32; c++) { hs[(long)t * 62 + 2 * (c - 1)] = A & ~lin[c]; hs[(long)t * 62 + 2 * (c - 1) + 1] = A & lin[c]; }
    uint32_t *H = hs + (long)t * 62; int hmin = 99; for (int i = 0; i < 62; i++) { int pc = __builtin_popcount(H[i]); if (pc < hmin) hmin = pc; }
    uint8_t *mw = minwt + (long)t * 62 * 32; uint16_t *ag = argG + (long)t * 62 * 32;
    for (int G = 1; G < nf; G++) {
      if (G == F) continue;
      uint32_t I = mask[G] & A; int v = vid[G], w = wt[G];
      if (__builtin_popcount(I) < hmin) continue;
      for (int i = 0; i < 62; i++) if ((H[i] & ~I) == 0 && w < mw[i * 32 + v]) { mw[i * 32 + v] = w; ag[i * 32 + v] = G; }
    }
  }
}

static inline int imin(int a, int b) { return a < b ? a : b; }

/* returns 0 pass, 1 deletion fail, 2 excluded (cert printed into buf) */
static int star = 0;
static int test(const int *cnt, char *buf) {
  int w[32];
  for (int v = 0; v < 32; v++) { int s = 0; for (int k = 0; k < 31; k++) if (pk[v] >> k & 1) s += cnt[k]; w[v] = s; }
  for (int v = 1; v < 32; v++) if (w[v] > 8) return 1;   /* every class v != 0 has flats (weight 8): s_G < 0 */
  int h[62], hg[62];
  for (int t = 0; t < nT; t++) {
    int F = Tf[t]; int sF = wt[F] - w[vid[F]]; if (sF <= 0) continue;
    uint8_t *mw = minwt + (long)t * 62 * 32; uint16_t *ag = argG + (long)t * 62 * 32;
    for (int i = 0; i < 62; i++) {
      int b = 1000, g = -1;
      for (int v = 0; v < 32; v++) if (mw[i * 32 + v] != 255) { int s = mw[i * 32 + v] - w[v]; if (s < b) { b = s; g = ag[i * 32 + v]; } }
      h[i] = b; hg[i] = g;
    }
    uint32_t A = mask[F];
    if (star && sF == __builtin_popcount(A)) {
      /* star lemma (e3-slice6 7.2, k = 0): E[D_F,D_F] nonsingular; Psi = D_F n {phi=c1, psi=c2}; the three affine
         hyperplanes through that 3-flat split the rest; 2|Psi| <= sum of the three half-set ranks. */
      for (int ph = 1; ph < 32; ph++) for (int ps = ph + 1; ps < 32; ps++) {
        int x = ph ^ ps;
        for (int c1 = 0; c1 < 2; c1++) for (int c2 = 0; c2 < 2; c2++) {
          uint32_t P = A & (c1 ? lin[ph] : ~lin[ph]) & (c2 ? lin[ps] : ~lin[ps]);
          int n = __builtin_popcount(P); if (n == 0) continue;
          int ia = 2 * (ph - 1) + c1, ib = 2 * (ps - 1) + c2, ic = 2 * (x - 1) + (c1 ^ c2);
          int sa = imin(__builtin_popcount(hs[(long)t * 62 + ia]), h[ia]);
          int sb = imin(__builtin_popcount(hs[(long)t * 62 + ib]), h[ib]);
          int sc = imin(__builtin_popcount(hs[(long)t * 62 + ic]), h[ic]);
          if (2 * n > sa + sb + sc) {
            sprintf(buf, "Y %d %d %d %d %d %d %d %d %d %d", F, ph, ps, c1, c2, 2 * n, sa + sb + sc, hg[ia], hg[ib], hg[ic]);
            return 3;
          }
        }
      }
    }
    for (int ph = 1; ph < 32; ph++) {
      int i0 = 2 * (ph - 1), i1 = i0 + 1;
      uint32_t P0 = A & ~lin[ph], P1 = A & lin[ph];
      int n0 = __builtin_popcount(P0), n1 = __builtin_popcount(P1);
      int d0 = imin(n0, h[i0]), d1 = imin(n1, h[i1]);
      if (d0 + d1 >= sF) continue;
      for (int ps = 1; ps < 32; ps++) {
        if (ps == ph) continue;
        int q00 = __builtin_popcount(P0 & ~lin[ps]), q01 = __builtin_popcount(P0 & lin[ps]);
        int q10 = __builtin_popcount(P1 & ~lin[ps]), q11 = __builtin_popcount(P1 & lin[ps]);
        int j0 = 2 * (ps - 1), x0 = 2 * ((ph ^ ps) - 1);
        int c1 = imin(imin(q00, q10), h[j0]);
        int c2 = imin(imin(q00, q11), h[x0]);
        int c3 = imin(imin(q01, q10), h[x0 + 1]);
        int c4 = imin(imin(q01, q11), h[j0 + 1]);
        int bound = d0 + d1 + 2 * (c1 + c2 + c3 + c4);
        if (bound < sF) {
          sprintf(buf, "X %d %d %d %d %d %d %d %d %d %d %d", F, ph, ps, bound, sF, hg[i0], hg[i1], hg[j0], hg[x0], hg[x0 + 1], hg[j0 + 1]);
          return 2;
        }
      }
    }
  }
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

static void report(const int *cnt, int res, const char *buf) {
  printf("S "); int first = 1;
  for (int k = 0; k < 31; k++) for (int r = 0; r < cnt[k]; r++) { printf(first ? "%d" : ",%d", k); first = 0; }
  if (res == 0) printf(" P\n"); else if (res == 1) printf(" D\n"); else printf(" %s\n", buf);
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage\n"); return 1; }
  load(argv[1]); mktab();
  const char *af = getenv("ALLF"); precompute(af && atoi(af));
  const char *st = getenv("STAR"); star = st && atoi(st);
  fprintf(stderr, "test forms %d\n", nT);
  char buf[512]; int cnt[31];
  long ntot = 0, nP = 0, nD = 0, nX = 0, nY = 0;
  if (!strcmp(argv[2], "list")) {
    char line[1024];
    while (fgets(line, sizeof line, stdin)) {
      memset(cnt, 0, sizeof cnt); char *p = line; int any = 0;
      while (*p) { if (*p >= '0' && *p <= '9') { int k = strtol(p, &p, 10); cnt[k % 31]++; any = 1; } else p++; }
      if (!any) continue;
      int res = test(cnt, buf); report(cnt, res, buf); ntot++; if (res == 0) nP++; else if (res == 1) nD++; else if (res == 2) nX++; else nY++;
      fflush(stdout);
    }
  } else {
    int W = atoi(argv[3]), pp = atoi(argv[4]), PP = atoi(argv[5]);
    long idx = 0; int r = W - 1;
    /* S = {0} u (r-subset of 1..30) as bits; Gosper over 30 bits */
    uint32_t c = (1u << r) - 1;
    while (c < (1u << 30)) {
      uint32_t S = 1u | (c << 1);
      if (canonical(S)) {
        if (idx % PP == pp) {
          memset(cnt, 0, sizeof cnt); for (int k = 0; k < 31; k++) if (S >> k & 1) cnt[k] = 1;
          int res = test(cnt, buf); report(cnt, res, buf); ntot++; if (res == 0) nP++; else if (res == 1) nD++; else if (res == 2) nX++; else nY++;
        }
        idx++;
      }
      if (r == 0) break;
      uint32_t u = c & -c, v = c + u; c = v + (((v ^ c) / u) >> 2);
    }
    fprintf(stderr, "classes total %ld\n", idx);
  }
  fprintf(stderr, "tested %ld: pass %ld, deletion %ld, excluded ten-piece %ld, star %ld\n", ntot, nP, nD, nX, nY);
  printf("SUMMARY tested %ld pass %ld deletion %ld excluded %ld star %ld\n", ntot, nP, nD, nX, nY);
  return 0;
}
