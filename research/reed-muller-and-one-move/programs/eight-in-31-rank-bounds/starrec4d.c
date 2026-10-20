/* starrec4.c (r3-classes, 15 Oct 2026): starrec3 with LB=2: the top-down lower bounds use the zero cells:
   rank E[X,X] <= rank E[S,S] + (largest matching of X x X minus S x S avoiding derived zeros), for S = S_A inside X = D_F
   (4- and 3-flats) and X = S_B for the flats B above A; and odd forms G with D_G inside D_F: rank E[S_A] >= s_G - matching(D_G^2 \ (D_G n S_A)^2).
   starrec3.c: starrec2 + LOWER-BOUND propagation (env LB=1, default on):
   ell(A) <= rank E[S_A,S_A] from (i) the top, rank E[D_F,D_F] = |D_F| - nu, and deletion: rank E[S_A] >= rank E[S_B] - 2|S_B \ S_A|;
   (ii) monotonicity rank E[S_B] >= rank E[S_A] for A in B; (iii) STAR with k >= ell(A): sum_j r(B_j) >= (g-2) ell(A) + 2|S_A| - 2 nu, so
   r(B_j) >= that minus sum_{i != j} u(B_i).  A flat with ell > u is a contradiction (output 'L').  Then as starrec2.
   Original header: starrec2.c (rm2-exact, 15 Oct 2026): star-lemma recursion + pooled zero blocks + clique and matching tests, for the
multi-target system at block m (3, 4, 5) with ONE move per nonzero direction and canonical Singer targets
T_k = mult(alpha^(2k)) (tables<m>.bin from gen_tables.py; k = 0 is the identity).

Facts used (all proved earlier; see NOTE.md):
 (T10) rank E[D_Q, D_Q] = |D_Q| - r_Q for every form Q (Theorem 10 with a target per array).
 (STAR) e3-slice6 7.2: E nonsingular on X, X = Psi + Q_1 + ... + Q_g, k = rank E[Psi,Psi]  =>
        (g - 2) k <= sum_j rank E[Psi u Q_j, Psi u Q_j] - 2|Psi|.
 (CLQ)  E[Y,Y] = 0  =>  |Y| <= W   (E = I + Gamma - K, Gamma strictly lower triangular, rank K <= W: ex2-families Thm 1).
 (MAT)  rank E[D,D] <= the largest matching of D x D avoiding cells known to be 0 (Corollary 10').
Procedure, per target multiset: for each form F with r_F = 0 (E[D_F,D_F] nonsingular) and each affine flat A of
F_2^m, u(A) bounds rank E[D_F n A, D_F n A]: start from |S_A|, s_G of forms G != F containing S_A, and (round >= 2)
the matching bound with the pooled zero cells; then STAR from dimension m-2 down to 0.  A negative right side is a
contradiction.  Every u(A) = 0 is a zero block of the GLOBAL matrix E; zero blocks are pooled over all F, then CLQ
and MAT (MAT on every form Q: matching >= s_Q needed).  Rounds repeat while new zero cells appear (max ROUNDS).
usage: starrec2 <tables.bin> list < sets      (sets: "k1,k2,..." per line, repeats allowed)
env: MAXD (max |D_F| of test forms, default 20), ROUNDS (default 3), CERT=1 (print certificate details)
Output: "S ks P"  |  "S ks D"  |  "S ks Z round F d flat tot 2|S|"  |  "S ks C round |Y| W"  |  "S ks M round F mm s_F"
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int m, nd, npt, nf, nv; static uint32_t *mask; static int *wt, *vid; static uint32_t *pk;
static int nfl[6]; static uint32_t *fl[6]; static int *up[6]; static int G_[6]; static int off[6]; static int NFL;
static int nT; static int *Tf; static uint8_t *mwf; static uint8_t *szf; static uint32_t *smf; /* S_A as direction mask */
static int minwv[64];

static void load(const char *fn) {
  FILE *f = fopen(fn, "rb"); if (!f) { perror(fn); exit(1); }
  if (fread(&nf, 4, 1, f) != 1) exit(1);
  mask = malloc(4 * nf); wt = malloc(4 * nf); vid = malloc(4 * nf);
  if (fread(mask, 4, nf, f) != (size_t)nf || fread(wt, 4, nf, f) != (size_t)nf || fread(vid, 4, nf, f) != (size_t)nf) exit(1);
  if (fread(&nv, 4, 1, f) != 1) exit(1);
  pk = malloc(4 * nv); if (fread(pk, 4, nv, f) != (size_t)nv) exit(1);
  fclose(f);
  m = 0; while ((1 << (m * (m + 1) / 2)) < nf) m++;
  nd = (1 << m) - 1; npt = 1 << m;
  if (nv != npt) { fprintf(stderr, "nv %d != 2^m\n", nv); exit(1); }
  for (int v = 0; v < nv; v++) minwv[v] = 1000;
  for (int G = 1; G < nf; G++) if (wt[G] < minwv[vid[G]]) minwv[vid[G]] = wt[G];
}

static void mkflats(void) {
  for (int d = 0; d < m; d++) {
    fl[d] = malloc(4 * 8192); nfl[d] = 0;
    uint32_t *subs = malloc(4 * 8192); int ns = 0; int gens[6];
    for (long code = 0; code < (1L << (m * d)); code++) {
      int ok = 1; long c = code;
      for (int i = 0; i < d; i++) { gens[i] = c & nd; c >>= m; if (gens[i] == 0 || (i && gens[i] <= gens[i - 1])) { ok = 0; break; } }
      if (!ok) continue;
      uint32_t S = 1u;
      for (int i = 0; i < d; i++) { uint32_t T = S; for (int p = 0; p < npt; p++) if (S >> p & 1) T |= 1u << (p ^ gens[i]); S = T; }
      if (__builtin_popcount(S) != (1 << d)) continue;
      int dup = 0; for (int i = 0; i < ns; i++) if (subs[i] == S) { dup = 1; break; }
      if (!dup) subs[ns++] = S;
    }
    for (int i = 0; i < ns; i++) for (int a = 0; a < npt; a++) {
      uint32_t T = 0; for (int p = 0; p < npt; p++) if (subs[i] >> p & 1) T |= 1u << (p ^ a);
      int dup = 0; for (int j = 0; j < nfl[d]; j++) if (fl[d][j] == T) { dup = 1; break; }
      if (!dup) fl[d][nfl[d]++] = T;
    }
    free(subs);
  }
  for (int d = 0; d < m; d++) G_[d] = (1 << (m - d)) - 1;
  for (int d = 0; d + 1 < m; d++) {
    int g = G_[d]; up[d] = malloc(4 * nfl[d] * g);
    for (int i = 0; i < nfl[d]; i++) {
      int c = 0;
      for (int j = 0; j < nfl[d + 1]; j++) if ((fl[d][i] & ~fl[d + 1][j]) == 0) { if (c < g) up[d][i * g + c] = j; c++; }
      if (c != g) { fprintf(stderr, "incidence error d=%d c=%d g=%d\n", d, c, g); exit(3); }
    }
  }
  NFL = 0; for (int d = 0; d < m; d++) { off[d] = NFL; NFL += nfl[d]; }
  fprintf(stderr, "m=%d flats:", m); for (int d = 0; d < m; d++) fprintf(stderr, " %d", nfl[d]); fprintf(stderr, "\n");
}

static void precompute(int maxd, uint32_t vmask_lo, uint32_t vmask_hi) {
  Tf = malloc(4 * nf); nT = 0;
  for (int F = 1; F < nf; F++) {
    int v = vid[F]; int inV = v < 32 ? (vmask_lo >> v & 1) : (vmask_hi >> (v - 32) & 1);
    if (inV && wt[F] <= maxd) Tf[nT++] = F;
  }
  mwf = malloc((long)nT * NFL * nv); szf = malloc((long)nT * NFL); smf = malloc(4L * nT * NFL);
  int dirs[32];
  for (int t = 0; t < nT; t++) {
    int F = Tf[t]; uint32_t A = mask[F]; int d = 0;
    for (int z = 1; z <= nd; z++) if (A >> (z - 1) & 1) dirs[d++] = z;
    long N = 1L << d; uint8_t *T = malloc(N * nv); memset(T, 255, N * nv);
    for (int G = 1; G < nf; G++) {
      if (G == F) continue;
      uint32_t I = mask[G] & A; uint32_t u = 0;
      for (int k = 0; k < d; k++) if (I >> (dirs[k] - 1) & 1) u |= 1u << k;
      uint8_t *e = T + (long)u * nv + vid[G]; if (wt[G] < *e) *e = wt[G];
    }
    for (int b = 0; b < d; b++) for (long u = 0; u < N; u++) if (!(u >> b & 1)) {
      uint8_t *x = T + u * nv, *y = T + (u | (1L << b)) * nv;
      for (int v = 0; v < nv; v++) x[v] = x[v] < y[v] ? x[v] : y[v];
    }
    for (int dd = 0; dd < m; dd++) for (int i = 0; i < nfl[dd]; i++) {
      uint32_t S = (fl[dd][i] >> 1) & A; uint32_t u = 0;
      for (int k = 0; k < d; k++) if (S >> (dirs[k] - 1) & 1) u |= 1u << k;
      long idx = (long)t * NFL + off[dd] + i;
      memcpy(mwf + idx * nv, T + (long)u * nv, nv);
      szf[idx] = __builtin_popcount(S); smf[idx] = S;
    }
    free(T);
  }
}

/* zero cells: Zr[x] = bitmask of y (directions, bit y-1) with E_xy known 0 */
static uint32_t Zr[33]; static int srcF[33][33], srcD[33][33], srcI[33][33], srcR[33][33]; static uint32_t srcS[33][33];
static int add_zero_block(uint32_t S, int F, int d, int i, int round) {
  int added = 0;
  for (int x = 1; x <= nd; x++) if (S >> (x - 1) & 1) {
    uint32_t nw = S & ~Zr[x];
    if (nw) {
      added = 1;
      for (int y = 1; y <= nd; y++) if (nw >> (y - 1) & 1) { srcF[x][y] = F; srcD[x][y] = d; srcI[x][y] = i; srcR[x][y] = round; srcS[x][y] = S; }
      Zr[x] |= S;
    }
  }
  return added;
}
/* maximum matching of rows R x cols C avoiding zero cells */
static int mt_match[33]; static uint32_t mt_seen;
static int aug(int x, uint32_t C) {
  uint32_t cand = C & ~Zr[x] & ~mt_seen;
  while (cand) { int y = __builtin_ctz(cand) + 1; cand &= cand - 1; mt_seen |= 1u << (y - 1);
    if (mt_match[y] == 0 || aug(mt_match[y], C)) { mt_match[y] = x; return 1; } }
  return 0;
}
static int maxmatch(uint32_t R, uint32_t C) {
  memset(mt_match, 0, sizeof mt_match); int r = 0;
  uint32_t T = R; while (T) { int x = __builtin_ctz(T) + 1; T &= T - 1; mt_seen = 0; r += aug(x, C); }
  return r;
}
/* maximum matching of X x X minus S x S (S inside X), avoiding zero cells */
static uint32_t ml_S; static uint32_t ml_X;
static int augL(int x) {
  uint32_t C = ml_X & ~Zr[x]; if (ml_S >> (x - 1) & 1) C &= ~ml_S;
  uint32_t cand = C & ~mt_seen;
  while (cand) { int y = __builtin_ctz(cand) + 1; cand &= cand - 1; mt_seen |= 1u << (y - 1);
    if (mt_match[y] == 0 || augL(mt_match[y])) { mt_match[y] = x; return 1; } }
  return 0;
}
static int maxmatchL(uint32_t X, uint32_t S) {
  memset(mt_match, 0, sizeof mt_match); int r = 0; ml_S = S; ml_X = X;
  uint32_t T = X; while (T) { int x = __builtin_ctz(T) + 1; T &= T - 1; mt_seen = 0; r += augL(x); }
  return r;
}
/* max clique among x with E_xx = 0, edges x~y if E_xy = E_yx = 0 */
static int bestq; static uint32_t bestY;
static void clique(uint32_t R, uint32_t P, int sz) {
  if (!P) { if (sz > bestq) { bestq = sz; bestY = R; } return; }
  if (sz + __builtin_popcount(P) <= bestq) return;
  while (P) {
    if (sz + __builtin_popcount(P) <= bestq) return;
    int v = __builtin_ctz(P) + 1; P &= P - 1;
    uint32_t nb = 0; for (int y = 1; y <= nd; y++) if ((Zr[v] >> (y - 1) & 1) && (Zr[y] >> (v - 1) & 1)) nb |= 1u << (y - 1);
    clique(R | (1u << (v - 1)), P & nb, sz + 1);
  }
}

static int lbmode = 1; static int dbgmode = 0; static long dbg_hist[200]; static int certmode = 0, rounds_max = 3, qsmax = -1, numax = 0; static int *cand, ncand; static uint32_t qS[40000]; static int qs[40000], qG[40000];
static int test(const int *cnt, char *buf) {
  int w[64], W = 0;
  for (int k = 0; k < nd; k++) W += cnt[k];
  for (int v = 0; v < nv; v++) { int s = 0; for (int k = 0; k < nd; k++) if (pk[v] >> k & 1) s += cnt[k]; w[v] = s; }
  for (int v = 0; v < nv; v++) if (minwv[v] < 1000 && minwv[v] < w[v]) return 1;
  uint8_t w8[64]; for (int v = 0; v < nv; v++) w8[v] = (uint8_t)w[v];   /* mw >= w (deletion) or mw = 255 */
  memset(Zr, 0, sizeof Zr);
  if (qsmax >= 0) { ncand = 0; for (int G = 1; G < nf; G++) if (wt[G] - w[vid[G]] <= qsmax) cand[ncand++] = G; }
  static int u[8192];
  for (int round = 1; round <= rounds_max; round++) {
    int newz = 0;
    for (int t = 0; t < nT; t++) {
      int F = Tf[t]; int nu = w[vid[F]]; if (nu > numax) continue;  /* rank E[D_F,D_F] = |D_F| - nu (T10) */
      uint8_t *mw = mwf + (long)t * NFL * nv; uint8_t *sz = szf + (long)t * NFL; uint32_t *sm = smf + (long)t * NFL;
      for (int i = 0; i < NFL; i++) {
        uint8_t bb = sz[i]; const uint8_t *row = mw + (long)i * nv;
        for (int v = 0; v < nv; v++) { uint8_t s = (uint8_t)(row[v] - w8[v]); bb = s < bb ? s : bb; }
        int b = bb;
        if (round > 1 && b > 0) { int mm = maxmatch(sm[i], sm[i]); if (mm < b) b = mm; }
        u[i] = b;
      }
      /* base zero blocks (tight forms or matching) */
      for (int d = 0; d < m; d++) for (int i = 0; i < nfl[d]; i++) if (u[off[d] + i] == 0 && sz[off[d] + i]) newz |= add_zero_block(sm[off[d] + i], F, d, i, round);
      for (int d = m - 2; d >= 0; d--) {
        int g = G_[d];
        for (int i = 0; i < nfl[d]; i++) {
          int n = sz[off[d] + i]; if (n == 0) continue;
          int tot = 0, mn = 1000;
          for (int j = 0; j < g; j++) { int x = u[off[d + 1] + up[d][i * g + j]]; tot += x; if (x < mn) mn = x; }
          int val = tot - 2 * n + 2 * nu;   /* STAR with nullity nu: (g-2) k <= sum u(B) - 2|S_A| + 2 nu */
          if (val < 0) { sprintf(buf, "Z %d %d %d %d %d %d S=%u", round, F, d, i, tot, 2 * n, sm[off[d] + i]); return 2; }
          int *ui = &u[off[d] + i]; int k = val / (g - 2);
          if (k < *ui) *ui = k;
          if (mn < *ui) *ui = mn;
          if (*ui == 0) newz |= add_zero_block(sm[off[d] + i], F, d, i, round);
        }
      }
      if (lbmode) {
        static int ell[8192];
        int ntop = __builtin_popcount(mask[F]) - nu;
        for (int i = 0; i < NFL; i++) ell[i] = 0;
        for (int i = 0; i < nfl[m - 1]; i++) { int n = sz[off[m - 1] + i]; int v = ntop - 2 * (__builtin_popcount(mask[F]) - n); if (v > ell[off[m - 1] + i]) ell[off[m - 1] + i] = v; }
        if (lbmode >= 2) {
          uint32_t DF = mask[F];
          for (int d = m - 1; d >= m - 2; d--) for (int i = 0; i < nfl[d]; i++) { int A = off[d] + i; if (!sz[A]) continue;
            int v = ntop - maxmatchL(DF, sm[A]); if (v > ell[A]) ell[A] = v; }
          /* odd forms G with D_G inside D_F */
          for (int G = 1; G < nf; G++) { if ((mask[G] & ~DF) || G == F) continue; int sG = wt[G] - w[vid[G]]; if (sG <= 0) continue;
            for (int i = 0; i < NFL; i++) { if (!sz[i]) continue; uint32_t I = mask[G] & sm[i]; int out = __builtin_popcount(mask[G] & ~sm[i]);
              if (out > 3 || sG - 2 * out <= ell[i] - 0 && out > 0 && sG - out <= ell[i]) continue;
              int v = sG - maxmatchL(mask[G], I); if (v > ell[i]) ell[i] = v; } }
        }
        for (int it = 0; it < 6; it++) {
          int changed = 0;
          for (int d = m - 2; d >= 0; d--) { int g = G_[d];
            for (int i = 0; i < nfl[d]; i++) { int n = sz[off[d] + i]; if (n == 0) continue; int *e = &ell[off[d] + i];
              for (int j = 0; j < g; j++) { int B = off[d + 1] + up[d][i * g + j]; int v = ell[B] - 2 * (sz[B] - n); if (v > *e) { *e = v; changed = 1; } } } }
          for (int d = m - 2; d >= 0; d--) { int g = G_[d];
            for (int i = 0; i < nfl[d]; i++) { int n = sz[off[d] + i]; if (n == 0) continue; int A = off[d] + i;
              int need = (g - 2) * ell[A] + 2 * n - 2 * nu; int tot = 0;
              for (int j = 0; j < g; j++) tot += u[off[d + 1] + up[d][i * g + j]];
              if (tot < need) { sprintf(buf, "L %d %d %d %d %d %d S=%u", round, F, d, i, tot, need, sm[A]); return 5; }
              for (int j = 0; j < g; j++) { int B = off[d + 1] + up[d][i * g + j]; int v = need - (tot - u[B]); if (ell[A] > v) v = ell[A];
                if (v > ell[B]) { ell[B] = v; changed = 1; } } } }
          if (dbgmode && it == 0) { int bg = -99, bi = -1; for (int i = 0; i < NFL; i++) if (sz[i] && ell[i] - u[i] > bg) { bg = ell[i] - u[i]; bi = i; } dbg_hist[bg + 99]++; if (bg >= -1) fprintf(stderr, "DBG round %d F %d gap %d flat %d |S| %d u %d ell %d\n", round, F, bg, bi, sz[bi], u[bi], ell[bi]); }
          for (int i = 0; i < NFL; i++) if (sz[i] && ell[i] > u[i]) { sprintf(buf, "L %d %d x %d %d %d %d S=%u", round, F, i, u[i], ell[i], 0, sm[i]); return 5; }
          if (!changed) break;
        }
      }
      if (qsmax >= 0) {
        /* quadratic star: Psi = D_F n D_G1 n D_G2; the three sets D_F n D_G1, D_F n D_G2, D_F n D_(F+G1+G2)
           (on D_F, F+G1+G2 = 1 exactly where G1 = G2) are Psi plus a partition of D_F \ Psi; STAR with g = 3. */
        uint32_t A = mask[F]; int nc = 0;
        for (int c = 0; c < ncand; c++) { int G = cand[c]; if (G == F) continue; uint32_t S = mask[G] & A; int sg = wt[G] - w[vid[G]];
          if (__builtin_popcount(S) > sg) { qS[nc] = S; qs[nc] = sg; qG[nc] = G; nc++; } }
        for (int a = 0; a < nc; a++) for (int b = a + 1; b < nc; b++) {
          uint32_t P = qS[a] & qS[b]; int n = __builtin_popcount(P); if (n == 0) continue;
          int G3 = F ^ qG[a] ^ qG[b]; uint32_t S3 = G3 ? (mask[G3] & A) : 0; int s3 = __builtin_popcount(S3);
          if (G3 && G3 != F) { int t = wt[G3] - w[vid[G3]]; if (t < s3) s3 = t; }
          int val = qs[a] + qs[b] + s3 - 2 * n + 2 * nu;
          if (val < 0) { sprintf(buf, "Q %d %d %d %d %d %d", round, F, qG[a], qG[b], qs[a] + qs[b] + s3, 2 * n); return 2; }
          if (val == 0) newz |= add_zero_block(P, F, -1, qG[a] * 65536 + qG[b], round);
        }
      }
    }
    /* clique test */
    uint32_t V = 0; for (int x = 1; x <= nd; x++) if (Zr[x] >> (x - 1) & 1) V |= 1u << (x - 1);
    bestq = 0; bestY = 0; clique(0, V, 0);
    if (bestq > W) {
      int p = sprintf(buf, "C %d %d %d Y=%u", round, bestq, W, bestY);
      if (certmode) {
        for (int x = 1; x <= nd; x++) if (bestY >> (x - 1) & 1) for (int y = 1; y <= nd; y++) if (bestY >> (y - 1) & 1)
          p += sprintf(buf + p, " %d,%d:%d/%u/%d/%d/%d", x, y, srcF[x][y], srcS[x][y], srcR[x][y], srcD[x][y], srcI[x][y]);
      }
      return 3;
    }
    /* E + K = I + Gamma is unit triangular, so rank E[Y,Y] >= |Y| - W for EVERY set Y of moves: test Y = all moves
       and Y = all moves but one */
    {
      uint32_t ALL = (nd == 32) ? 0xffffffffu : ((1u << nd) - 1);
      int mm = maxmatch(ALL, ALL);
      if (mm < nd - W) { sprintf(buf, "A %d %d %d", round, mm, nd - W); return 4; }
      for (int x = 1; x <= nd; x++) { uint32_t Y = ALL & ~(1u << (x - 1)); int m2 = maxmatch(Y, Y);
        if (m2 < nd - 1 - W) { sprintf(buf, "A %d %d %d x=%d", round, m2, nd - 1 - W, x); return 4; } }
    }
    /* matching test on every form */
    for (int Q = 1; Q < nf; Q++) {
      int sQ = wt[Q] - w[vid[Q]]; uint32_t D = mask[Q];
      int mm = maxmatch(D, D);
      if (mm < sQ) {
        int p = sprintf(buf, "M %d %d %d %d", round, Q, mm, sQ);
        if (certmode) for (int x = 1; x <= nd; x++) if (D >> (x - 1) & 1) for (int y = 1; y <= nd; y++) if ((D >> (y - 1) & 1) && (Zr[x] >> (y - 1) & 1))
          p += sprintf(buf + p, " %d,%d:%d/%u/%d/%d/%d", x, y, srcF[x][y], srcS[x][y], srcR[x][y], srcD[x][y], srcI[x][y]);
        return 4; }
    }
    if (!newz) break;
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
static inline uint32_t rotsub(uint32_t S, int e) { uint64_t x = S; return (uint32_t)(((x >> e) | (x << (31 - e))) & 0x7fffffffu); }
static int canonical(uint32_t S) {
  for (int j = 0; j < 5; j++) {
    uint32_t M = mulj(j, S), T = M;
    while (T) { int e = __builtin_ctz(T); T &= T - 1; if (rotsub(M, e) < S) return 0; }
  }
  return 1;
}

int main(int argc, char **argv) {
  load(argv[1]); mkflats();
  const char *md = getenv("MAXD"); const char *vm = getenv("VMASK"); const char *ro = getenv("ROUNDS"); const char *ce = getenv("CERT");
  uint64_t vmask = vm ? strtoull(vm, 0, 0) : 1ull;     /* parity classes of test forms (default: class 0 = even) */
  if (ro) rounds_max = atoi(ro); if (ce) certmode = atoi(ce);
  const char *qq = getenv("QSTAR"); if (qq) qsmax = atoi(qq); cand = malloc(4 * nf);
  const char *nn = getenv("NU"); if (nn) numax = atoi(nn);
  const char *lbe = getenv("LB"); if (lbe) lbmode = atoi(lbe);
  const char *dbe = getenv("DBG"); if (dbe) dbgmode = atoi(dbe);
  precompute(md ? atoi(md) : 20, (uint32_t)vmask, (uint32_t)(vmask >> 32));
  fprintf(stderr, "test forms %d\n", nT);
  static char buf[1 << 16]; int cnt[32]; long ntot = 0, nres[6] = {0};
  char line[1024];
  int enumW = 0, pp = 0, PP = 1; long eidx = 0; int er = 0; uint32_t ec = 0;
  const char *es = getenv("START"); long estart = es ? atol(es) : 0;
  const char *ee = getenv("END"); long eend = ee ? atol(ee) : (1L << 40);
  if (argc > 2 && !strcmp(argv[2], "enum")) { enumW = atoi(argv[3]); pp = atoi(argv[4]); PP = atoi(argv[5]); mktab(); er = enumW - 1; ec = (1u << er) - 1; }
  while (1) {
    memset(cnt, 0, sizeof cnt);
    if (enumW) {
      int got = 0;
      while (ec < (1u << 30)) {
        uint32_t S = 1u | (ec << 1); int take = 0;
        if (canonical(S)) { if (eidx >= eend) break; if (eidx % PP == pp && eidx >= estart) take = 1; eidx++; }
        uint32_t uu = ec & -ec, vv = ec + uu; ec = vv + (((vv ^ ec) / uu) >> 2);
        if (take) { for (int k = 0; k < 31; k++) if (S >> k & 1) cnt[k] = 1; got = 1; break; }
      }
      if (!got) break;
    } else {
      if (!fgets(line, sizeof line, stdin)) break;
      char *p = line; int any = 0;
      while (*p) { if (*p >= '0' && *p <= '9') { int k = strtol(p, &p, 10); cnt[k % nd]++; any = 1; } else p++; }
      if (!any) continue;
    }
    int res = test(cnt, buf);
    printf("S "); int first = 1;
    for (int k = 0; k < nd; k++) for (int r = 0; r < cnt[k]; r++) { printf(first ? "%d" : ",%d", k); first = 0; }
    if (res == 0) printf(" P\n"); else if (res == 1) printf(" D\n"); else printf(" %s\n", buf);
    fflush(stdout); ntot++; nres[res]++;
  }
  if (dbgmode) { fprintf(stderr, "gap histogram:"); for (int i = 0; i < 200; i++) if (dbg_hist[i]) fprintf(stderr, " %d:%ld", i - 99, dbg_hist[i]); fprintf(stderr, "\n"); }
  printf("SUMMARY tested %ld pass %ld deletion %ld star %ld clique %ld matching %ld lower %ld\n", ntot, nres[0], nres[1], nres[2], nres[3], nres[4], nres[5]);
  return 0;
}
