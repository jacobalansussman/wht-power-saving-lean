/* vr3.c (v-r3-classes, 16 Oct 2026): a skeptic's from-scratch re-implementation of r3-classes' stages 1-3.
   Shares no code with starrec*.c / bpm.c / vrec*.py.  Tables from mytab.py (own field, forms, parities).
   Facts used (each re-derived in ../NOTE.md):
    T10   rank E[D_Q,D_Q] = |D_Q| - r_Q.
    MONO  T in S  =>  rank E[T,T] <= rank E[S,S].
    FL    T in S  =>  rank E[S,S] <= rank E[T,T] + termrank(S x S minus T x T, avoiding known zero cells)   (Fact L)
          (with no zeros this contains the deletion bound, 2 per removed move).
    TR    rank <= term rank avoiding zero cells.
    WB    rank E[Y,Y] >= |Y| - W  (E = I + Gamma - K, rank K <= W).
    STAR  E[X,X] of nullity nu, X = P + Q_1..Q_g, k = rank E[P,P]: (g-2) k <= sum_j rank E[P u Q_j] - 2|P| + 2 nu.
    CLQ   E[Y,Y] = 0 => |Y| <= W.
   Per target list: for every test form F (r_F <= NU, |D_F| <= MAXW), bounds l(A) <= rank E[S_A,S_A] <= u(A) on the
   parts S_A = D_F n A of the affine flats A; iterate STAR / MONO / FL to a fixed point; l > u or a star sum below
   (g-2) l(A) is a contradiction.  Zero blocks (u = 0) are pooled; with POOL=1 every nontrivial u and every l >= 2 is
   pooled too (facts about the one matrix E).  After each round: CLQ and MAT (= TR on every D_Q).
   env: NU (0), POOL (0), MAXW (16), ROUNDS (3), SUBOUT (3), DUMP (0: print every final bound 'B mask l u').
   usage: vr3 tab.bin < lists ("k1,k2,..."); output: "<list> <verdict> <details>"  verdict P D L S C M            */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define PC __builtin_popcount
static int m, n, nd, nq; static uint32_t *qmask, *qpar; static int *qr, *qs;
static int W, NU, POOL, MAXW, ROUNDS, SUBOUT, DUMP;
static int nfl[7], *fstart; static uint32_t flats[4096]; static int fdim[4096], NF; static int up[4096][32], nup[4096];
static uint32_t zero[32];

static void mkflats(void) {
  /* linear subspaces by span closure (dedup on masks), then all cosets */
  static uint32_t lin[512]; static int ldim[512]; int nl = 1; lin[0] = 1u; ldim[0] = 0;
  for (int i = 0; i < nl; i++) for (int v = 1; v < n; v++) if (!(lin[i] >> v & 1)) {
    uint32_t L = lin[i], T = L; for (int p = 0; p < n; p++) if (L >> p & 1) T |= 1u << (p ^ v);
    int dup = 0; for (int j = 0; j < nl; j++) if (lin[j] == T) { dup = 1; break; }
    if (!dup) { lin[nl] = T; ldim[nl] = ldim[i] + 1; nl++; } }
  NF = 0;
  for (int d = 0; d <= m; d++) {
    int s0 = NF;
    for (int i = 0; i < nl; i++) if (ldim[i] == d) for (int a = 0; a < n; a++) {
      uint32_t T = 0; for (int p = 0; p < n; p++) if (lin[i] >> p & 1) T |= 1u << (p ^ a);
      int dup = 0; for (int j = s0; j < NF; j++) if (flats[j] == T) { dup = 1; break; }
      if (!dup) { flats[NF] = T; fdim[NF] = d; NF++; } }
    nfl[d] = NF - s0;
  }
  for (int i = 0; i < NF; i++) {
    nup[i] = 0;
    for (int j = 0; j < NF; j++) if (fdim[j] == fdim[i] + 1 && (flats[j] & flats[i]) == flats[i]) up[i][nup[i]++] = j;
    if (fdim[i] < m && nup[i] != (1 << (m - fdim[i])) - 1) { fprintf(stderr, "up count wrong\n"); exit(1); }
  }
  if (getenv("FLATINFO")) { for (int d = 0; d <= m; d++) fprintf(stderr, "dim %d: %d flats\n", d, nfl[d]); }
}

/* term rank of rows x cols minus ex x ex, avoiding zero cells */
static uint32_t gadj[32]; static int matchc[32]; static uint32_t vis;
static int aug(int i) {
  uint32_t a = gadj[i] & ~vis;
  while (a) { int j = __builtin_ctz(a); a &= a - 1; vis |= 1u << j;
    if (matchc[j] < 0 || aug(matchc[j])) { matchc[j] = i; return 1; } }
  return 0;
}
static int termrank(uint32_t rows, uint32_t cols, uint32_t ex) {
  for (int j = 0; j < 32; j++) matchc[j] = -1;
  int r = 0; uint32_t R = rows;
  while (R) { int i = __builtin_ctz(R); R &= R - 1;
    gadj[i] = cols & ~zero[i]; if (ex >> i & 1) gadj[i] &= ~ex; }
  R = rows;
  while (R) { int i = __builtin_ctz(R); R &= R - 1; vis = 0; r += aug(i); }
  return r;
}

/* slack bound with a hash cache: min s_G over forms G with S inside D_G */
#define HS (1 << 22)
static uint32_t *hkey; static int8_t *hval; static uint16_t **plist; static int *plen; /* forms containing pair */
static int pidx(int a, int b) { return a * 32 + b; }
static int uslack(uint32_t S) {
  if (!S) return 0;
  uint32_t h = (S * 2654435761u) >> 10; h &= HS - 1;
  while (hkey[h]) { if (hkey[h] == S) return hval[h]; h = (h + 1) & (HS - 1); }
  int best = PC(S);
  int a = __builtin_ctz(S); uint32_t S2 = S & (S - 1); int b = S2 ? __builtin_ctz(S2) : a;
  uint16_t *L = plist[pidx(a, b)]; int len = plen[pidx(a, b)];
  for (int t = 0; t < len; t++) { int q = L[t]; if ((qmask[q] & S) == S && qs[q] < best) best = qs[q]; }
  hkey[h] = S; hval[h] = best; return best;
}

/* pooled bounds */
typedef struct { uint32_t S; int v; } PB;
static PB *pu, *pl; static int npu, npl, cpu_, cpl_;
static int poolu(uint32_t S, int cur) { for (int i = 0; i < npu; i++) if (pu[i].v < cur && (pu[i].S & S) == S) cur = pu[i].v; return cur; }
static int pooll(uint32_t S, int cur) { for (int i = 0; i < npl; i++) if (pl[i].v > cur && (pl[i].S & S) == pl[i].S) cur = pl[i].v; return cur; }
static void addpu(uint32_t S, int v) { for (int i = 0; i < npu; i++) if (pu[i].S == S) { if (v < pu[i].v) pu[i].v = v; return; }
  if (npu == cpu_) { cpu_ = cpu_ ? 2 * cpu_ : 1024; pu = realloc(pu, cpu_ * sizeof(PB)); } pu[npu].S = S; pu[npu].v = v; npu++; }
static void addpl(uint32_t S, int v) { for (int i = 0; i < npl; i++) if (pl[i].S == S) { if (v > pl[i].v) pl[i].v = v; return; }
  if (npl == cpl_) { cpl_ = cpl_ ? 2 * cpl_ : 1024; pl = realloc(pl, cpl_ * sizeof(PB)); } pl[npl].S = S; pl[npl].v = v; npl++; }

static int lo[4096], hi[4096], su[4096]; static uint32_t SA[4096]; static int fl_pair_tr[4096][32][2];
static char why[256];

/* returns 0 = fine, 1 = contradiction (why filled) */
static int runform(int F, int *newfacts) {
  uint32_t X = qmask[F]; int nu = qr[F]; int top = PC(X) - nu;
  /* sub-forms of F */
  static int subs[40000]; int ns = 0;
  for (int q = 1; q < nq; q++) if ((qmask[q] & ~X) == 0) subs[ns++] = q;
  for (int i = 0; i < NF; i++) {
    uint32_t S = X & flats[i]; SA[i] = S; su[i] = PC(S);
    if (fdim[i] == m) { lo[i] = hi[i] = top; continue; }
    if (!S) { lo[i] = hi[i] = 0; continue; }
    int u = su[i]; int t = uslack(S); if (t < u) u = t; t = termrank(S, S, 0); if (t < u) u = t;
    if (POOL) u = poolu(S, u);
    int l = su[i] - W; if (l < 0) l = 0;
    t = top - termrank(X, X, S); if (t > l) l = t;            /* FL with X = D_F */
    for (int k = 0; k < ns; k++) { uint32_t D = qmask[subs[k]]; if (PC(D & ~S) > SUBOUT || !(D & S)) continue;
      t = qs[subs[k]] - termrank(D, D, D & S); if (t > l) l = t; }   /* FL with X = D_G, then MONO */
    if (POOL) l = pooll(S, l);
    lo[i] = l; hi[i] = u;
  }
  for (int i = 0; i < NF; i++) for (int j = 0; j < nup[i]; j++) {
    int B = up[i][j]; fl_pair_tr[i][j][0] = termrank(SA[B], SA[B], SA[i]); }
  int start_u[4096]; memcpy(start_u, hi, sizeof(int) * NF);
  for (int it = 0; it < 50; it++) {
    int ch = 0;
    for (int d = m - 2; d >= 0; d--) for (int i = 0; i < NF; i++) if (fdim[i] == d) {
      int g = nup[i], tot = 0; for (int j = 0; j < g; j++) tot += hi[up[i][j]];
      int rhs = tot - 2 * su[i] + 2 * nu;
      if (rhs < (g - 2) * lo[i]) { sprintf(why, "S F=%d dim=%d S=%u tot=%d need=%d", F, d, SA[i], tot, (g - 2) * lo[i] + 2 * su[i] - 2 * nu); return 1; }
      int nu_ = rhs / (g - 2); if (nu_ < hi[i]) { hi[i] = nu_; ch = 1; }
      for (int j = 0; j < g; j++) { int B = up[i][j]; int lb = (g - 2) * lo[i] + 2 * su[i] - 2 * nu - (tot - hi[B]);
        if (lb > lo[B]) { lo[B] = lb; ch = 1; } }
    }
    for (int i = 0; i < NF; i++) for (int j = 0; j < nup[i]; j++) {
      int B = up[i][j], t = fl_pair_tr[i][j][0];
      if (hi[B] < hi[i]) { hi[i] = hi[B]; ch = 1; }                 /* MONO */
      if (lo[i] > lo[B]) { lo[B] = lo[i]; ch = 1; }                 /* MONO */
      if (lo[B] - t > lo[i]) { lo[i] = lo[B] - t; ch = 1; }         /* FL inside B */
      if (hi[i] + t < hi[B]) { hi[B] = hi[i] + t; ch = 1; }         /* FL inside B */
    }
    for (int i = 0; i < NF; i++) if (lo[i] > hi[i]) { sprintf(why, "L F=%d dim=%d S=%u l=%d u=%d", F, fdim[i], SA[i], lo[i], hi[i]); return 1; }
    if (!ch) break;
  }
  for (int i = 0; i < NF; i++) {
    if (fdim[i] == m || !SA[i]) continue;
    if (DUMP) printf("B %u %d %d\n", SA[i], lo[i], hi[i]);
    if (hi[i] == 0) { uint32_t S = SA[i]; uint32_t R = S;
      while (R) { int z = __builtin_ctz(R); R &= R - 1; if ((zero[z] & S) != S) { zero[z] |= S; *newfacts = 1; } } }
    if (POOL) {
      if (hi[i] < start_u[i] && hi[i] > 0) { int before = npu; addpu(SA[i], hi[i]); if (npu > before) *newfacts = 1; }
      if (lo[i] >= 2 && lo[i] > su[i] - W) { int before = npl; addpl(SA[i], lo[i]); if (npl > before) *newfacts = 1; }
    }
  }
  return 0;
}

static int maxclq(uint32_t R, uint32_t P, int sz, int best, uint32_t *adj) {
  if (!P) return sz > best ? sz : best;
  if (sz + PC(P) <= best) return best;
  while (P) { if (sz + PC(P) <= best) break; int v = __builtin_ctz(P); P &= P - 1;
    int b = maxclq(R | (1u << v), P & adj[v], sz + 1, best, adj); if (b > best) best = b; }
  return best;
}

static int *order; static int norder;
static int cmpw(const void *a, const void *b) { int x = *(int *)a, y = *(int *)b; int wx = PC(qmask[x]), wy = PC(qmask[y]);
  return wx != wy ? wx - wy : x - y; }

int main(int argc, char **argv) {
  NU = getenv("NU") ? atoi(getenv("NU")) : 0; POOL = getenv("POOL") ? atoi(getenv("POOL")) : 0;
  MAXW = getenv("MAXW") ? atoi(getenv("MAXW")) : 16; ROUNDS = getenv("ROUNDS") ? atoi(getenv("ROUNDS")) : 3;
  SUBOUT = getenv("SUBOUT") ? atoi(getenv("SUBOUT")) : 3; DUMP = getenv("DUMP") ? atoi(getenv("DUMP")) : 0;
  FILE *f = fopen(argv[1], "rb"); uint32_t hd[2]; if (fread(hd, 4, 2, f) != 2) return 1; m = hd[0]; nq = hd[1];
  n = 1 << m; nd = n - 1;
  qmask = malloc(4 * nq); qpar = malloc(4 * nq); qr = malloc(4 * nq); qs = malloc(4 * nq);
  for (int q = 0; q < nq; q++) { uint32_t t[2]; if (fread(t, 4, 2, f) != 2) return 1; qmask[q] = t[0]; qpar[q] = t[1]; }
  fclose(f);
  mkflats();
  hkey = calloc(HS, 4); hval = malloc(HS);
  plist = calloc(32 * 32, sizeof(uint16_t *)); plen = calloc(32 * 32, sizeof(int));
  for (int a = 1; a < n; a++) for (int b = a; b < n; b++) {
    int c = 0; for (int q = 0; q < nq; q++) if ((qmask[q] >> a & 1) && (qmask[q] >> b & 1)) c++;
    plist[pidx(a, b)] = malloc(2 * c + 2); plen[pidx(a, b)] = 0;
    for (int q = 0; q < nq; q++) if ((qmask[q] >> a & 1) && (qmask[q] >> b & 1)) plist[pidx(a, b)][plen[pidx(a, b)]++] = q; }
  order = malloc(4 * nq);
  char line[512];
  while (fgets(line, sizeof line, stdin)) {
    int ks[64]; int nk = 0; char *p = line; line[strcspn(line, "\r\n")] = 0; if (!*line) continue;
    char copy[512]; strcpy(copy, line);
    while (*p) { ks[nk++] = strtol(p, &p, 10); if (*p == ',') p++; else break; }
    W = nk;
    int del = 0;
    for (int q = 0; q < nq; q++) { int r = 0; for (int l = 0; l < nk; l++) r += qpar[q] >> ks[l] & 1; qr[q] = r; qs[q] = PC(qmask[q]) - r; if (qs[q] < 0) del = 1; }
    if (del) { printf("%s D\n", copy); fflush(stdout); continue; }
    memset(hkey, 0, 4 * (size_t)HS); memset(zero, 0, sizeof zero); npu = npl = 0;
    norder = 0; for (int q = 1; q < nq; q++) if (qr[q] <= NU && PC(qmask[q]) <= MAXW) order[norder++] = q;
    qsort(order, norder, 4, cmpw);
    int done = 0;
    for (int rd = 1; rd <= ROUNDS && !done; rd++) {
      int newf = 0;
      for (int t = 0; t < norder && !done; t++) if (runform(order[t], &newf)) { printf("%s %s round=%d\n", copy, why, rd); done = 1; }
      if (done) break;
      uint32_t adj[32]; uint32_t V = 0;
      for (int z = 1; z < n; z++) if (zero[z] >> z & 1) V |= 1u << z;
      for (int z = 1; z < n; z++) { adj[z] = 0; if (V >> z & 1) for (int y = 1; y < n; y++) if (y != z && (V >> y & 1) && (zero[z] >> y & 1) && (zero[y] >> z & 1)) adj[z] |= 1u << y; }
      int cq = maxclq(0, V, 0, 0, adj);
      if (cq > W) { printf("%s C clique=%d round=%d\n", copy, cq, rd); done = 1; break; }
      for (int q = 1; q < nq; q++) { int t = termrank(qmask[q], qmask[q], 0); if (t < qs[q]) { printf("%s M Q=%d tr=%d s=%d round=%d\n", copy, q, t, qs[q], rd); done = 1; break; } }
      if (done) break;
      if (!newf) break;
    }
    if (!done) { int zc = 0; for (int z = 1; z < n; z++) zc += PC(zero[z]); printf("%s P zeros=%d pu=%d pl=%d\n", copy, zc, npu, npl); }
    fflush(stdout);
  }
  return 0;
}
