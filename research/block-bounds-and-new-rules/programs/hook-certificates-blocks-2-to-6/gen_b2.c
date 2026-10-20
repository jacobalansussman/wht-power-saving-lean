/* gen_b2.c (pf-logging, 14 Oct 2026): GENERATOR (untrusted) for f(dir z1 dir z2) <= dim N - 16 r at f odd.  Builds the exact
 * vectors c of Y(z1) cap Y(z2) (as in check_yint.c), the spanning set {D1 c, D2 c} of U = D1 Z12 + D2 Z12 (inside
 * Y(B2) = {v in N : D1 D2 v in N}), projects one representative per X_s-orbit to the trivial isotypic block (alpha = 0,
 * coordinates (0, R) with 0 in R), and finds r = rank mod p with a nonsingular minor.  Writes the minor as
 * "rows: member-index j" and "cols: R-mask".  Usage: ./gen_b2 m k z1 z2 OUT
 * (Header of the original check_yint.c follows.)
 * Exact vectors in the intersection  Y(z_1) cap ... cap Y(z_r),  Y(z) = {v in N : dir(z)_W v in N},
 * for the hook N in W = V (x) Lambda^k V (PROOF.md section 2b).
 *
 * Construction (Lemma A'): Y(z) is spanned by Lemma A's vectors u^z_{a,w}; on a two-term support {p1 = (a, a+w), p2 = (b, b+w)}
 * (b = a xor z, a, b not in w) the coordinates of any v in Y(z) satisfy v[p2] = -s_a s_b v[p1], where s_x = +-1 is the sign of
 * e_x ^ e_w.  Put an edge p1 -- p2 with that ratio for every two-term pair of every z in the list; a connected component whose
 * ratios agree around every cycle carries exactly one vector c (entries +-1, supported on the component) lying in every Y(z_j),
 * and different components have disjoint supports.
 * This program builds those vectors and CHECKS EACH ONE EXACTLY (Gaussian integers, no modulus): c in N, and
 * (2 dir z_j)_W c in N for every j.  Output: the number of consistent components = a certified lower bound for
 * dim(Y(z_1) cap ... cap Y(z_r)) (equality by Lemma A').
 * Build: cc -O2 -o check_yint check_yint.c      Run: ./check_yint m k z1 [z2 ...]      (n = 2^m <= 16)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef int64_t i64; typedef uint64_t u64;
typedef struct { i64 re, im; } gi;
static gi gmul(gi a, gi b) { gi r = { a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re }; return r; }
static int n, k;
#define HS (1 << 16)
static u64 hkey[HS]; static int hy[HS]; static gi hval[HS]; static int hused[HS], nused; static char hocc[HS];
static void hadd(int y, u64 key, gi v) {
    u64 h = ((key ^ ((u64)y << 58) ^ (u64)y) * 0x9E3779B97F4A7C15ULL) >> 48;
    while (hocc[h] && (hkey[h] != key || hy[h] != y)) h = (h + 1) & (HS - 1);
    if (!hocc[h]) { hocc[h] = 1; hkey[h] = key; hy[h] = y; hval[h].re = 0; hval[h].im = 0; hused[nused++] = (int)h; if (nused > HS / 2) { fprintf(stderr, "hash full\n"); exit(3); } }
    hval[h].re += v.re; hval[h].im += v.im;
}
static void hclear(void) { for (int i = 0; i < nused; i++) hocc[hused[i]] = 0; nused = 0; }
static int wedge_sign(const int *c, int len, u64 *mask) {
    u64 m = 0; int inv = 0;
    for (int i = 0; i < len; i++) { if (m >> c[i] & 1) return 0; m |= 1ULL << c[i]; }
    for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) if (c[i] > c[j]) inv++;
    *mask = m; return (inv & 1) ? -1 : 1;
}
static void add_image(int x, u64 T, gi coef, int z) {   /* coef * (2 dir z)_W (e_x (x) e_T) */
    int t[64], L = 0; for (int j = 0; j < n; j++) if (T >> j & 1) t[L++] = j;
    gi al = { 1, 1 }, be = { 1, -1 };
    for (int cx = 0; cx < 2; cx++) {
        int y = cx ? (x ^ z) : x; gi cy = gmul(coef, cx ? be : al);
        for (u64 cm = 0; cm < (1ULL << L); cm++) {
            int c[64]; gi cc = cy;
            for (int j = 0; j < L; j++) { int bit = (int)(cm >> j & 1); c[j] = bit ? (t[j] ^ z) : t[j]; cc = gmul(cc, bit ? be : al); }
            u64 M; int s = wedge_sign(c, L, &M); if (!s) continue;
            if (s < 0) { cc.re = -cc.re; cc.im = -cc.im; }
            hadd(y, M, cc);
        }
    }
}
/* basis of N: index i <-> (x, T), x in T */
static int *bx; static u64 *bT; static int nb; static int *idx_of;   /* idx_of[x * 2^n + T] */
static int find_idx(int x, u64 T) { return idx_of[((size_t)x << n) | T]; }
/* union-find with signs: par[i], sg[i] = sign of i relative to par[i] (value[i] = sg[i] * value[par[i]]) */
static int *par, *sg, *bad;
static int findr(int i, int *s) { int acc = 1; while (par[i] != i) { acc *= sg[i]; i = par[i]; } *s = acc; return i; }
static void compress(int i) { int s; int r = findr(i, &s); while (par[i] != r) { int nx = par[i]; int s2; findr(i, &s2); sg[i] = s2; par[i] = r; i = nx; } }

#define P 2147483629ULL
static u64 pw(u64 b, u64 e) { u64 r = 1; b %= P; while (e) { if (e & 1) r = r * b % P; b = b * b % P; e >>= 1; } return r; }
int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: gen_b2 m k z1 z2 OUT\n"); return 2; }
    int m = atoi(argv[1]); k = atoi(argv[2]); n = 1 << m; int zs[2] = { atoi(argv[3]), atoi(argv[4]) }, nz = 2;
    u64 I = 0; for (u64 g = 2; g < 1000; g++) if (pw(g, (P - 1) / 2) == P - 1) { I = pw(g, (P - 1) / 4); break; }
    idx_of = malloc(sizeof(int) * ((size_t)n << n)); for (size_t i = 0; i < ((size_t)n << n); i++) idx_of[i] = -1;
    bx = malloc(sizeof(int) * ((size_t)n << n)); bT = malloc(sizeof(u64) * ((size_t)n << n)); nb = 0;
    for (u64 T = 0; T < (1ULL << n); T++) if (__builtin_popcountll(T) == k) for (int x = 0; x < n; x++) if (T >> x & 1) { idx_of[((size_t)x << n) | T] = nb; bx[nb] = x; bT[nb] = T; nb++; }
    int *colof = malloc(sizeof(int) << n); int ncol = 0; u64 *colmask = malloc(sizeof(u64) << n);
    for (u64 T = 0; T < (1ULL << n); T++) { colof[T] = -1; if (__builtin_popcountll(T) == k && (T & 1)) { colof[T] = ncol; colmask[ncol++] = T; } }
    par = malloc(sizeof(int) * nb); sg = malloc(sizeof(int) * nb); bad = calloc(nb, sizeof(int));
    for (int i = 0; i < nb; i++) { par[i] = i; sg[i] = 1; }
    for (int iz = 0; iz < nz; iz++) { int z = zs[iz];
        for (int i = 0; i < nb; i++) {
            int a = bx[i], b = a ^ z; u64 T = bT[i]; if (b < a || (T >> b & 1)) continue;
            u64 w = T & ~(1ULL << a); int c1[64], c2[64], L = 0; c1[L] = a; c2[L] = b; L++;
            for (int j = 0; j < n; j++) if (w >> j & 1) { c1[L] = j; c2[L] = j; L++; }
            u64 M1, M2; int sa = wedge_sign(c1, L, &M1), sb = wedge_sign(c2, L, &M2);
            int j2 = find_idx(b, M2); int ratio = -sa * sb; int s1, s2; int r1 = findr(i, &s1), r2 = findr(j2, &s2);
            if (r1 == r2) { if (s2 != ratio * s1) bad[r1] = 1; } else { par[r2] = r1; sg[r2] = ratio * s1 * s2; bad[r1] |= bad[r2]; }
            compress(i); compress(j2);
        } }
    int *root = malloc(sizeof(int) * nb), *rsg = malloc(sizeof(int) * nb);
    for (int i = 0; i < nb; i++) { int s; root[i] = findr(i, &s); rsg[i] = s; }
    int *head = malloc(sizeof(int) * nb), *next = malloc(sizeof(int) * nb), *minm = malloc(sizeof(int) * nb);
    for (int i = 0; i < nb; i++) { head[i] = -1; minm[i] = nb; }
    for (int i = nb - 1; i >= 0; i--) { next[i] = head[root[i]]; head[root[i]] = i; if (i < minm[root[i]]) minm[root[i]] = i; }
    /* representatives: components whose minimum member index is minimal over their X_s-images */
    int nrep = 0; int *reps = malloc(sizeof(int) * nb);
    for (int r = 0; r < nb; r++) if (root[r] == r && !bad[r]) {
        int isrep = 1;
        for (int s = 1; s < n && isrep; s++) {
            int mn = nb; for (int i = head[r]; i >= 0; i = next[i]) { u64 T2 = 0; for (int t = 0; t < n; t++) if (bT[i] >> t & 1) T2 |= 1ULL << (t ^ s); int j = find_idx(bx[i] ^ s, T2); if (j < mn) mn = j; }
            if (mn < minm[r]) isrep = 0;
        }
        if (isrep) reps[nrep++] = r;
    }
    fprintf(stderr, "components: reps %d, columns %d\n", nrep, ncol);
    int nrows = 2 * nrep; uint32_t *A = malloc(sizeof(uint32_t) * (size_t)nrows * ncol);
    for (int q = 0; q < nrep; q++) for (int j = 0; j < 2; j++) {
        int r = reps[q]; uint32_t *row = A + (size_t)(2 * q + j) * ncol; memset(row, 0, sizeof(uint32_t) * ncol);
        hclear(); for (int i = head[r]; i >= 0; i = next[i]) { gi cf = { rsg[i], 0 }; add_image(bx[i], bT[i], cf, zs[j]); }
        for (int h = 0; h < nused; h++) {
            gi v = hval[hused[h]]; if (!v.re && !v.im) continue; int y = hy[hused[h]]; u64 T = hkey[hused[h]];
            if (!(T >> y & 1)) { fprintf(stderr, "image outside N!\n"); return 3; }
            int c[64], L = 0; for (int t = 0; t < n; t++) if (T >> t & 1) c[L++] = t ^ y;
            u64 R; int s = wedge_sign(c, L, &R); int col = colof[R];
            u64 val = ((u64)((v.re % (i64)P + (i64)P) % (i64)P) + (u64)((v.im % (i64)P + (i64)P) % (i64)P) * I) % P;
            if (s < 0) val = (P - val) % P;
            row[col] = (uint32_t)((row[col] + val) % P);
        }
    }
    fprintf(stderr, "rows built: %d\n", nrows);
    /* incremental elimination: basis rows reduced, pivot = first nonzero column */
    u64 *B = malloc(sizeof(u64) * (size_t)ncol * ncol); int *pc = malloc(sizeof(int) * ncol), *prow = malloc(sizeof(int) * ncol); int rank = 0;
    u64 *t = malloc(sizeof(u64) * ncol);
    for (int i = 0; i < nrows; i++) {
        for (int c = 0; c < ncol; c++) t[c] = A[(size_t)i * ncol + c];
        for (int b = 0; b < rank; b++) { u64 f = t[pc[b]]; if (!f) continue; u64 *br = B + (size_t)b * ncol; u64 g = P - f;
            for (int c = 0; c < ncol; c++) if (br[c]) t[c] = (t[c] + g * br[c]) % P; }
        int piv = -1; for (int c = 0; c < ncol; c++) if (t[c]) { piv = c; break; }
        if (piv < 0) continue;
        u64 iv = pw(t[piv], P - 2); u64 *br = B + (size_t)rank * ncol; for (int c = 0; c < ncol; c++) br[c] = t[c] * iv % P;
        pc[rank] = piv; prow[rank] = i; rank++;
        if (rank % 500 == 0) fprintf(stderr, "  rank %d after %d rows\n", rank, i + 1);
    }
    printf("m=%d k=%d z=%d,%d: reps %d, rows %d, rank mod p of the alpha=0 projection of U = D1 Z12 + D2 Z12: %d -> dim U >= %d x %d = %d (via Z-symmetry), f(B2) <= %d - %d = %d\n",
           m, k, zs[0], zs[1], nrep, nrows, rank, n, rank, n * rank, nb, n * rank, nb - n * rank);
    FILE *fo = fopen(argv[5], "w");
    fprintf(fo, "m %d k %d z1 %d z2 %d p %llu I %llu size %d\n", m, k, zs[0], zs[1], (unsigned long long)P, (unsigned long long)I, rank);
    for (int b = 0; b < rank; b++) { int q = prow[b] / 2, j = prow[b] % 2; fprintf(fo, "%d %llu %d %llu\n", bx[minm[reps[q]]], (unsigned long long)bT[minm[reps[q]]], j, (unsigned long long)colmask[pc[b]]); }
    fclose(fo);
    return 0;
}
