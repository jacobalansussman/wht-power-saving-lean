/* pf-logging, 14 Oct 2026.  Exact vectors in the intersection  Y(z_1) cap ... cap Y(z_r),  Y(z) = {v in N : dir(z)_W v in N},
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
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: check_yint m k z1 [z2 ...]\n"); return 2; }
    int m = atoi(argv[1]); k = atoi(argv[2]); n = 1 << m; if (m < 1 || m > 4) { fprintf(stderr, "m <= 4\n"); return 2; }
    int zs[16], nz = 0; for (int i = 3; i < argc; i++) zs[nz++] = atoi(argv[i]);
    idx_of = malloc(sizeof(int) * ((size_t)n << n)); for (size_t i = 0; i < ((size_t)n << n); i++) idx_of[i] = -1;
    bx = malloc(sizeof(int) * ((size_t)n << n)); bT = malloc(sizeof(u64) * ((size_t)n << n)); nb = 0;
    for (u64 T = 0; T < (1ULL << n); T++) if (__builtin_popcountll(T) == k) for (int x = 0; x < n; x++) if (T >> x & 1) { idx_of[((size_t)x << n) | T] = nb; bx[nb] = x; bT[nb] = T; nb++; }
    par = malloc(sizeof(int) * nb); sg = malloc(sizeof(int) * nb); bad = calloc(nb, sizeof(int));
    for (int i = 0; i < nb; i++) { par[i] = i; sg[i] = 1; }
    long edges = 0;
    for (int iz = 0; iz < nz; iz++) {
        int z = zs[iz];
        for (int i = 0; i < nb; i++) {
            int a = bx[i], b = a ^ z; u64 T = bT[i]; if (b < a || (T >> b & 1)) continue;   /* two-term pairs, counted once */
            u64 w = T & ~(1ULL << a); int c1[64], c2[64], L = 0; c1[L] = a; c2[L] = b; L++;
            for (int j = 0; j < n; j++) if (w >> j & 1) { c1[L] = j; c2[L] = j; L++; }
            u64 M1, M2; int sa = wedge_sign(c1, L, &M1), sb = wedge_sign(c2, L, &M2);
            int j2 = find_idx(b, M2); if (M1 != T || j2 < 0) { fprintf(stderr, "internal\n"); return 3; }
            int ratio = -sa * sb;   /* u = sa e_p1 - sb e_p2 : v[p2] / v[p1] = -sb / sa */
            edges++;
            int s1, s2; int r1 = findr(i, &s1), r2 = findr(j2, &s2);
            /* value[j2] must equal ratio * value[i]:  s2 v[r2] = ratio s1 v[r1] */
            if (r1 == r2) { if (s2 != ratio * s1) bad[r1] = 1; }
            else { par[r2] = r1; sg[r2] = ratio * s1 * s2; bad[r1] |= bad[r2]; }
            compress(i); compress(j2);
        }
    }
    /* components */
    int ncomp = 0, nbad = 0; int *root = malloc(sizeof(int) * nb), *rsg = malloc(sizeof(int) * nb);
    for (int i = 0; i < nb; i++) { int s; root[i] = findr(i, &s); rsg[i] = s; }
    for (int i = 0; i < nb; i++) if (root[i] == i) { ncomp++; if (bad[i]) nbad++; }
    /* group members by root */
    int *head = malloc(sizeof(int) * nb), *next = malloc(sizeof(int) * nb); for (int i = 0; i < nb; i++) head[i] = -1;
    for (int i = nb - 1; i >= 0; i--) { next[i] = head[root[i]]; head[root[i]] = i; }
    long checked = 0, fails = 0, maxsize = 0;
    for (int r = 0; r < nb; r++) if (root[r] == r && !bad[r]) {
        long sz = 0; for (int i = head[r]; i >= 0; i = next[i]) sz++; if (sz > maxsize) maxsize = sz;
        for (int iz = 0; iz < nz; iz++) {
            hclear();
            for (int i = head[r]; i >= 0; i = next[i]) { gi cf = { rsg[i], 0 }; if (!(bT[i] >> bx[i] & 1)) fails++; add_image(bx[i], bT[i], cf, zs[iz]); }
            for (int h = 0; h < nused; h++) { gi v = hval[hused[h]]; if (!v.re && !v.im) continue; if (!(hkey[hused[h]] >> hy[hused[h]] & 1)) { fails++; break; } }
        }
        checked++;
    }
    printf("m=%d k=%d directions", m, k); for (int i = 0; i < nz; i++) printf(" %d", zs[i]);
    printf(": dim N = %d, two-term edges %ld, components %d (sign-inconsistent %d), consistent vectors checked %ld (max support %ld), exact failures %ld\n",
           nb, edges, ncomp, nbad, checked, maxsize, fails);
    printf("RESULT: %s -> dim( cap_j Y(z_j) ) >= %ld   (and = by Lemma A')\n", fails ? "FAIL" : "PASS", checked);
    return fails ? 1 : 0;
}
