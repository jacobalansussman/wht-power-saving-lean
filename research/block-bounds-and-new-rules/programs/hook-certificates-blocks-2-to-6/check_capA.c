/* pf-logging, 14 Oct 2026.  EXACT check of Lemma A's kernel (PROOF.md section 2) for the hook N in W = V (x) Lambda^k V.
 *
 * For every direction z != 0 of F_2^m (or the ones given on the command line) it builds, from the closed form,
 *   u_{a,w} = e_a (x) (e_a ^ e_w) - e_b (x) (e_b ^ e_w),   b = a xor z,  {a,b} an unordered pair, w a (k-1)-subset,
 * and checks, in exact Gaussian-integer arithmetic (no modulus anywhere):
 *   (1) every u lies in N = span{e_x (x) e_T : x in T};
 *   (2) G_W u lies in N, where G = 2 dir z = (1+i) Id + (1-i) X_z and G_W = G (x) Lambda^k G
 *       (computed by expanding G e_x (x) G e_{t_1} ^ ... ^ G e_{t_k} term by term, sorting each wedge with its sign);
 *   (3) the nonzero u have pairwise disjoint supports (so they are linearly independent).
 * Then {v in N : dir_W v in N} has dimension >= #nonzero u, and by Lemma B' f(dir z) <= dim N - #nonzero u.
 * dim N is counted directly (pairs (x,T) with x in T).  Shares no code with gen_m3.py / check_m3.py / gen_m0.c.
 * Build: cc -O2 -o check_capA check_capA.c      Run: ./check_capA m k [z ...]       (n = 2^m <= 64, k <= n)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef int64_t i64; typedef uint64_t u64;
typedef struct { i64 re, im; } gi;
static gi gmul(gi a, gi b) { gi r = { a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re }; return r; }
static int n, k;
static u64 *SUB; static size_t nsub, csub;
static void subs(int lo, int need, u64 cur) {   /* all need-subsets of [lo, n) added to cur */
    if (!need) { if (nsub == csub) { csub = csub ? 2 * csub : 1024; SUB = realloc(SUB, csub * sizeof(u64)); } SUB[nsub++] = cur; return; }
    for (int x = lo; x <= n - need; x++) subs(x + 1, need - 1, cur | (1ULL << x));
}
typedef struct { u64 T; int x; } sk;
static int skcmp(const void *a, const void *b) { const sk *p = a, *q = b; if (p->T != q->T) return p->T < q->T ? -1 : 1; return (p->x > q->x) - (p->x < q->x); }
/* hash table for one image vector: key = (y << 32) | Tmask, value = Gaussian integer */
#define HS (1 << 14)
static u64 hkey[HS]; static int hy[HS]; static gi hval[HS]; static int hused[HS], nused; static char hocc[HS];
static void hadd(int y, u64 key, gi v) {
    u64 h = ((key ^ ((u64)y << 58) ^ (u64)y) * 0x9E3779B97F4A7C15ULL) >> 50;
    while (hocc[h] && (hkey[h] != key || hy[h] != y)) h = (h + 1) & (HS - 1);
    if (!hocc[h]) { hocc[h] = 1; hkey[h] = key; hy[h] = y; hval[h].re = 0; hval[h].im = 0; hused[nused++] = (int)h; if (nused > HS / 2) { fprintf(stderr, "hash full\n"); exit(3); } }
    hval[h].re += v.re; hval[h].im += v.im;
}
static void hclear(void) { for (int i = 0; i < nused; i++) hocc[hused[i]] = 0; nused = 0; }
/* wedge of the basis vectors e_{c_1} ^ ... ^ e_{c_k} (c in the given order): returns 0 if a repeat, else the sign
 * of the sorting permutation, and the sorted set as a mask */
static int wedge_sign(const int *c, int len, u64 *mask) {
    u64 m = 0; int inv = 0;
    for (int i = 0; i < len; i++) { if (m >> c[i] & 1) return 0; m |= 1ULL << c[i]; }
    for (int i = 0; i < len; i++) for (int j = i + 1; j < len; j++) if (c[i] > c[j]) inv++;
    *mask = m; return (inv & 1) ? -1 : 1;
}
/* add coef * G_W (e_x (x) e_T) to the table, G = (1+i) Id + (1-i) X_z */
static void add_image(int x, u64 T, gi coef, int z) {
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
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: check_capA m k [z ...]\n"); return 2; }
    int m = atoi(argv[1]); k = atoi(argv[2]); n = 1 << m;
    if (m < 1 || m > 6 || k < 1 || k > n) { fprintf(stderr, "range\n"); return 2; }
    /* dim N, counted directly: x in [0,n), T a k-subset containing x */
    nsub = 0; subs(0, k, 0); u64 dimN = (u64)k * nsub;
    int zs[64], nz = 0;
    if (argc > 3) for (int i = 3; i < argc; i++) zs[nz++] = atoi(argv[i]); else for (int z = 1; z < n; z++) zs[nz++] = z;
    /* (k-1)-subsets of [0,n) */
    nsub = 0; subs(0, k - 1, 0); u64 *W = SUB; size_t nw = nsub; SUB = NULL; csub = 0;
    /* supports of the nonzero u, collected and sorted to detect overlaps */
    size_t cap_sk = 2 * (size_t)(n / 2) * nw + 16; sk *sup = malloc(cap_sk * sizeof(sk)); size_t nsk;
    int allpass = 1;
    for (int iz = 0; iz < nz; iz++) {
        int z = zs[iz]; if (z <= 0 || z >= n) { fprintf(stderr, "bad z\n"); return 2; }
        nsk = 0;
        u64 nonzero = 0, notinN = 0, imgout = 0, overlap = 0;
        for (int a = 0; a < n; a++) {
            int b = a ^ z; if (b < a) continue;            /* unordered pair {a, b}, a < b */
            for (size_t iw = 0; iw < nw; iw++) {
                u64 w = W[iw];
                /* the two terms: e_x (x) (e_x ^ e_w), x = a with +1, x = b with -1 */
                int tx[2] = { a, b }; int tsg[2] = { 1, -1 }; u64 tT[2]; int tc[2];
                for (int q = 0; q < 2; q++) {
                    int x = tx[q]; int c[64], L = 0; c[L++] = x; for (int j = 0; j < n; j++) if (w >> j & 1) c[L++] = j;
                    u64 M = 0; int s = wedge_sign(c, L, &M); tc[q] = s * tsg[q]; tT[q] = M;
                }
                if (!tc[0] && !tc[1]) continue;              /* u = 0 */
                nonzero++;
                for (int q = 0; q < 2; q++) if (tc[q]) {
                    if (!(tT[q] >> tx[q] & 1)) notinN++;
                    sup[nsk].T = tT[q]; sup[nsk].x = tx[q]; nsk++;
                }
                hclear();
                for (int q = 0; q < 2; q++) if (tc[q]) { gi cf = { tc[q], 0 }; add_image(tx[q], tT[q], cf, z); }
                for (int i = 0; i < nused; i++) {
                    int h = hused[i]; gi v = hval[h]; if (!v.re && !v.im) continue;
                    int y = hy[h]; u64 T = hkey[h];
                    if (!(T >> y & 1)) imgout++;
                }
            }
        }
        qsort(sup, nsk, sizeof(sk), skcmp);
        for (size_t i = 1; i < nsk; i++) if (sup[i].T == sup[i - 1].T && sup[i].x == sup[i - 1].x) overlap++;
        int pass = notinN == 0 && imgout == 0 && overlap == 0;
        if (!pass) allpass = 0;
        printf("[%s] m=%d k=%d dir %d: nonzero u %llu, u outside N %llu, image coordinates outside N %llu, support overlaps %llu"
               " -> f(dir %d) <= dim N - %llu = %llu - %llu = %llu\n", pass ? "PASS" : "FAIL", m, k, z,
               (unsigned long long)nonzero, (unsigned long long)notinN, (unsigned long long)imgout, (unsigned long long)overlap,
               z, (unsigned long long)nonzero, (unsigned long long)dimN, (unsigned long long)nonzero, (unsigned long long)(dimN - nonzero));
        fflush(stdout);
    }
    printf("SUMMARY m=%d k=%d: dim N = %llu; %s\n", m, k, (unsigned long long)dimN, allpass ? "all directions PASS" : "FAIL");
    return allpass ? 0 : 1;
}
