/* pf-logging, 13 Oct 2026.  GENERATOR (untrusted) of nonsingular-minor certificates for M_alpha (Lemma C of PROOF.md)
 * at larger block sizes (n = 2^m <= 64, subsets as 64-bit masks), with a random sample of the rows.
 *   M_alpha[U][T] = sum_s (-1)^{alpha.s} Q[0][s] det( Q[u_a][t_b xor s] ),  Q[u][v] = (-i)^{popcount(u xor v)}
 *   rows: k-subsets U of {1..n-1}; cols: k-subsets T of {0..n-1} with 0 in T.
 * Build: cc -O2 -DP=2147483629ULL -o gen_big gen_big.c     Run: ./gen_big m k alpha NSAMPLE SEED OUTFILE
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#ifndef P
#define P 2147483629ULL
#endif
typedef uint64_t u64;
static u64 pw(u64 b, u64 e) { u64 r = 1; b %= P; while (e) { if (e & 1) r = r * b % P; b = b * b % P; e >>= 1; } return r; }
static int pc(u64 x) { return __builtin_popcountll(x); }
static u64 rng = 88172645463325252ULL;
static u64 xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static int K, N;
static u64 pmi[4];
static void subsets(int lo, int need, u64 cur, u64 **out, int *cnt, int *cap) {
    if (need == 0) { if (*cnt == *cap) { *cap *= 2; *out = realloc(*out, sizeof(u64) * *cap); } (*out)[(*cnt)++] = cur; return; }
    for (int x = lo; x <= N - need; x++) subsets(x + 1, need - 1, cur | (1ULL << x), out, cnt, cap);
}
static u64 det3(const u64 *A, int k, int r0, int c0, int c1, int c2) {   /* rows r0..r0+2 of A, columns c0,c1,c2 */
    const u64 *a = A + r0 * k, *b = A + (r0 + 1) * k, *c = A + (r0 + 2) * k;
    u64 t1 = a[c0] * ((b[c1] * c[c2] + (P - b[c2]) * c[c1]) % P) % P;
    u64 t2 = a[c1] * ((b[c0] * c[c2] + (P - b[c2]) * c[c0]) % P) % P;
    u64 t3 = a[c2] * ((b[c0] * c[c1] + (P - b[c1]) * c[c0]) % P) % P;
    return (t1 + (P - t2) + t3) % P;
}
static u64 detk(u64 *A, int k) {   /* destroys A (general case) */
    if (k == 3) return det3(A, 3, 0, 0, 1, 2);
    if (k == 4) {                  /* Laplace along row 0 */
        u64 r = 0; const int cols[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
        for (int j = 0; j < 4; j++) { u64 t = A[j] * det3(A, 4, 1, cols[j][0], cols[j][1], cols[j][2]) % P; r += (j & 1) ? P - t : t; }
        return r % P;
    }
    u64 det = 1;
    for (int c = 0; c < k; c++) {
        int pr = -1; for (int r = c; r < k; r++) if (A[r * k + c]) { pr = r; break; }
        if (pr < 0) return 0;
        if (pr != c) { for (int j = 0; j < k; j++) { u64 t = A[c * k + j]; A[c * k + j] = A[pr * k + j]; A[pr * k + j] = t; } det = (P - det) % P; }
        det = det * A[c * k + c] % P; u64 iv = pw(A[c * k + c], P - 2);
        for (int r = c + 1; r < k; r++) if (A[r * k + c]) {
            u64 f = A[r * k + c] * iv % P;
            for (int j = c; j < k; j++) A[r * k + j] = (A[r * k + j] + (P - f) * A[c * k + j]) % P;
        }
    }
    return det;
}
int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: gen_big m k alpha NSAMPLE SEED OUTFILE\n"); return 2; }
    int m = atoi(argv[1]); K = atoi(argv[2]); u64 alpha = strtoull(argv[3], 0, 10); int ns = atoi(argv[4]);
    rng ^= strtoull(argv[5], 0, 10) * 0x9E3779B97F4A7C15ULL; N = 1 << m;
    u64 I = 0; for (u64 g = 2; g < 1000; g++) if (pw(g, (P - 1) / 2) == P - 1) { I = pw(g, (P - 1) / 4); break; }
    pmi[0] = 1; pmi[1] = (P - I) % P; pmi[2] = P - 1; pmi[3] = I;
    int cr = 0, cc = 0, capr = 1024, capc = 1024; u64 *R = malloc(sizeof(u64) * capr), *C = malloc(sizeof(u64) * capc);
    subsets(1, K, 0, &R, &cr, &capr);            /* rows: k-subsets of {1..n-1} */
    subsets(1, K - 1, 1ULL, &C, &cc, &capc);     /* cols: 0 plus a (k-1)-subset of {1..n-1} */
    if (ns > 0 && ns < cr) {                     /* random sample of the rows (Fisher-Yates prefix) */
        for (int i = 0; i < ns; i++) { int j = i + (int)(xr() % (u64)(cr - i)); u64 t = R[i]; R[i] = R[j]; R[j] = t; }
        cr = ns;
    }
    fprintf(stderr, "m=%d k=%d alpha=%llu: %d rows (sampled), %d cols\n", m, K, (unsigned long long)alpha, cr, cc);
    u64 *M = malloc(sizeof(u64) * (size_t)cr * cc); u64 A[64];
    for (int r = 0; r < cr; r++) {
        int u[64], L = 0; for (int x = 0; x < N; x++) if (R[r] >> x & 1) u[L++] = x;
        for (int c = 0; c < cc; c++) {
            int t[64], L2 = 0; for (int x = 0; x < N; x++) if (C[c] >> x & 1) t[L2++] = x;
            u64 acc = 0;
            for (int s = 0; s < N; s++) {
                for (int a = 0; a < K; a++) for (int b = 0; b < K; b++) A[a * K + b] = pmi[pc((u64)(u[a] ^ t[b] ^ s)) & 3];
                u64 term = pmi[pc((u64)s) & 3] * detk(A, K) % P;
                acc += (pc(alpha & (u64)s) & 1) ? (P - term) % P : term;
            }
            M[(size_t)r * cc + c] = acc % P;
        }
    }
    fprintf(stderr, "matrix built\n");
    char *used = calloc(cr, 1); int *pr = malloc(sizeof(int) * cc), *pcl = malloc(sizeof(int) * cc), rank = 0;
    for (int c = 0; c < cc; c++) {
        int p0 = -1; for (int r = 0; r < cr; r++) if (!used[r] && M[(size_t)r * cc + c]) { p0 = r; break; }
        if (p0 < 0) continue;
        used[p0] = 1; pr[rank] = p0; pcl[rank] = c; rank++;
        u64 *rp = M + (size_t)p0 * cc; u64 iv = pw(rp[c], P - 2);
        for (int r = 0; r < cr; r++) {
            if (used[r]) continue; u64 *rr = M + (size_t)r * cc; if (!rr[c]) continue;
            u64 f = (P - rr[c] * iv % P) % P;
            for (int j = c; j < cc; j++) rr[j] = (rr[j] + f * rp[j]) % P;
        }
    }
    printf("m=%d k=%d alpha=%llu: rank mod %llu of the sampled M_alpha (%d x %d) = %d\n", m, K, (unsigned long long)alpha,
           (unsigned long long)P, cr, cc, rank);
    FILE *fo = fopen(argv[6], "w");
    fprintf(fo, "m %d k %d alpha %llu p %llu I %llu size %d\n", m, K, (unsigned long long)alpha, (unsigned long long)P, (unsigned long long)I, rank);
    for (int i = 0; i < rank; i++) fprintf(fo, "%llu %llu\n", (unsigned long long)R[pr[i]], (unsigned long long)C[pcl[i]]);
    fclose(fo);
    return 0;
}
