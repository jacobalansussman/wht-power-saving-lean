/* pf-logging, 13 Oct 2026.  check_big.c = check_m0.c with 64-bit subset masks (n = 2^m <= 64) and distinctness by
 * sorting instead of a 2^n table; the mathematics and the determinant routine are unchanged.
 * INDEPENDENT CHECKER of a nonsingular-minor certificate for the block M_alpha
 * (shares no code with gen_m0.c; different determinant routine, no subset tables, no sorting signs).
 * Definitions (PROOF.md, Lemma C):  n = 2^m, Q[u][v] = (-i)^{popcount(u xor v)} with i -> I (I^2 = -1 mod p);
 *   for a row set U = {u_1 < ... < u_k} (0 not in U) and a column set T = {t_1 < ... < t_k} (0 in T):
 *   M_alpha[U][T] = sum over s in F_2^m of (-1)^{alpha.s} * Q[0][s] * det( Q[u_a][t_b xor s] )_{a,b=1..k}.
 * The certificate file lists p, I, and `size` pairs (row mask, col mask).  The checker verifies: p is prime,
 * I^2 = -1 mod p, the row masks are distinct k-subsets of {1..n-1}, the col masks distinct k-subsets containing 0,
 * then builds the size x size matrix from the formula above and computes its determinant mod p.
 * Nonzero determinant mod p  =>  the same minor of M_alpha over Z[i] (scaled entries) is nonzero  =>
 * rank M_alpha over Q(i) >= size.
 * Build: cc -O2 -o check_big check_big.c        Run: ./check_big CERTFILE
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64;
static u64 p;
static u64 mulm(u64 a, u64 b) { return a * b % p; }   /* exact: p < 2^31 is checked, so a*b < 2^62 */
static u64 powm(u64 b, u64 e) { u64 r = 1; b %= p; while (e) { if (e & 1) r = mulm(r, b); b = mulm(b, b); e >>= 1; } return r; }
static int bits(uint64_t x) { int c = 0; while (x) { c += (int)(x & 1); x >>= 1; } return c; }
static int cmpu(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return x < y ? -1 : x > y; }
/* determinant mod p of an N x N matrix (row-major, destroyed), Gaussian elimination with row swaps */
static u64 detm(u64 *A, int N) {
    u64 det = 1;
    for (int c = 0; c < N; c++) {
        int piv = -1;
        for (int r = c; r < N; r++) if (A[(size_t)r * N + c] % p) { piv = r; break; }
        if (piv < 0) return 0;
        if (piv != c) {
            for (int j = 0; j < N; j++) { u64 t = A[(size_t)c * N + j]; A[(size_t)c * N + j] = A[(size_t)piv * N + j]; A[(size_t)piv * N + j] = t; }
            det = (p - det) % p;
        }
        u64 a = A[(size_t)c * N + c]; det = mulm(det, a);
        u64 ia = powm(a, p - 2);
        for (int r = c + 1; r < N; r++) {
            u64 b = A[(size_t)r * N + c]; if (!b) continue;
            u64 f = mulm(b, ia);
            u64 *rr = A + (size_t)r * N, *rc = A + (size_t)c * N;
            for (int j = c; j < N; j++) rr[j] = (rr[j] + p - mulm(f, rc[j])) % p;
        }
    }
    return det;
}
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: check_m0 CERT\n"); return 2; }
    FILE *fi = fopen(argv[1], "r"); if (!fi) { perror("cert"); return 2; }
    int m, k, size; unsigned long long alpha; unsigned long long pp, II;
    if (fscanf(fi, "m %d k %d alpha %llu p %llu I %llu size %d", &m, &k, &alpha, &pp, &II, &size) != 6) { fprintf(stderr, "bad header\n"); return 2; }
    p = pp; u64 I = II; int n = 1 << m; if (m < 1 || m > 6) { fprintf(stderr, "m out of range\n"); return 2; }
    /* p prime (trial division), I^2 = -1 */
    int prime = p > 2 && p < (1ULL << 31);
    for (u64 d = 2; d * d <= p && prime; d++) if (p % d == 0) prime = 0;
    int okI = mulm(I, I) == p - 1;
    u64 *R = malloc(sizeof(u64) * size), *C = malloc(sizeof(u64) * size);
    for (int i = 0; i < size; i++) if (fscanf(fi, "%llu %llu", (unsigned long long *)&R[i], (unsigned long long *)&C[i]) != 2) { fprintf(stderr, "short cert\n"); return 2; }
    fclose(fi);
    int okS = 1; u64 full = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
    for (int i = 0; i < size; i++)
        if ((R[i] & ~full) || (C[i] & ~full) || bits(R[i]) != k || bits(C[i]) != k || (R[i] & 1) || !(C[i] & 1)) okS = 0;
    u64 *Rs = malloc(sizeof(u64) * size), *Cs = malloc(sizeof(u64) * size);
    for (int i = 0; i < size; i++) { Rs[i] = R[i]; Cs[i] = C[i]; }
    qsort(Rs, size, sizeof(u64), cmpu); qsort(Cs, size, sizeof(u64), cmpu);
    for (int i = 1; i < size; i++) if (Rs[i] == Rs[i - 1] || Cs[i] == Cs[i - 1]) okS = 0;
    printf("cert: m=%d k=%d alpha=%llu p=%llu (prime: %s) I^2=-1: %s, %d distinct well-formed rows/cols: %s\n",
           m, k, alpha, pp, prime ? "yes" : "NO", okI ? "yes" : "NO", size, okS ? "yes" : "NO");
    if (!(prime && okI && okS)) { printf("RESULT: FAIL (format)\n"); return 1; }
    u64 mi[4] = {1, p - I, p - 1, I};                 /* (-i)^0, (-i)^1, (-i)^2, (-i)^3 */
    u64 *A = malloc(sizeof(u64) * (size_t)size * size); if (!A) { fprintf(stderr, "alloc\n"); return 1; }
    u64 *B = malloc(sizeof(u64) * k * k);
    for (int i = 0; i < size; i++) {
        int u[32], L = 0; for (int x = 0; x < n; x++) if (R[i] >> x & 1) u[L++] = x;
        for (int j = 0; j < size; j++) {
            int t[32], L2 = 0; for (int x = 0; x < n; x++) if (C[j] >> x & 1) t[L2++] = x;
            u64 acc = 0;
            for (int s = 0; s < n; s++) {
                for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) B[a * k + b] = mi[bits((u64)(u[a] ^ (t[b] ^ s))) & 3];
                u64 d = detm(B, k);
                u64 term = mulm(mi[bits((u64)s) & 3], d);
                if (bits(alpha & (u64)s) & 1) term = (p - term) % p;
                acc = (acc + term) % p;
            }
            A[(size_t)i * size + j] = acc;
        }
        if (i % 500 == 0) { fprintf(stderr, "  built row %d / %d\n", i, size); }
    }
    u64 d = detm(A, size);
    printf("det mod p of the %d x %d minor of M_%llu = %llu\n", size, size, alpha, (unsigned long long)d);
    printf("RESULT: %s  -> rank M_%llu over Q(i) >= %d%s\n", d ? "PASS" : "FAIL", alpha, d ? size : 0, d ? "" : " (not certified)");
    return d ? 0 : 1;
}
