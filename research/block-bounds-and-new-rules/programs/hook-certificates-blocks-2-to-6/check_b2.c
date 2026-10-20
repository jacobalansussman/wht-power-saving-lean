/* pf-logging, 14 Oct 2026.  INDEPENDENT CHECKER for the rank-2 block certificate (PROOF.md section 2c); shares no code with
 * gen_b2.c / check_yint.c (own wedge signs by insertion sort, own breadth-first component search, own expansion of g_W by
 * recursion, own determinant).
 * Certificate: header "m k z1 z2 p I size", then `size` lines "x Tmask j Rmask":
 *   (x, T) a basis element of N (x in T); c = the vector of the component of (x, T) in the graph whose edges are the
 *   two-term supports of Lemma A for z1 and for z2, with c[(x,T)] = 1 and c[p2] = -s_a s_b c[p1] along every edge;
 *   j in {0,1} picks D = 2 dir z1 or 2 dir z2; the row is the alpha = 0 projection of D_W c, read at the coordinate
 *   (0, R): sum over terms (y, T') of D_W c of sign(e_{t xor y} wedge, t in T' increasing) * coefficient, R = T' xor y.
 * Checks, all exact over Z[i] except the last: every edge relation is consistent around the component (so c is in
 * Y(z1) cap Y(z2), Lemma G), c in N, (2 dir z1)_W c in N and (2 dir z2)_W c in N; then the size x size matrix of the
 * rows at the listed columns R (all distinct k-subsets containing 0) has nonzero determinant mod p (p prime, I^2 = -1).
 * Conclusion (PROOF.md 2c): dim(D1 Z12 + D2 Z12) >= 2^m * size, hence f(dir z1 dir z2) <= dim N - 2^m * size.
 * Build: cc -O2 -o check_b2 check_b2.c      Run: ./check_b2 CERT
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef int64_t i64; typedef uint64_t u64;
static int n, k; static u64 p;
static int sortsign(int *c, int len) {           /* insertion sort in place; returns 0 on a repeat, else +-1 */
    int s = 1;
    for (int i = 1; i < len; i++) { int v = c[i], j = i - 1; while (j >= 0 && c[j] > v) { c[j + 1] = c[j]; j--; s = -s; } c[j + 1] = v; }
    for (int i = 1; i < len; i++) if (c[i] == c[i - 1]) return 0;
    return s;
}
static u64 tomask(const int *c, int len) { u64 M = 0; for (int i = 0; i < len; i++) M |= 1ULL << c[i]; return M; }
/* a sparse Gaussian-integer vector: parallel arrays, linear search merge (vectors are small) */
typedef struct { int *y; u64 *T; i64 *re, *im; int len, cap; } vec;
static void vadd(vec *v, int y, u64 T, i64 re, i64 im) {
    for (int i = 0; i < v->len; i++) if (v->y[i] == y && v->T[i] == T) { v->re[i] += re; v->im[i] += im; return; }
    if (v->len == v->cap) { v->cap = v->cap ? 2 * v->cap : 256; v->y = realloc(v->y, v->cap * sizeof(int)); v->T = realloc(v->T, v->cap * sizeof(u64));
        v->re = realloc(v->re, v->cap * sizeof(i64)); v->im = realloc(v->im, v->cap * sizeof(i64)); }
    v->y[v->len] = y; v->T[v->len] = T; v->re[v->len] = re; v->im[v->len] = im; v->len++;
}
/* recursion: coefficient (re, im) times wedge of the chosen images so far; G e_t = (1+i) e_t + (1-i) e_{t xor z} */
static void expand(vec *out, int y, const int *t, int L, int pos, int *chosen, i64 re, i64 im, int z) {
    if (pos == L) { int c[64]; memcpy(c, chosen, sizeof(int) * L); int s = sortsign(c, L); if (!s) return; vadd(out, y, tomask(c, L), s * re, s * im); return; }
    chosen[pos] = t[pos];       expand(out, y, t, L, pos + 1, chosen, re - im, re + im, z);   /* times (1+i) */
    chosen[pos] = t[pos] ^ z;   expand(out, y, t, L, pos + 1, chosen, re + im, im - re, z);   /* times (1-i) */
}
static void apply_G(vec *out, const vec *in, int z) {   /* out += (2 dir z)_W in */
    for (int i = 0; i < in->len; i++) {
        if (!in->re[i] && !in->im[i]) continue;
        int t[64], L = 0; for (int b = 0; b < n; b++) if (in->T[i] >> b & 1) t[L++] = b;
        int ch[64]; i64 a = in->re[i], b2 = in->im[i];
        expand(out, in->y[i], t, L, 0, ch, a - b2, a + b2, z);              /* e_x -> (1+i) e_x */
        expand(out, in->y[i] ^ z, t, L, 0, ch, a + b2, b2 - a, z);          /* e_x -> (1-i) e_{x xor z} */
    }
}
static int inN(const vec *v) { for (int i = 0; i < v->len; i++) if ((v->re[i] || v->im[i]) && !(v->T[i] >> v->y[i] & 1)) return 0; return 1; }
static u64 mulm(u64 a, u64 b) { return a * b % p; }
static u64 powm(u64 b, u64 e) { u64 r = 1; b %= p; while (e) { if (e & 1) r = mulm(r, b); b = mulm(b, b); e >>= 1; } return r; }
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: check_b2 CERT\n"); return 2; }
    FILE *fi = fopen(argv[1], "r"); if (!fi) { perror("cert"); return 2; }
    int m, z1, z2, size; unsigned long long pp, II;
    if (fscanf(fi, "m %d k %d z1 %d z2 %d p %llu I %llu size %d", &m, &k, &z1, &z2, &pp, &II, &size) != 7) { fprintf(stderr, "bad header\n"); return 2; }
    n = 1 << m; p = pp; u64 I = II; if (m < 1 || m > 5 || z1 <= 0 || z2 <= 0 || z1 >= n || z2 >= n || z1 == z2) { fprintf(stderr, "range\n"); return 2; }
    int prime = p > 2 && p < (1ULL << 31); for (u64 d = 2; d * d <= p && prime; d++) if (p % d == 0) prime = 0;
    int okI = mulm(I, I) == p - 1;
    int *X = malloc(sizeof(int) * size), *J = malloc(sizeof(int) * size); u64 *T0 = malloc(sizeof(u64) * size), *Rc = malloc(sizeof(u64) * size);
    for (int i = 0; i < size; i++) { unsigned long long a, b; if (fscanf(fi, "%d %llu %d %llu", &X[i], &a, &J[i], &b) != 4) { fprintf(stderr, "short\n"); return 2; } T0[i] = a; Rc[i] = b; }
    fclose(fi);
    int okfmt = 1;
    for (int i = 0; i < size; i++) {
        if (__builtin_popcountll(T0[i]) != k || !(T0[i] >> X[i] & 1) || (J[i] != 0 && J[i] != 1) || __builtin_popcountll(Rc[i]) != k || !(Rc[i] & 1)) okfmt = 0;
        for (int j = 0; j < i; j++) if (Rc[j] == Rc[i] || (X[j] == X[i] && T0[j] == T0[i] && J[j] == J[i])) okfmt = 0;
    }
    printf("cert: m=%d k=%d z=%d,%d p=%llu (prime %s) I^2=-1 %s size %d well-formed %s\n", m, k, z1, z2, pp, prime ? "yes" : "NO", okI ? "yes" : "NO", size, okfmt ? "yes" : "NO");
    if (!(prime && okI && okfmt)) { printf("RESULT: FAIL (format)\n"); return 1; }
    u64 *A = malloc(sizeof(u64) * (size_t)size * size); long exact_fail = 0, incons = 0; int maxcomp = 0;
    for (int r = 0; r < size; r++) {
        /* breadth-first search of the component of (X[r], T0[r]) with values */
        vec c = { 0 }; int qy[4096]; u64 qT[4096]; int qv[4096]; int qn = 0, qh = 0;
        qy[qn] = X[r]; qT[qn] = T0[r]; qv[qn] = 1; qn++;
        while (qh < qn) {
            int x = qy[qh]; u64 T = qT[qh]; int val = qv[qh]; qh++;
            int zz[2] = { z1, z2 };
            for (int e = 0; e < 2; e++) {
                int b = x ^ zz[e]; if (T >> b & 1) continue;                 /* singleton for this direction: no edge */
                int c1[64], c2[64], L = 0; c1[L] = x; c2[L] = b; L++;
                for (int t = 0; t < n; t++) if ((T >> t & 1) && t != x) { c1[L] = t; c2[L] = t; L++; }
                int sx = sortsign(c1, L), sb = sortsign(c2, L); u64 T2 = tomask(c2, L);
                int nv = -sx * sb * val;                                          /* u = s_x e_(x,T) - s_b e_(b,T2) */
                int found = -1; for (int i = 0; i < qn; i++) if (qy[i] == b && qT[i] == T2) { found = i; break; }
                if (found >= 0) { if (qv[found] != nv) incons++; }
                else { if (qn >= 4096) { fprintf(stderr, "component too large\n"); return 3; } qy[qn] = b; qT[qn] = T2; qv[qn] = nv; qn++; }
            }
        }
        if (qn > maxcomp) maxcomp = qn;
        for (int i = 0; i < qn; i++) vadd(&c, qy[i], qT[i], qv[i], 0);
        if (!inN(&c)) exact_fail++;
        vec d1 = { 0 }, d2 = { 0 }; apply_G(&d1, &c, z1); apply_G(&d2, &c, z2);
        if (!inN(&d1) || !inN(&d2)) exact_fail++;
        vec *d = J[r] ? &d2 : &d1;
        for (int j = 0; j < size; j++) A[(size_t)r * size + j] = 0;
        for (int i = 0; i < d->len; i++) {
            if (!d->re[i] && !d->im[i]) continue;
            int cc[64], L = 0; for (int t = 0; t < n; t++) if (d->T[i] >> t & 1) cc[L++] = t ^ d->y[i];
            int s = sortsign(cc, L); u64 R = tomask(cc, L);
            for (int j = 0; j < size; j++) if (Rc[j] == R) {
                i64 re = d->re[i] % (i64)p, im = d->im[i] % (i64)p; if (re < 0) re += p; if (im < 0) im += p;
                u64 v = ((u64)re + mulm((u64)im, I)) % p; if (s < 0) v = (p - v) % p;
                A[(size_t)r * size + j] = (A[(size_t)r * size + j] + v) % p;
            }
        }
        free(c.y); free(c.T); free(c.re); free(c.im); free(d1.y); free(d1.T); free(d1.re); free(d1.im); free(d2.y); free(d2.T); free(d2.re); free(d2.im);
        if (r % 500 == 0) fprintf(stderr, "  row %d / %d\n", r, size);
    }
    /* determinant mod p (Gaussian elimination with row swaps) */
    u64 det = 1;
    for (int c = 0; c < size && det; c++) {
        int piv = -1; for (int r = c; r < size; r++) if (A[(size_t)r * size + c]) { piv = r; break; }
        if (piv < 0) { det = 0; break; }
        if (piv != c) { for (int j = 0; j < size; j++) { u64 t = A[(size_t)c * size + j]; A[(size_t)c * size + j] = A[(size_t)piv * size + j]; A[(size_t)piv * size + j] = t; } det = (p - det) % p; }
        det = mulm(det, A[(size_t)c * size + c]); u64 iv = powm(A[(size_t)c * size + c], p - 2);
        for (int r = c + 1; r < size; r++) { u64 f = A[(size_t)r * size + c]; if (!f) continue; f = mulm(f, iv);
            for (int j = c; j < size; j++) A[(size_t)r * size + j] = (A[(size_t)r * size + j] + p - mulm(f, A[(size_t)c * size + j])) % p; }
    }
    long dimN = 0; { /* dim N = n * C(n-1, k-1) */ double b = 1; for (int i = 1; i <= k - 1; i++) b = b * (n - k + i) / i; dimN = (long)(n * b + 0.5); }
    printf("components: max size %d, sign inconsistencies %ld, exact membership failures %ld; det mod p of the %d x %d minor = %llu\n",
           maxcomp, incons, exact_fail, size, size, (unsigned long long)det);
    int ok = det && !incons && !exact_fail;
    printf("RESULT: %s -> dim(D1 Z12 + D2 Z12) >= %d x %d = %ld, f(dir %d dir %d) <= %ld - %ld = %ld\n", ok ? "PASS" : "FAIL",
           n, size, (long)n * size, z1, z2, dimN, (long)n * size, dimN - (long)n * size);
    return ok ? 0 : 1;
}
