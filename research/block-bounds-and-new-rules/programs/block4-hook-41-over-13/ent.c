// b56-compute (wave 16).  Entries of the Walsh-coefficient matrix C[u, t] = prod_i det( H[U_i, T_i] ),
// H[s, x] = (-1)^{popcount(s & x)}, for Walsh tuples u = (U_1..U_r) (rows) and standard tuples t = (T_1..T_r) (cols).
// Exact integer determinants (Bareiss, int64), returned as doubles (exact: |entries| < 2^40).
// Build: cc -O3 -shared -fPIC -o libent.so ent.c -lpthread
#include <stdint.h>
#include <stdlib.h>
#include <pthread.h>

static int64_t det_pm(int k, const uint8_t *U, const uint8_t *T) {
    int64_t M[12][12];
    if (k == 0) return 1;
    for (int a = 0; a < k; a++)
        for (int b = 0; b < k; b++)
            M[a][b] = (__builtin_popcount((unsigned)(U[a] & T[b])) & 1) ? -1 : 1;
    if (k == 1) return M[0][0];
    int64_t prev = 1; int sign = 1;
    for (int c = 0; c < k - 1; c++) {
        if (M[c][c] == 0) {
            int r = c + 1;
            while (r < k && M[r][c] == 0) r++;
            if (r == k) return 0;
            for (int j = 0; j < k; j++) { int64_t t = M[c][j]; M[c][j] = M[r][j]; M[r][j] = t; }
            sign = -sign;
        }
        for (int r = c + 1; r < k; r++)
            for (int j = c + 1; j < k; j++)
                M[r][j] = (M[r][j] * M[c][c] - M[r][c] * M[c][j]) / prev;
        prev = M[c][c];
    }
    return sign * M[k - 1][k - 1];
}

typedef struct {
    int r, f; const int *sz; const uint8_t *rows; const uint8_t *cols; int ncols; double *out; int r0, r1;
} job_t;

static void *work(void *arg) {
    job_t *J = (job_t *)arg;
    for (int i = J->r0; i < J->r1; i++) {
        const uint8_t *U = J->rows + (size_t)i * J->f;
        double *o = J->out + (size_t)(i - J->r0 + J->r0) * J->ncols;
        for (int j = 0; j < J->ncols; j++) {
            const uint8_t *T = J->cols + (size_t)j * J->f;
            int64_t v = 1; int off = 0;
            for (int q = 0; q < J->r && v != 0; q++) {
                v *= det_pm(J->sz[q], U + off, T + off);
                off += J->sz[q];
            }
            o[j] = (double)v;
        }
    }
    return 0;
}

// out: nrows x ncols row-major
void fill(int r, const int *sz, int nrows, const uint8_t *rows, int ncols, const uint8_t *cols, double *out, int nthreads) {
    int f = 0; for (int q = 0; q < r; q++) f += sz[q];
    if (nthreads < 1) nthreads = 1;
    pthread_t th[64]; job_t jobs[64];
    if (nthreads > 64) nthreads = 64;
    int per = (nrows + nthreads - 1) / nthreads;
    int nt = 0;
    for (int t = 0; t < nthreads; t++) {
        int a = t * per, b = a + per; if (b > nrows) b = nrows; if (a >= b) break;
        jobs[t] = (job_t){r, f, sz, rows, cols, ncols, out, a, b};
        pthread_create(&th[t], 0, work, &jobs[t]); nt++;
    }
    for (int t = 0; t < nt; t++) pthread_join(th[t], 0);
}
