// tb4-beam: dense reduced row echelon form mod p (p < 2^21) in place, uint64 storage, delayed reduction.
// rref(A, R, C, p, piv): returns rank r; rows 0..r-1 of A become the reduced pivot rows (entries < p),
// rows r..R-1 become zero mod p (not stored reduced); piv[0..r-1] = pivot columns.
#include <stdint.h>
#include <string.h>
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t p){uint64_t r=1;a%=p;while(e){if(e&1)r=r*a%p;a=a*a%p;e>>=1;}return r;}
long rref(uint64_t *A, long R, long C, uint64_t p, long *piv, long full) {
  long r = 0; long since = 0;
  for (long c = 0; c < C && r < R; c++) {
    long pr = -1;
    for (long i = r; i < R; i++) { uint64_t v = A[i*C+c] % p; A[i*C+c] = v; if (v) { pr = i; break; } }
    if (pr < 0) continue;
    if (pr != r) { for (long j = 0; j < C; j++) { uint64_t t = A[r*C+j]; A[r*C+j] = A[pr*C+j]; A[pr*C+j] = t; } }
    uint64_t *rr = A + r*C;
    for (long j = c; j < C; j++) rr[j] %= p;
    uint64_t inv = powmod(rr[c], p-2, p);
    for (long j = c; j < C; j++) rr[j] = rr[j] * inv % p;
    long i0 = full ? 0 : r+1;
    for (long i = i0; i < R; i++) {
      if (i == r) continue;
      uint64_t *ri = A + i*C;
      uint64_t f = ri[c] % p;
      if (!f) { ri[c] = 0; continue; }
      uint64_t g = p - f;
      for (long j = c; j < C; j++) ri[j] += g * rr[j];
    }
    piv[r] = c; r++; since++;
    if (since >= 1000) { // keep entries below 2^63: each update adds < 2^42
      for (long i = 0; i < R; i++) for (long j = 0; j < C; j++) A[i*C+j] %= p;
      since = 0;
    }
  }
  for (long i = 0; i < r; i++) for (long j = 0; j < C; j++) A[i*C+j] %= p;
  return r;
}
