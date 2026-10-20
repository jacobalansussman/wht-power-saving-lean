// tb4-tight: nullspace mod p (p < 2^25) of a dense R x C matrix. Input: uint32 R, C, p, then R*C uint32 row-major.
// Output file: uint32 nnull, C, then nnull*C uint32 (rows = basis of {x : M x = 0}). Lazy reduction: entries are uint64,
// each row update adds (p - f) * pivotrow[j] < 2^50, so up to 8192 updates fit before a reduction is forced.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static uint64_t P;
static uint64_t pw(uint64_t b, uint64_t e) { uint64_t r = 1; b %= P; while (e) { if (e & 1) r = r * b % P; b = b * b % P; e >>= 1; } return r; }
int main(int argc, char **argv) {
  FILE *fh = fopen(argv[1], "rb"); uint32_t hdr[3];
  if (fread(hdr, 4, 3, fh) != 3) return 1;
  long R = hdr[0], C = hdr[1]; P = hdr[2];
  uint32_t *tmp = malloc(sizeof(uint32_t) * R * C);
  if (fread(tmp, 4, R * C, fh) != (size_t)(R * C)) return 1;
  fclose(fh);
  uint64_t *M = malloc(sizeof(uint64_t) * R * C);
  for (long i = 0; i < R * C; i++) M[i] = tmp[i] % P;
  free(tmp);
  int *cnt = calloc(R, sizeof(int));          // updates since last reduction, per row
  long *pivcol = malloc(sizeof(long) * (R < C ? R : C)); long r = 0;
  for (long c = 0; c < C && r < R; c++) {
    long piv = -1;
    for (long i = r; i < R; i++) { M[i * C + c] %= P; if (M[i * C + c]) { piv = i; break; } }
    if (piv < 0) continue;
    if (piv != r) { for (long j = 0; j < C; j++) { uint64_t t = M[r * C + j]; M[r * C + j] = M[piv * C + j]; M[piv * C + j] = t; }
                    int t = cnt[r]; cnt[r] = cnt[piv]; cnt[piv] = t; }
    uint64_t *rr = M + r * C;
    for (long j = c; j < C; j++) rr[j] %= P;
    cnt[r] = 0;
    uint64_t inv = pw(rr[c], P - 2);
    for (long j = c; j < C; j++) rr[j] = rr[j] * inv % P;
    for (long i = r + 1; i < R; i++) {
      uint64_t *ri = M + i * C;
      uint64_t f = ri[c] % P;
      if (!f) { ri[c] = 0; continue; }
      if (cnt[i] >= 8000) { for (long j = c; j < C; j++) ri[j] %= P; cnt[i] = 0; }
      uint64_t g = P - f;
      for (long j = c; j < C; j++) ri[j] += g * rr[j];
      cnt[i]++;
    }
    pivcol[r] = c; r++;
    if (r % 500 == 0) { fprintf(stderr, "pivots %ld (col %ld)\n", r, c); }
  }
  for (long i = 0; i < R * C; i++) M[i] %= P;
  // echelon rows 0..r-1 with leading 1 at pivcol; free columns -> basis vectors by back substitution
  char *isp = calloc(C, 1); for (long i = 0; i < r; i++) isp[pivcol[i]] = 1;
  long nf = 0; for (long j = 0; j < C; j++) if (!isp[j]) nf++;
  fprintf(stderr, "rank %ld, nullity %ld\n", r, nf);
  FILE *fo = fopen(argv[2], "wb"); uint32_t oh[2] = {(uint32_t)nf, (uint32_t)C}; fwrite(oh, 4, 2, fo);
  uint64_t *x = malloc(sizeof(uint64_t) * C); uint32_t *xo = malloc(sizeof(uint32_t) * C);
  for (long fc = 0; fc < C; fc++) {
    if (isp[fc]) continue;
    memset(x, 0, sizeof(uint64_t) * C); x[fc] = 1;
    for (long i = r - 1; i >= 0; i--) {
      long pc = pivcol[i]; uint64_t s = 0; uint64_t *ri = M + i * C;
      for (long j = pc + 1; j < C; j++) if (x[j]) s = (s + ri[j] * x[j]) % P;
      x[pc] = (P - s) % P;
    }
    for (long j = 0; j < C; j++) xo[j] = (uint32_t)x[j];
    fwrite(xo, 4, C, fo);
  }
  fclose(fo);
  return 0;
}
