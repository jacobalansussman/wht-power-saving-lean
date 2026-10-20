// o3-ratio: rank mod p (p < 2^25) of a dense R x C matrix. Input file: uint32 R, C, p, then R*C uint32 row-major.
// Prints "rank r". Lazy reduction as in wave16/tb4-tight/nullsp.c (entries uint64; each update adds < 2^50).
// Multi-threaded row updates with OpenMP when available.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static uint64_t P;
static uint64_t pw(uint64_t b, uint64_t e) { uint64_t r = 1; b %= P; while (e) { if (e & 1) r = r * b % P; b = b * b % P; e >>= 1; } return r; }
int main(int argc, char **argv) {
  FILE *fh = fopen(argv[1], "rb"); uint32_t hdr[3];
  if (!fh || fread(hdr, 4, 3, fh) != 3) return 1;
  long R = hdr[0], C = hdr[1]; P = hdr[2];
  uint64_t *M = malloc(sizeof(uint64_t) * R * C); uint32_t *row = malloc(4 * C);
  for (long i = 0; i < R; i++) { if (fread(row, 4, C, fh) != (size_t)C) return 1; for (long j = 0; j < C; j++) M[i * C + j] = row[j] % P; }
  fclose(fh);
  int *cnt = calloc(R, sizeof(int)); long r = 0;
  for (long c = 0; c < C && r < R; c++) {
    long piv = -1;
    for (long i = r; i < R; i++) { M[i * C + c] %= P; if (M[i * C + c]) { piv = i; break; } }
    if (piv < 0) continue;
    if (piv != r) { for (long j = 0; j < C; j++) { uint64_t t = M[r * C + j]; M[r * C + j] = M[piv * C + j]; M[piv * C + j] = t; }
                    int t = cnt[r]; cnt[r] = cnt[piv]; cnt[piv] = t; }
    uint64_t *rr = M + r * C;
    for (long j = c; j < C; j++) rr[j] %= P;
    uint64_t inv = pw(rr[c], P - 2);
    for (long j = c; j < C; j++) rr[j] = rr[j] * inv % P;
    #pragma omp parallel for schedule(dynamic, 16)
    for (long i = r + 1; i < R; i++) {
      uint64_t *ri = M + i * C;
      uint64_t f = ri[c] % P;
      if (!f) { ri[c] = 0; continue; }
      if (cnt[i] >= 8000) { for (long j = c; j < C; j++) ri[j] %= P; cnt[i] = 0; }
      uint64_t g = P - f;
      for (long j = c; j < C; j++) ri[j] += g * rr[j];
      cnt[i]++;
    }
    r++;
    if (r % 1000 == 0) { fprintf(stderr, "pivots %ld (col %ld)\n", r, c); }
  }
  printf("rank %ld\n", r);
  return 0;
}
