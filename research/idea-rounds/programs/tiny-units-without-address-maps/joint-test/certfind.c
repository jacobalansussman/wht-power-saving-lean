/* joint/cert/certfind.c -- agent "cert", 10 Oct 2026.  Finds CUT CERTIFICATES for the joint routing test.
   No max-flow code: for every set B of sink blocks (inside one of the two families of blocks, bit 0: blocks 0..7,
   bit 1: blocks 8..14) the cheapest cut is found by exhaustive dynamic programming over nested sets
   U_0 >= U_1 >= ... >= U_K = B:  f_0(U) = W * (labels the kernel carries from U \ B into B),
   f_t(U) = min over U' >= U of f_{t-1}(U') + sum_{s not in U', d in U} cap[w_t][s][d].
   (U_K = B loses nothing: a smaller last set only removes terms.)  A certificate exists for B iff f_K(B) < need(B),
   need(B) = W * (labels the kernel carries into B from outside B).
   Usage: certfind model.txt W K all LO HI out.bin stats.txt     (all 7^K ordered patterns whose first direction is LO..HI)
          certfind model.txt W K list patterns.txt out.bin stats.txt
   Record in out.bin: K bytes (the directions), uint16 B, uint16 U_0 .. U_K (bit j = block j).  B = 0: no certificate. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifndef MAXK
#define MAXK 12
#endif
#define MAXB 400
static int cap[8][15][15], dem[15][15];
static const int OFF[2] = {0, 8}, NC[2] = {8, 7};
static unsigned char Cst[2][8][256][256];
typedef struct { int c, b, need; int dinto[8]; } Bset;
static Bset BS[MAXB]; static int NBS;
static short f[MAXK + 1][MAXB][256];
static int W, K, pat[MAXK + 1];
static long long nwith, nnone, cntFirst[MAXB], cntBest[MAXB], histBest[200], histNrej[MAXB + 1];
static FILE *fout, *fstat;
static void init(void) {
    for (int c = 0; c < 2; c++) for (int w = 1; w < 8; w++) for (int U1 = 0; U1 < (1 << NC[c]); U1++) for (int U = 0; U < (1 << NC[c]); U++) {
        int v = 0;
        if ((U & ~U1) == 0) for (int s = 0; s < NC[c]; s++) if (!((U1 >> s) & 1)) for (int d = 0; d < NC[c]; d++) if ((U >> d) & 1) v += cap[w][OFF[c] + s][OFF[c] + d];
        Cst[c][w][U1][U] = (unsigned char)v;
    }
    NBS = 0;
    for (int pc = 1; pc < 8; pc++) for (int c = 0; c < 2; c++) for (int b = 1; b < (1 << NC[c]) - 1; b++) if (__builtin_popcount(b) == pc) {
        Bset *B = &BS[NBS]; B->c = c; B->b = b; B->need = 0;
        for (int s = 0; s < NC[c]; s++) { B->dinto[s] = 0; if (!((b >> s) & 1)) for (int d = 0; d < NC[c]; d++) if ((b >> d) & 1) B->dinto[s] += dem[OFF[c] + s][OFF[c] + d]; B->need += W * B->dinto[s]; }
        int full = (1 << NC[c]) - 1, rest = full ^ b;
        for (int sub = rest;; sub = (sub - 1) & rest) { int v = 0; for (int s = 0; s < NC[c]; s++) if ((sub >> s) & 1) v += W * B->dinto[s]; f[0][NBS][b | sub] = (short)v; if (!sub) break; }
        NBS++;
    }
}
static void step(int t, int w) {            /* f[t] from f[t-1], move t in direction w, all sets B */
    for (int j = 0; j < NBS; j++) {
        int c = BS[j].c, b = BS[j].b, full = (1 << NC[c]) - 1, rest = full ^ b;
        const short *g = f[t - 1][j]; short *h = f[t][j];
        for (int sub = rest;; sub = (sub - 1) & rest) {
            int U = b | sub, r2 = full ^ U, best = 1 << 20;
            for (int s2 = r2;; s2 = (s2 - 1) & r2) { int v = g[U | s2] + Cst[c][w][U | s2][U]; if (v < best) best = v; if (!s2) break; }
            h[U] = (short)best; if (!sub) break;
        }
    }
}
static int last(int j, int w) {             /* f_K(B) for set j, last move w */
    int c = BS[j].c, b = BS[j].b, full = (1 << NC[c]) - 1, rest = full ^ b, best = 1 << 20;
    const short *g = f[K - 1][j];
    for (int s2 = rest;; s2 = (s2 - 1) & rest) { int v = g[b | s2] + Cst[c][w][b | s2][b]; if (v < best) best = v; if (!s2) break; }
    return best;
}
static void leaf(void) {
    int first = -1, bestj = -1, bestdef = 0, nrej = 0, w = pat[K];
    for (int j = 0; j < NBS; j++) { int def = BS[j].need - last(j, w); if (def > 0) { nrej++; if (first < 0) first = j; if (def > bestdef) { bestdef = def; bestj = j; } } }
    unsigned char rec[MAXK]; uint16_t u[MAXK + 2];
    for (int t = 1; t <= K; t++) rec[t - 1] = (unsigned char)pat[t];
    memset(u, 0, sizeof u); histNrej[nrej]++;
    if (first < 0) { nnone++; fprintf(fstat, "NONE"); for (int t = 1; t <= K; t++) fprintf(fstat, " %d", pat[t]); fprintf(fstat, "\n"); }
    else {
        nwith++; cntFirst[first]++; cntBest[bestj]++; histBest[bestdef]++;
        int j = bestj, c = BS[j].c, full = (1 << NC[c]) - 1, U = BS[j].b;
        u[0] = (uint16_t)(BS[j].b << OFF[c]); u[K + 1] = (uint16_t)(U << OFF[c]);
        for (int t = K; t >= 1; t--) {      /* walk back through the minimising sets */
            int r2 = full ^ U, best = 1 << 20, arg = -1; const short *g = f[t - 1][j];
            for (int s2 = r2;; s2 = (s2 - 1) & r2) { int v = g[U | s2] + Cst[c][pat[t]][U | s2][U]; if (v < best) { best = v; arg = U | s2; } if (!s2) break; }
            U = arg; u[t] = (uint16_t)(U << OFF[c]);
        }
    }
    fwrite(rec, 1, K, fout); fwrite(u, 2, K + 2, fout);
}
static void rec_all(int t, int lo, int hi) {
    for (int w = (t == 1 ? lo : 1); w <= (t == 1 ? hi : 7); w++) { pat[t] = w; if (t < K) { step(t, w); rec_all(t + 1, lo, hi); } else leaf(); }
}
static void name(char *o, int j) { int c = BS[j].c; *o = 0; for (int s = 0; s < NC[c]; s++) if ((BS[j].b >> s) & 1) sprintf(o + strlen(o), "(%d,%d)", c == 0 ? s : s + 1, c); }
int main(int argc, char **argv) {
    FILE *fm = fopen(argv[1], "r"); W = atoi(argv[2]); K = atoi(argv[3]);
    for (int w = 1; w < 8; w++) for (int s = 0; s < 15; s++) for (int d = 0; d < 15; d++) if (fscanf(fm, "%d", &cap[w][s][d]) != 1) return 2;
    for (int s = 0; s < 15; s++) for (int d = 0; d < 15; d++) if (fscanf(fm, "%d", &dem[s][d]) != 1) return 2;
    for (int s = 0; s < 15; s++) for (int d = 0; d < 15; d++) if ((s < 8) != (d < 8)) { if (dem[s][d]) return 3; for (int w = 1; w < 8; w++) if (cap[w][s][d]) return 3; }
    if (K < 1 || K > MAXK) return 4;
    init();
    if (!strcmp(argv[4], "all")) { fout = fopen(argv[7], "wb"); fstat = fopen(argv[8], "w"); rec_all(1, atoi(argv[5]), atoi(argv[6])); }
    else {
        FILE *fl = fopen(argv[5], "r"); fout = fopen(argv[6], "wb"); fstat = fopen(argv[7], "w"); int d, t = 0;
        while (fscanf(fl, "%d", &d) == 1) { pat[++t] = d; if (t < K) step(t, d); else { leaf(); t = 0; } }
    }
    fprintf(fstat, "W %d K %d SETS %d\nWITH %lld\nNOCERT %lld\n", W, K, NBS, nwith, nnone);
    for (int j = 0; j < NBS; j++) if (cntFirst[j] || cntBest[j]) { char nm[200]; name(nm, j); fprintf(fstat, "SET %d %s need %d first %lld best %lld\n", j, nm, BS[j].need, cntFirst[j], cntBest[j]); }
    for (int d = 0; d < 200; d++) if (histBest[d]) fprintf(fstat, "MARGIN %d %lld\n", d, histBest[d]);
    for (int r = 0; r <= NBS; r++) if (histNrej[r]) fprintf(fstat, "NREJ %d %lld\n", r, histNrej[r]);
    fclose(fout); fclose(fstat);
    printf("W=%d K=%d sets=%d patterns with a certificate: %lld ; without: %lld\n", W, K, NBS, nwith, nnone);
    return 0;
}
