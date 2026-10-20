/* joint/cert/certcheck.c -- agent "cert", 10 Oct 2026.  Checker of cut certificates: plain sums, no flow, no search.
   Usage: certcheck model.txt W K all|list file.bin [file.bin ...]
   A record is K directions, a set B and sets U_0..U_K of blocks (bit masks over 15 blocks).  It is accepted when
   B is not empty, B <= U_K <= ... <= U_0, and
     W * (labels the kernel carries from U_0 \ B into B)  +  sum over moves t of (labels dir(w_t) carries from
     outside U_(t-1) into U_t)   <   W * (labels the kernel carries into B from outside B).
   With "all", the records must be the 7^K ordered patterns in counting order, each exactly once. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
int main(int argc, char **argv) {
    static int cap[8][15][15], dem[15][15];
    FILE *fm = fopen(argv[1], "r"); int W = atoi(argv[2]), K = atoi(argv[3]), all = !strcmp(argv[4], "all");
    for (int w = 1; w < 8; w++) for (int s = 0; s < 15; s++) for (int d = 0; d < 15; d++) if (fscanf(fm, "%d", &cap[w][s][d]) != 1) return 2;
    for (int s = 0; s < 15; s++) for (int d = 0; d < 15; d++) if (fscanf(fm, "%d", &dem[s][d]) != 1) return 2;
    long long n = 0, good = 0, bad = 0, none = 0, order_bad = 0; int minmargin = 1 << 20;
    for (int a = 5; a < argc; a++) {
        FILE *fi = fopen(argv[a], "rb"); unsigned char p[32]; uint16_t u[32];
        while (fread(p, 1, K, fi) == (size_t)K && fread(u, 2, K + 2, fi) == (size_t)(K + 2)) {
            if (all) { long long i = n; for (int t = K - 1; t >= 0; t--) { if (p[t] != 1 + i % 7) order_bad++; i /= 7; } }
            n++;
            int B = u[0], ok = 1, need = 0, src = 0, cut = 0;
            if (B == 0) { none++; continue; }
            for (int t = 0; t < K; t++) if (p[t] < 1 || p[t] > 7) ok = 0;
            if (u[K + 1] & ~0x7fff || u[1] & ~0x7fff || (B & ~u[K + 1])) ok = 0;
            for (int t = 1; t <= K; t++) if (u[t + 1] & ~u[t]) ok = 0;
            for (int s = 0; s < 15; s++) for (int d = 0; d < 15; d++) if (((B >> d) & 1) && !((B >> s) & 1)) { need += W * dem[s][d]; if ((u[1] >> s) & 1) src += W * dem[s][d]; }
            for (int t = 1; t <= K && ok; t++) for (int s = 0; s < 15; s++) if (!((u[t] >> s) & 1)) for (int d = 0; d < 15; d++) if ((u[t + 1] >> d) & 1) cut += cap[p[t - 1]][s][d];
            if (ok && src + cut < need) { good++; if (need - src - cut < minmargin) minmargin = need - src - cut; } else bad++;
        }
        fclose(fi);
    }
    long long expect = 1; for (int t = 0; t < K; t++) expect *= 7;
    printf("W=%d K=%d records %lld%s ; certificate verified %lld ; certificate WRONG %lld ; no certificate %lld ; smallest margin (need - cut) %d\n",
           W, K, n, all ? (n == expect && order_bad == 0 ? " = 7^K, every ordered pattern once, in order" : " ORDER OR COUNT WRONG") : "", good, bad, none, good ? minmargin : -1);
    return bad || (all && (n != expect || order_bad)) ? 1 : 0;
}
