/* ctr-window / onesided.c (own code). One-sided clique floor for P112: for each d = 1..9 enumerate every d-dim
   subspace W of F2^9 (reduced row echelon form, pivots = leading bits), let L(W) = lines whose label lies in W and
   T(W) = number of orthogonal 8-frames meeting L(W); floor(d) = 8 * F * d / max_W T(W). Input p112.txt. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define NW 20
static int n, F, lab[112], lineof[512];
static uint64_t fm[112][NW];
static long cnt; static int Tmax; static int bestbasis[9]; static int bestd;
static int popc(const uint64_t *m){int s=0;for(int i=0;i<NW;i++)s+=__builtin_popcountll(m[i]);return s;}
static void eval(int *rows, int d){
  uint64_t m[NW]; memset(m,0,sizeof m);
  int els = 1<<d;
  for(int s=1;s<els;s++){int x=0;for(int r=0;r<d;r++)if(s>>r&1)x^=rows[r];int t=lineof[x];if(t>=0)for(int i=0;i<NW;i++)m[i]|=fm[t][i];}
  int T=popc(m); cnt++;
  if(T>Tmax){Tmax=T;for(int r=0;r<d;r++)bestbasis[r]=rows[r];bestd=d;}
}
/* RREF: choose pivot bits p0>p1>...; row r has bit p_r set, no other pivot bits, and free bits only below p_r */
static int piv[9], rows[9], d0;
static void fill(int r, int d){
  if(r==d){eval(rows,d);return;}
  int freebits[9], nf=0;
  for(int b=0;b<piv[r];b++){int isp=0;for(int q=0;q<d;q++)if(piv[q]==b)isp=1;if(!isp)freebits[nf++]=b;}
  for(int mask=0;mask<(1<<nf);mask++){int v=1<<piv[r];for(int i=0;i<nf;i++)if(mask>>i&1)v|=1<<freebits[i];rows[r]=v;fill(r+1,d);}
}
static void pivots(int r, int d, int below){
  if(r==d){fill(0,d);return;}
  for(int p=below-1;p>=d-r-1;p--){piv[r]=p;pivots(r+1,d,p);}
}
int main(int argc,char**argv){
  FILE*f=fopen(argc>1?argv[1]:"p112.txt","r"); if(!f){perror("open");return 1;}
  if(fscanf(f,"%d %d",&n,&F)!=2)return 1;
  memset(lineof,-1,sizeof lineof);
  for(int i=0;i<n;i++){if(fscanf(f,"%d",&lab[i])!=1)return 1;lineof[lab[i]]=i;}
  memset(fm,0,sizeof fm);
  for(int c=0;c<F;c++)for(int j=0;j<8;j++){int t;if(fscanf(f,"%d",&t)!=1)return 1;fm[t][c/64]|=1ULL<<(c%64);}
  double best=1e9;
  for(int d=1;d<=9;d++){
    cnt=0;Tmax=0; pivots(0,d,9);
    double fl=8.0*F*d/Tmax; if(fl<best)best=fl;
    printf("d=%d subspaces=%ld Tmax=%d of %d floor(d)=%.4f basis:",d,cnt,Tmax,F,fl);
    for(int r=0;r<d;r++)printf(" %d",bestbasis[r]); printf("\n"); fflush(stdout);
  }
  printf("one-sided clique floor = %.4f per call\n",best);
  return 0;
}
