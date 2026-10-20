/* ctr-window / verify_mu.c (own code). Exact integer re-check of a weighted covering bound.
   Input: pool file (n P F / labels / P sets of 8 lines) and a weight file (P doubles, mu). Weights are rounded to
   integers c_S = llround(mu_S * 2^30) (any nonnegative weights give a valid bound). Every set needs >= 8 terms meeting
   it, so for every single-layer decoder (one-sided cost = dim of gather span W; term meets S iff S has a line with
   label in W):  cost * max_W [H(W)/dim W] >= 8 * C,  C = sum c_S,  H(W) = sum of c_S over S meeting L(W).
   All sums are exact int64; the maximum of H/d is found by exact cross-multiplication over EVERY subspace W of F2^9.
   Output: the floor 8 C d_best / H_best as a fraction and decimal. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define MAXP 8192
#define PW (MAXP/64)
static int n,P,F,lab[112],lineof[512],pw; static uint64_t sm[112][PW]; static int64_t c[MAXP]; static int64_t tab[MAXP/8][256];
static int64_t Hb=0; static int db=1; static int bb[9]; static long nsub=0;
static int piv[9],rows[9];
static void evalW(int d){ uint64_t l0=0,l1=0; int els=1<<d;
  for(int s=1;s<els;s++){int x=0;for(int r=0;r<d;r++)if(s>>r&1)x^=rows[r];int t=lineof[x];if(t>=0){if(t<64)l0|=1ULL<<t;else l1|=1ULL<<(t-64);}}
  nsub++; if(!l0&&!l1)return; uint64_t m[PW]; memset(m,0,sizeof(uint64_t)*pw);
  for(int t=0;t<n;t++) if((t<64?(l0>>t):(l1>>(t-64)))&1) for(int i=0;i<pw;i++) m[i]|=sm[t][i];
  int64_t H=0; int nb=(P+7)/8; for(int b=0;b<nb;b++) H+=tab[b][(m[b/8]>>((b%8)*8))&255];
  if((__int128)H*db > (__int128)Hb*d){Hb=H;db=d;for(int q=0;q<d;q++)bb[q]=rows[q];} }
static void fill(int r,int d){ if(r==d){evalW(d);return;} int fb[9],nf=0; for(int b=0;b<piv[r];b++){int isp=0;for(int q=0;q<d;q++)if(piv[q]==b)isp=1;if(!isp)fb[nf++]=b;}
  for(int mask=0;mask<(1<<nf);mask++){int v=1<<piv[r];for(int i=0;i<nf;i++)if(mask>>i&1)v|=1<<fb[i];rows[r]=v;fill(r+1,d);} }
static void pivots(int r,int d,int below){ if(r==d){fill(0,d);return;} for(int p=below-1;p>=d-r-1;p--){piv[r]=p;pivots(r+1,d,p);} }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); FILE*g=fopen(argv[2],"r"); if(!f||!g){perror("open");return 1;}
  if(fscanf(f,"%d %d %d",&n,&P,&F)!=3)return 1; pw=(P+63)/64; memset(lineof,-1,sizeof lineof);
  for(int i=0;i<n;i++){if(fscanf(f,"%d",&lab[i])!=1)return 1; if(lineof[lab[i]]!=-1){printf("dup label\n");return 1;} lineof[lab[i]]=i;}
  for(int s=0;s<P;s++)for(int j=0;j<8;j++){int t;if(fscanf(f,"%d",&t)!=1)return 1;sm[t][s/64]|=1ULL<<(s%64);}
  int64_t C=0; for(int s=0;s<P;s++){double x; if(fscanf(g,"%lf",&x)!=1)return 1; c[s]=llround(x*1073741824.0); if(c[s]<0)return 1; C+=c[s];}
  int nb=(P+7)/8; for(int b=0;b<nb;b++){tab[b][0]=0;for(int v=1;v<256;v++){int low=__builtin_ctz(v);int s=b*8+low;tab[b][v]=tab[b][v&(v-1)]+(s<P?c[s]:0);}}
  for(int d=1;d<=9;d++) pivots(0,d,9);
  /* floor = 8*C*db/Hb */
  __int128 num=(__int128)8*C*db; int64_t den=Hb;
  printf("subspaces checked %ld (all nonzero subspaces of F2^9: 8,283,457); C=%lld; worst W dim %d, H=%lld\n",nsub,(long long)C,db,(long long)Hb);
  printf("floor = %lld*8*%d/%lld = %.6f per call (integer cost >= %lld)\n",(long long)C,db,(long long)Hb,(double)num/den,(long long)((num+den-1)/den));
  printf("worst W basis:"); for(int q=0;q<db;q++)printf(" %d",bb[q]); printf("\n");
  return 0;
}
