/* sa4.c -- (= sa3 + exhaustive 3-row and 3-column flat tests in the objective) annealing over call sequences on 16 wires, objective = violations of NECESSARY conditions
   for computing H^(x)4 up to monomials with g calls of 2x2 maps monomially equal to H:
     for every row set I in the future-block family and column set J in the past-block family
     (blocks of the components of calls t..g and of calls 1..s, over all s,t), the GENERIC rank of
     W[I,J] must be >= ceil(|I||J|/16)  (every a x b submatrix of a Hadamard matrix of order 16
     has rank >= ab/16; monomials do not change ranks);  plus density (1x1 tests).
   Generic rank is certified from below by exact rank mod p = 2^31-1 at random gate parameters
   (a specialisation can only lower the rank; a shortfall can be bad luck with prob ~ deg/p).
   Usage: sa2 g iters restarts seed [T0 T1].  Memory: a few KB. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define N 16
#define MAXG 40
#define P 2147483647ULL
typedef unsigned long long u64; typedef unsigned short u16;
static int G; static int cu[MAXG], cv[MAXG]; static u64 gp[MAXG][4];
static u64 rs;
static inline u64 rnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static inline double urand(void){ return (rnd()>>11)*(1.0/9007199254740992.0); }
static int pc(unsigned x){ return __builtin_popcount(x); }
static inline u64 mm(u64 a,u64 b){ u64 x=a*b; x=(x&P)+(x>>31); x=(x&P)+(x>>31); return x>=P?x-P:x; }
static u64 pw(u64 a,u64 e){ u64 r=1; while(e){ if(e&1) r=mm(r,a); a=mm(a,a); e>>=1;} return r; }
static void newgate(int k){ u64 x1=rnd()%(P-1)+1,x2=rnd()%(P-1)+1,y1=rnd()%(P-1)+1,y2=rnd()%(P-1)+1;
  gp[k][0]=mm(x1,y1); gp[k][1]=mm(x1,y2); gp[k][2]=mm(x2,y1); gp[k][3]=(P-mm(x2,y2))%P; }
static u64 W[N][N];
static void buildW(void){
  for(int i=0;i<N;i++) for(int j=0;j<N;j++) W[i][j]=(i==j);
  for(int k=0;k<G;k++){ int u=cu[k],v=cv[k];
    for(int j=0;j<N;j++){ u64 a=W[u][j],b=W[v][j];
      u64 s1=mm(gp[k][0],a)+mm(gp[k][1],b), s2=mm(gp[k][2],a)+mm(gp[k][3],b); W[u][j]=s1>=P?s1-P:s1; W[v][j]=s2>=P?s2-P:s2; } }
}
static int rankmod(u16 I,u16 J){
  u64 M[N][N]; int r=0,c=0; int ri[N],cj[N],nr=0,nc=0;
  for(int i=0;i<N;i++) if(I>>i&1) ri[nr++]=i; for(int j=0;j<N;j++) if(J>>j&1) cj[nc++]=j;
  for(int a=0;a<nr;a++) for(int b=0;b<nc;b++) M[a][b]=W[ri[a]][cj[b]];
  for(c=0;c<nc && r<nr;c++){
    int piv=-1; for(int a=r;a<nr;a++) if(M[a][c]){piv=a;break;}
    if(piv<0) continue;
    if(piv!=r) for(int b=0;b<nc;b++){u64 t=M[r][b];M[r][b]=M[piv][b];M[piv][b]=t;}
    u64 pv=M[r][c];
    for(int a=r+1;a<nr;a++) if(M[a][c]){ u64 f=M[a][c]; for(int b=c;b<nc;b++) { u64 t=mm(M[a][b],pv)+P-mm(f,M[r][b]); M[a][b]=t>=P?t-P:t; } }
    r++;
  }
  return r;
}
static int verbose=0;
static u16 fam[2][8*MAXG]; static int nfam[2];
static void addfam0(int which,u16 m){ for(int i=0;i<nfam[which];i++) if(fam[which][i]==m) return; fam[which][nfam[which]++]=m; }
static void addfam(int which,u16 m){ addfam0(which,m); if(m!=0xFFFF) addfam0(which,(u16)~m); }
static int twotest(int rows){ /* rows=1: pairs of rows; 0: pairs of columns. returns violations */
  int viol=0;
  for(int a=0;a<N;a++) for(int b=a+1;b<N;b++){
    u64 x[N],y[N]; for(int j=0;j<N;j++){ x[j]= rows? W[a][j]:W[j][a]; y[j]= rows? W[b][j]:W[j][b]; }
    int zero=0, best=0;
    for(int j=0;j<N;j++){ if(!x[j]&&!y[j]){zero++;continue;} int c=0; for(int k=0;k<N;k++){ if(!x[k]&&!y[k]) continue; if(mm(x[j],y[k])==mm(x[k],y[j])) c++; } if(c>best) best=c; }
    int r1=best+zero; if(r1>8){ viol+=1; if(verbose) printf("  2-%s %d,%d: rank<=1 on %d >8\n",rows?"row":"col",a,b,r1);} }
  return viol;
}

/* all 3-row (rows=1) or 3-column sets: largest rank-1 flat must be <= 5 columns, largest rank-2 flat <= 10 */
static int threetest(int rows){
  int viol=0;
  for(int a=0;a<N;a++) for(int b=a+1;b<N;b++) for(int c=b+1;c<N;c++){
    u64 V[N][3]; int nz[N], nn=0, z=0;
    for(int j=0;j<N;j++){ V[j][0]= rows? W[a][j]:W[j][a]; V[j][1]= rows? W[b][j]:W[j][b]; V[j][2]= rows? W[c][j]:W[j][c];
      if(!V[j][0]&&!V[j][1]&&!V[j][2]) z++; else nz[nn++]=j; }
    int best1=0, best2=0;
    for(int p=0;p<nn;p++){ int cnt1=0; u64 *x=V[nz[p]];
      for(int q=0;q<nn;q++){ u64 *y=V[nz[q]];
        /* parallel iff cross product zero */
        if(mm(x[1],y[2])==mm(x[2],y[1]) && mm(x[0],y[2])==mm(x[2],y[0]) && mm(x[0],y[1])==mm(x[1],y[0])) cnt1++; }
      if(cnt1>best1) best1=cnt1;
      for(int q=p+1;q<nn;q++){ u64 *y=V[nz[q]];
        u64 n0=(mm(x[1],y[2])+P-mm(x[2],y[1]))%P, n1=(mm(x[2],y[0])+P-mm(x[0],y[2]))%P, n2=(mm(x[0],y[1])+P-mm(x[1],y[0]))%P;
        if(!n0&&!n1&&!n2) continue;
        int cnt2=0; for(int r=0;r<nn;r++){ u64 *w=V[nz[r]]; if((mm(n0,w[0])+mm(n1,w[1])+mm(n2,w[2]))%P==0) cnt2++; }
        if(cnt2>best2) best2=cnt2; } }
    if(best1+z>5){ viol++; if(verbose) printf("  3-%s %d,%d,%d: %d columns of rank<=1\n",rows?"row":"col",a,b,c,best1+z); }
    if(best2+z>10){ viol++; if(verbose) printf("  3-%s %d,%d,%d: %d columns of rank<=2\n",rows?"row":"col",a,b,c,best2+z); }
  }
  return viol;
}
static void families(void){
  nfam[0]=nfam[1]=0;
  /* past blocks: forward merges */
  u16 blk[N]; for(int w=0;w<N;w++){ blk[w]=1u<<w; addfam(0,blk[w]); }
  for(int k=0;k<G;k++){ int u=cu[k],v=cv[k]; if(blk[u]!=blk[v]){ u16 m=blk[u]|blk[v]; for(int w=0;w<N;w++) if(m>>w&1) blk[w]=m; addfam(0,m);} }
  for(int w=0;w<N;w++) blk[w]=1u<<w;
  for(int w=0;w<N;w++) addfam(1,blk[w]);
  for(int k=G-1;k>=0;k--){ int u=cu[k],v=cv[k]; if(blk[u]!=blk[v]){ u16 m=blk[u]|blk[v]; for(int w=0;w<N;w++) if(m>>w&1) blk[w]=m; addfam(1,m);} }
}
static double evaluate(void){
  buildW(); families();
  int viol=0;
  for(int i=0;i<N;i++) for(int j=0;j<N;j++) if(!W[i][j]) viol+=1;
  viol+=twotest(1)+twotest(0); if(viol==0) viol+=threetest(1)+threetest(0);
  for(int a=0;a<nfam[1];a++) for(int b=0;b<nfam[0];b++){
    u16 I=fam[1][a], J=fam[0][b]; int na=pc(I), nb=pc(J); if(na*nb<=16) continue; /* covered by density */
    int need=(na*nb+15)/16; int r=rankmod(I,J);
    if(r<need){ viol+=need-r; if(verbose) printf("  I=%04x(%d) J=%04x(%d) rank %d < %d\n",I,na,J,nb,r,need); }
  }
  return viol;
}
int main(int argc,char**argv){
  G=atoi(argv[1]); long iters=atol(argv[2]); int restarts=atoi(argv[3]); rs=strtoull(argv[4],0,10)*2654435761ULL+12345;
  double T0= argc>5? atof(argv[5]):1.5, T1= argc>6? atof(argv[6]):0.05;
  if(iters==0){ /* verify mode: read calls from stdin */
    for(int k=0;k<G;k++){ if(scanf("%d %d",&cu[k],&cv[k])!=2) return 1; newgate(k);} verbose=1;
    printf("violations %.0f\n",evaluate()); return 0; }
  double best_all=1e9;
  for(int r=0;r<restarts;r++){
    if(getenv("SEED_STD")){ /* start from the standard 32-call butterfly with 32-G random calls deleted */
      int su[32],sv[32],n=0; for(int b=0;b<4;b++) for(int x=0;x<16;x++) if(!((x>>b)&1)){ su[n]=x; sv[n]=x|(1<<b); n++; }
      int keep[32]; for(int i=0;i<32;i++) keep[i]=1; for(int del=0;del<32-G;){ int i=rnd()%32; if(keep[i]){keep[i]=0;del++;} }
      int k=0; for(int i=0;i<32 && k<G;i++) if(keep[i]){ cu[k]=su[i]; cv[k]=sv[i]; newgate(k); k++; }
    } else
    for(int k=0;k<G;k++){ int u=rnd()%N, v; do v=rnd()%N; while(v==u); cu[k]=u; cv[k]=v; newgate(k); }
    double cur=evaluate(), best=cur; int bu[MAXG],bv[MAXG]; memcpy(bu,cu,sizeof cu); memcpy(bv,cv,sizeof cv);
    clock_t c0=clock();
    for(long it=0;it<iters;it++){
      if((it&63)==0 && (double)(clock()-c0)/CLOCKS_PER_SEC>120.0){ printf("  (time limit at it=%ld)\n",it); break; }
      double fr=(double)it/iters, fe=(double)(clock()-c0)/CLOCKS_PER_SEC/120.0; if(fe>fr) fr=fe;
      double T=T0*pow(T1/T0,fr);
      int k=rnd()%G, mv=rnd()%3; int su=cu[k], sv=cv[k], k2=-1;
      if(mv==0){ int u=rnd()%N,v; do v=rnd()%N; while(v==u); cu[k]=u; cv[k]=v; }
      else if(mv==1){ if(rnd()&1){ int u; do u=rnd()%N; while(u==cv[k]); cu[k]=u; } else { int v; do v=rnd()%N; while(v==cu[k]); cv[k]=v; } }
      else { k2=rnd()%G; int tu=cu[k2],tv=cv[k2]; cu[k2]=cu[k]; cv[k2]=cv[k]; cu[k]=tu; cv[k]=tv; }
      double nw=evaluate();
      if(nw<=cur || urand()<exp((cur-nw)/T)){ cur=nw; if(cur<best){best=cur; memcpy(bu,cu,sizeof cu); memcpy(bv,cv,sizeof cv);}
        if(cur==0) break; }
      else { if(mv==2){ int tu=cu[k2],tv=cv[k2]; cu[k2]=cu[k]; cv[k2]=cv[k]; cu[k]=tu; cv[k]=tv; } else { cu[k]=su; cv[k]=sv; } }
    }
    if(best<best_all) best_all=best;
    printf("g=%d restart %d best %.0f :",G,r,best); for(int q=0;q<G;q++) printf(" %d-%d",bu[q],bv[q]); printf("\n"); fflush(stdout);
  }
  printf("g=%d overall best %.0f\n",G,best_all);
  return 0;
}
