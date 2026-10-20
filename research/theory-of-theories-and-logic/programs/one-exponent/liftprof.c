/* meta-spectrum liftprof.c: full distance profile of a lift family on the Lagrangian frames of Sp(2m,2).
   Input: a lift file of b56-methods/gen.py (header "NF M D DEG"; per frame: fixed-flag, M*D integers, DEG neighbours).
   For every unordered pair of frames {U,V}: d = rank_F2([B_U;B_V]) - M (frame distance), delta = rank_Q([B_U;B_V]) - M.
   rank_Q exactly by fraction-free (Bareiss) elimination with 128-bit products (minors fit in 64 bits: Hadamard bound
   for |entries| <= 4 and 2M <= 12 rows, as in b56-methods/anneal.c); rank_F2 by the same elimination on entries mod 2.
   Output: for each d, the number of unordered pairs and the sum of delta, plus the (d, delta) histogram.
   Ordered-pair averages abar_d = sum/count (delta is symmetric because every Y_U has dimension M).
   usage: liftprof FILE NTHREADS [ilo ihi]   (rows i in [ilo, ihi) of the upper triangle; default all)
          liftprof FILE 1 sample NPAIRS SEED  (uniform random pairs U != V; for an estimate only) */
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
typedef long long ll;
static int NF, M, D, DEG, RL; static int *R; static int *NB; static int ilo, ihi, NT;
static ll cnt[64][8][8];   /* per thread: [t][d][delta] */
static int rankq(const int *a, const int *b){
  ll A[12][13]; int rows=2*M;
  for(int i=0;i<M;i++)for(int c=0;c<D;c++){A[i][c]=a[D*i+c];A[M+i][c]=b[D*i+c];}
  int rk=0; ll prev=1;
  for(int c=0;c<D&&rk<rows;c++){
    int p=-1;for(int r=rk;r<rows;r++)if(A[r][c]){p=r;break;}
    if(p<0)continue;
    if(p!=rk)for(int k=0;k<D;k++){ll t=A[rk][k];A[rk][k]=A[p][k];A[p][k]=t;}
    for(int r=rk+1;r<rows;r++){
      for(int k=c+1;k<D;k++)A[r][k]=(ll)(((__int128)A[r][k]*A[rk][c]-(__int128)A[r][c]*A[rk][k])/prev);
      A[r][c]=0;
    }
    prev=A[rk][c];rk++;
  }
  return rk;
}
static unsigned *R2;  /* rows mod 2 as bitmasks */
static int rank2(const unsigned *a, const unsigned *b){
  unsigned A[12]; int rows=2*M; for(int i=0;i<M;i++){A[i]=a[i];A[M+i]=b[i];}
  int rk=0;
  for(int c=0;c<D&&rk<rows;c++){
    int p=-1;for(int r=rk;r<rows;r++)if(A[r]>>c&1){p=r;break;}
    if(p<0)continue;
    unsigned t=A[rk];A[rk]=A[p];A[p]=t;
    for(int r=0;r<rows;r++)if(r!=rk&&(A[r]>>c&1))A[r]^=A[rk];
    rk++;
  }
  return rk;
}
static void *work(void *arg){
  int t=(int)(long)arg;
  for(int i=ilo+t;i<ihi;i+=NT){
    for(int j=i+1;j<NF;j++){
      int d=rank2(R2+i*M,R2+j*M)-M, q=rankq(R+i*RL,R+j*RL)-M;
      if(q<d||q>M){fprintf(stderr,"IMPOSSIBLE pair %d %d d=%d delta=%d\n",i,j,d,q);exit(9);}
      cnt[t][d][q]++;
    }
  }
  return 0;
}
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: liftprof FILE NTHREADS [ilo ihi]\n");return 2;}
  FILE*f=fopen(argv[1],"r");if(!f)return 3;
  if(fscanf(f,"%d %d %d %d",&NF,&M,&D,&DEG)!=4)return 3;
  RL=M*D; NB=malloc(sizeof(int)*NF*DEG); R=malloc(sizeof(int)*NF*RL); R2=malloc(sizeof(unsigned)*NF*M);
  for(int i=0;i<NF;i++){int fl;if(fscanf(f,"%d",&fl)!=1)return 3;
    for(int c=0;c<RL;c++)if(fscanf(f,"%d",&R[i*RL+c])!=1)return 3;
    for(int k=0;k<DEG;k++){int x;if(fscanf(f,"%d",&x)!=1)return 3;NB[i*DEG+k]=x;}
    for(int r=0;r<M;r++){unsigned w=0;for(int c=0;c<D;c++)if(R[i*RL+r*D+c]&1)w|=1u<<c;R2[i*M+r]=w;}}
  fclose(f);
  NT=atoi(argv[2]); if(NT<1||NT>64)NT=1;
  ilo=0; ihi=NF;
  if(argc>=4 && argv[3][0]=='e'){ /* all edges (neighbour lists), exact */
    for(int i=0;i<NF;i++)for(int k=0;k<DEG;k++){int j=NB[i*DEG+k]; if(j<=i)continue;
      int d=rank2(R2+i*M,R2+j*M)-M, q=rankq(R+i*RL,R+j*RL)-M; cnt[0][d][q]++;}
    printf("file %s NF %d M %d D %d ALL EDGES of the neighbour lists (unordered)\n",argv[1],NF,M,D);
    for(int d=0;d<=M;d++){ll c=0,s=0;for(int q=0;q<=M;q++){c+=cnt[0][d][q];s+=cnt[0][d][q]*q;}
      if(!c)continue; printf("d=%d pairs %lld sumdelta %lld  hist",d,c,s);
      for(int q=0;q<=M;q++)if(cnt[0][d][q])printf(" %d:%lld",q,cnt[0][d][q]); printf("\n");}
    return 0; }
  if(argc>=6 && argv[3][0]=='s'){ /* sampling mode */
    long long ns=atoll(argv[4]); unsigned long long x=strtoull(argv[5],0,10)*0x9E3779B97F4A7C15ULL+12345;
    for(long long s=0;s<ns;s++){
      x^=x<<13;x^=x>>7;x^=x<<17; int i=(int)(x%NF); x^=x<<13;x^=x>>7;x^=x<<17; int j=(int)(x%NF);
      if(i==j){s--;continue;}
      int d=rank2(R2+i*M,R2+j*M)-M, q=rankq(R+i*RL,R+j*RL)-M; cnt[0][d][q]++; }
    printf("file %s NF %d M %d D %d SAMPLE of %lld uniform random ordered pairs (seed %s)\n",argv[1],NF,M,D,ns,argv[5]);
    for(int d=0;d<=M;d++){ll c=0,s=0,s2=0;for(int q=0;q<=M;q++){c+=cnt[0][d][q];s+=cnt[0][d][q]*q;s2+=cnt[0][d][q]*q*q;}
      if(!c)continue; double mu=(double)s/c, var=(double)s2/c-mu*mu;
      printf("d=%d pairs %lld sumdelta %lld mean %.6f stderr %.6f  hist",d,c,s,mu,c>1?__builtin_sqrt(var/(c-1)):0);
      for(int q=0;q<=M;q++)if(cnt[0][d][q])printf(" %d:%lld",q,cnt[0][d][q]); printf("\n");}
    return 0; }
  if(argc>=5){ilo=atoi(argv[3]);ihi=atoi(argv[4]);}
  /* sanity: every frame's rows have rank M over F2 and over Q */
  for(int i=0;i<NF;i++){ if(rank2(R2+i*M,R2+i*M)!=M){fprintf(stderr,"frame %d not of F2-rank M\n",i);return 4;} }
  pthread_t th[64]; for(long t=0;t<NT;t++)pthread_create(&th[t],0,work,(void*)t);
  for(int t=0;t<NT;t++)pthread_join(th[t],0);
  ll tot[8][8]={{0}}; for(int t=0;t<NT;t++)for(int a=0;a<8;a++)for(int b=0;b<8;b++)tot[a][b]+=cnt[t][a][b];
  printf("file %s NF %d M %d D %d rows [%d,%d) unordered pairs\n",argv[1],NF,M,D,ilo,ihi);
  for(int d=0;d<=M;d++){ll c=0,s=0;for(int q=0;q<=M;q++){c+=tot[d][q];s+=tot[d][q]*q;}
    if(!c)continue; printf("d=%d pairs %lld sumdelta %lld  hist",d,c,s);
    for(int q=0;q<=M;q++)if(tot[d][q])printf(" %d:%lld",q,tot[d][q]); printf("\n");}
  return 0;
}
