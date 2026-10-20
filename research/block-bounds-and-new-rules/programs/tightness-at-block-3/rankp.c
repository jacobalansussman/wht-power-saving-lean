/* rankp.c: rank of a matrix modulo a prime p < 2^31, in place (destroys the input).
   Row-major int64 entries already reduced to [0,p). Forward elimination only. Barrett reduction. */
#include <stdint.h>
static inline uint64_t redp(uint64_t x, uint64_t P, uint64_t mu){
  uint64_t q=(uint64_t)(((__uint128_t)x*mu)>>64); uint64_t r=x-q*P; while(r>=P) r-=P; return r;}
static uint64_t pw(uint64_t a, uint64_t e, uint64_t p){uint64_t r=1;a%=p;while(e){if(e&1)r=r*a%p;a=a*a%p;e>>=1;}return r;}
long rank_modp(int64_t *A, long rows, long cols, int64_t p){
  long r=0; uint64_t P=(uint64_t)p; uint64_t mu=(uint64_t)(((__uint128_t)1<<64)/P);
  for(long c=0;c<cols && r<rows;c++){
    long piv=-1;
    for(long i=r;i<rows;i++) if(A[i*cols+c]){piv=i;break;}
    if(piv<0) continue;
    if(piv!=r){ for(long j=c;j<cols;j++){int64_t t=A[r*cols+j];A[r*cols+j]=A[piv*cols+j];A[piv*cols+j]=t;} }
    uint64_t inv=pw((uint64_t)A[r*cols+c],P-2,P);
    uint64_t *R=(uint64_t*)(A+r*cols);
    for(long j=c;j<cols;j++) R[j]=redp(R[j]*inv,P,mu);
    for(long i=r+1;i<rows;i++){
      uint64_t *X=(uint64_t*)(A+i*cols); uint64_t fct=X[c];
      if(!fct) continue;
      uint64_t nf=P-fct;
      for(long j=c;j<cols;j++) X[j]=redp(X[j]+nf*R[j],P,mu);
    }
    r++;
  }
  return r;
}
