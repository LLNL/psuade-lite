/* Declarations for dynamic allocation */

#ifndef __GALOIS__
#define __GALOIS__

int *ivector(), **imatrix();

/* Definitions for Galois Fields */

struct GF {
  int n,p,q;
  int *xton;
  int **plus;
  int **times;
  int *inv;
  int *neg;
  int *root;
  int **poly;
};

int bush(struct GF *gf, int **A, int str, int ncol);
int GF_getfield( int q, struct GF *gf);
int bosecheck( int q, int ncol);
int ipow( int a, int b);
int GF_ready(struct GF *gf, int  n,int p,int *xton);
int isprime(int);
void primepow(int q,int *p,int *n,int *isit);
void free_ivector(int *v, int nl, int nh);
int  akodd(struct GF *gf, int *kay, int *b, int *c, int *k);
int  akeven(struct GF *gf, int *kay, int *b, int *c, int *k);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
int  bose( struct GF *gf, int **A, int ncol);

#endif

