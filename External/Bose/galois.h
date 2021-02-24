/* Declarations for dynamic allocation */

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

