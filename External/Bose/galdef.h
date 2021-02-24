/* Definitions for Galois Field code */

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
