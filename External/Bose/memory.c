/*  Use numerical recipes in C styled memory management methods */

/*#include <malloc.h>*/
#include <stdlib.h>
#include <stdio.h>

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
  int i;
  int **m;

  m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
  if (!m){
    fprintf(stderr,"Unable to allocate %d int* s.\n",nrh-nrl+1);
    return m;
  }
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
    if (!m[i]) {
      fprintf(stderr,"Unable to allocate %d'th row in an integer matrix.\n",
	      i-nrl+1);
      return NULL;
    }
    m[i] -= ncl;
  }
  return m;
}

int *ivector(nl,nh)
int nl,nh;
{
  int *v;

  v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int) );
  if (!v){
    fprintf(stderr,"Unable to allocate %d int s.\n",nh-nl+1);
    return v;
  }
  return(v-nl);
}

/* ARGSUSED */
void free_imatrix(m,nrl,nrh,ncl,nch)

int nrl,nrh,ncl,nch;
int **m;

{
int i;

for(  i=nrh;i>=nrl;i--  )free( (char*) (m[i]+ncl) );
free((char*) (m+nrl));
}

/* ARGSUSED */
void free_ivector(v,nl,nh)
int *v;
int nl,nh;
{
  free((char*) (v+nl));
}



double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
  int      i;
  double **m;

  m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
  if (!m){
    fprintf(stderr,"Unable to allocate %d double* s.\n",nrh-nrl+1);
    return m;
  }
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    if (!m[i]) {
      fprintf(stderr,"Unable to allocate %d'th row in an double matrix.\n",
	      i-nrl+1);
      return NULL;
    }
    m[i] -= ncl;
  }
  return m;
}

double *dvector(nl,nh)
int nl,nh;
{
  double *v;

  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double) );
  if (!v){
    fprintf(stderr,"Unable to allocate %d double s.\n",nh-nl+1);
    return v;
  }
  return(v-nl);
}

/* ARGSUSED */
void free_dmatrix(m,nrl,nrh,ncl,nch)

int nrl,nrh,ncl,nch;
double **m;

{
int i;

for(  i=nrh;i>=nrl;i--  )free( (char*) (m[i]+ncl) );
free((char*) (m+nrl));
}

/* ARGSUSED */
void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
  free((char*) (v+nl));
}










/*  Add rows to a matrix:

      imat            pointer to that matrix
      oldrowsize      highest index+1 of old rows
      newrowsize      highest index+1 of new rows
      colsize         number of cols

*/

int grow_imatrix_byrows( imat, oldrowsize, newrowsize, colsize )
int ***imat, oldrowsize, newrowsize, colsize;
{
  int i;

  imat[0] = (int **) realloc(imat[0],
		  (unsigned) (newrowsize)*sizeof(int*));
  if (!(imat[0])){
    fprintf(stderr,"Unable to reallocate %d int* s.\n",newrowsize);
    return 0;
  }
  for(i=oldrowsize; i<newrowsize; i++) {
    imat[0][i] = (int *) malloc((unsigned) (colsize)*sizeof(int));
    if (!imat[0][i]) {
      fprintf(stderr,"Unable to reallocate %d'th row in an integer matrix.\n",
	      i);
      return 0;
    }
  }
  return 1;
}
