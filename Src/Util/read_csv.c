#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LBUFSIZ 10000

int read_csv(char *filename, int *nSamples, int *nInps, double **XX,
             int *nOuts, char ***inames, char ***onames)
{
  FILE *myfile, *outfile;
  char line1[LBUFSIZ];
  char line2[LBUFSIZ];
  char line3[LBUFSIZ]; /* Each field in the line */
  char cword[1001];
  char *stptr, **ionames, **tcarray;
  int flag = 0, counter=0, first=1, csize, nclines, ncols, nInOuts;
  int idx = 0, xsize, ii, jj, kk, leng, nInputs, nOutputs, ipos;
  int lcount = 0; /* Cell Separator */
  int numLines=0, isFloat;
  double ddata;
  double *X, *Xt;

  /* check that the file exists */
  if (!(myfile = fopen(filename, "r")))
  {
    printf("read_csv ERROR: Could not open file for reading\n");
    return -1;
  }
  /* read the first line, if there is one */
  if (fgets(line1,sizeof line1,myfile) == NULL)
  {
    fclose(myfile);
    return -1;
  }
  /* see if the line has any non-numeric character and non-comma in it*/
  leng = strlen(line1);
  for (kk = 0; kk < leng; kk++)
    if (((line1[kk]-'0') < 0 || (line1[kk]-'0') > 9) && 
          line1[kk] != ',' && line1[kk] != ' ' &&
          line1[kk] != 'E' && line1[kk] != '+' &&
          line1[kk] != 'e' && line1[kk] != '.' &&
          line1[kk] != '-' && line1[kk] != 10  && line1[kk] != 13)
       break;
  /* if the line has any non-numeric character */
  (*nInps) = 0;
  (*nOuts) = 0;
  (*inames) = NULL;
  (*onames) = NULL;
  nclines = nInputs = nOutputs = 0;
  ionames = NULL;
  csize   = 0;
  if (kk < leng)
  {
    printf("READ_CSV: the first line contains non-numerics.\n");
    printf("Use the first line as the variable definition line? (y or n) ");
    scanf("%s", cword);
    if (cword[0] == 'y') nclines++;
    fgets(cword,10,stdin);
  }
  if (kk < leng && nclines > 0)
  {
    /* read in the variable names or input/output */
    csize = 10;
    ionames = (char **) malloc(csize * sizeof(char*));
    for (ii = 0; ii < csize; ii++) 
       ionames[ii] = (char *) malloc(1000 * sizeof(char));
    ncols = ipos = kk = flag = 0;
    while (kk < leng)
    {
      if (line1[kk] != ' ' && line1[kk] != ',' && line1[kk] != 10 &&
          line1[kk] != 13)
        ionames[ncols][ipos++] = line1[kk];
      else if (line1[kk] == ',' || line1[kk] == 10 || line1[kk] == 13)
      {
        ionames[ncols][ipos] = '\0';
        if (!strcmp(ionames[ncols],"INPUT") || 
            !strcmp(ionames[ncols],"input") || 
            !strcmp(ionames[ncols],"Input")) nInputs++;
        else if (!strcmp(ionames[ncols],"OUTPUT") || 
                 !strcmp(ionames[ncols],"output") || 
                 !strcmp(ionames[ncols],"Output")) nOutputs++;
        else flag = 1;

        ncols++; ipos = 0;
        if (ncols >= csize)
        {
          tcarray = ionames;
          ionames = (char **) malloc((csize+10) * sizeof(char*));
          for (ii = 0; ii < ncols; ii++) ionames[ii] = tcarray[ii];
            for (ii = ncols; ii < csize+10; ii++) 
              ionames[ii] = (char *) malloc(1000 * sizeof(char));
          csize += 10;
          free(tcarray);
        }
      }
      kk++;
      if (line1[kk] == 10 || line1[kk] == 13)
      {
        for (ii = 0; ii < ipos; ii++) if (ionames[ncols][ii] != ' ') break;
        if (ipos == 0 || (ipos != 0 && ii == ipos))  break;
      }
    }
    if (ipos > 0) ncols++;
    nclines = 1;
  }  
  if (flag == 0)
  {
    if (ionames != NULL)
    {
      for (ii = 0; ii < csize; ii++) free(ionames[ii]);
      free(ionames);
      ionames = NULL;
    }
  }
#if 0
  /* read the second line, if the first line is non-numeric */
  if (nclines > 0 && flag == 0)
  {
    if (fgets(line1,sizeof line1,myfile) == NULL) 
    {
      fclose(myfile);
      for (ii = 0; ii < csize; ii++) 
        if (ionames[ii] != NULL) free(ionames[ii]);
      return -1;
    }
    leng = strlen(line1);
    /* see if the line has any non-numeric character and non-comma in it*/
    leng = strlen(line1);
    for (kk = 0; kk < leng; kk++)
      if (((line1[kk]-'0') < 0 || (line1[kk]-'0') > 9) && 
            line1[kk] != ',' && line1[kk] != ' ' &&
            line1[kk] != 'E' && line1[kk] != '+' &&
            line1[kk] != 'e' && line1[kk] != '.' &&
            line1[kk] != '-' && line1[kk] != 10  && line1[kk] != 13)
        break;
    /* INPUT/OUTPUT line has been read, so read in variable names now */
    if (kk < leng)
    {
      /* read in the variable names */
      csize = 10;
      ionames = (char **) malloc(csize * sizeof(char*));
      for (ii = 0; ii < csize; ii++) 
        ionames[ii] = (char *) malloc(1000 * sizeof(char));
      ncols = ipos = kk = flag = 0;
      while (kk < leng)
      {
        if (line1[kk] != ' ' && line1[kk] != ',' && line1[kk] != 10 &&
            line1[kk] != 13)
          ionames[ncols][ipos++] = line1[kk];
        else if (line1[kk] == ',' || line1[kk] == 10 || line1[kk] == 13)
        {
          ionames[ncols][ipos] = '\0';
          ncols++; ipos = 0;
          if (ncols >= csize)
          {
            tcarray = ionames;
            ionames = (char **) malloc((csize+10) * sizeof(char*));
            for (ii = 0; ii < ncols; ii++) ionames[ii] = tcarray[ii];
            for (ii = ncols; ii < csize+10; ii++) 
              ionames[ii] = (char *) malloc(1000 * sizeof(char));
            csize += 10;
            free(tcarray);
          }
        }
        kk++;
        if (line1[kk] == 10 || line1[kk] == 13)
        {
          for (ii = 0; ii < ipos; ii++) if (ionames[ncols][ii] != ' ') break;
          if (ipos == 0 || (ipos != 0 && ii == ipos))  break;
        }
      }
      if (ipos > 0) ncols++;
      if (ncols != nInputs+nOutputs)
      {
        printf("read_csr ERROR (1).\n");
        if (ionames != NULL)
        {
          for (ii = 0; ii < csize; ii++) free(ionames[ii]);
          free(ionames);
          ionames = NULL;
        }
        return -1;
      }
    }
    nclines++;
  }
#endif
  fclose(myfile);

  /* read the file again */
  myfile = fopen(filename, "r");
  for (kk = 0; kk < nclines; kk++) fgets(line1,sizeof line1,myfile);
  xsize = 1000;
  X = (double *) malloc(xsize * sizeof(double));
  /* Get a line from file */
  numLines = 0;
  while (fgets(line1,sizeof line1,myfile) != NULL)
  { 
    if (line1[0] == '#') continue;
    lcount = 0;
    strcpy(line2,line1);
    stptr = line2;
    numLines++;

    /* start going character by character thro the line */
    while (*stptr != '\0' && *stptr != 10)
    { 
      lcount++;
      /* If field begins with " */
      if (*stptr == '"')
      {
        int flag = 0;
        idx = 0;
        while (flag == 0)
        {
          stptr++;
          /* Find corresponding closing " */
          while (*stptr != '"')
          { 
            line3[idx] = *stptr;
            idx++;
            stptr++;
          }
          stptr++;
          if (*stptr != '\0' && *stptr == ',')
          {
            line3[idx] = '\0';
            flag = 1;
          }
          else if (*stptr != '\0' && *stptr == '"')
          { 
            line3[idx] = *stptr;
            idx++;
          }
          else
          {
            line3[idx] = '\0';
            flag = 1;
          }
        }
        X[counter++] = 0;
        if (numLines <= 1)
          printf("INFO: Non-numeric data read - set the field %d to 0.\n",
                 counter);
      }
      else
      { 
        idx = 0;
        while (*stptr != '\0' && *stptr != ',')
        {  
          line3[idx] = *stptr;
          idx++;
          stptr++;
        }
        line3[idx] = '\0';
        isFloat = sscanf(line3, "%lg", &ddata);
        if (isFloat != 1) ddata = 0;
        X[counter++] = ddata;
        if (counter >= xsize)
        {
          Xt = X;
          xsize += 10000;
          X = (double *) malloc(xsize*sizeof(double));
          for (kk = 0; kk < counter; kk++) X[kk] = Xt[kk];
          free(Xt);
        }
      }
      if (*stptr != '\0' && *stptr == ',')
        stptr++;
    }
    if (first == 1)
    {
      nInOuts = counter;
      first = 0;
    }
  }
  fclose(myfile);
  (*XX) = X;
  (*nSamples) = counter / nInOuts;
  if (nInputs == 0) (*nInps) = nInOuts;
  if (nInputs+nOutputs > 0)
  {
    if (nInputs > 0)
    {
      (*inames) = (char **) malloc(nInputs * sizeof(char *));
      for (ii = 0; ii < nInputs; ii++)
        (*inames)[ii] = ionames[ii]; 
    }
    if (nOutputs > 0)
    {
      (*onames) = (char **) malloc(nOutputs * sizeof(char *));
      for (ii = 0; ii < nOutputs; ii++)
        (*onames)[ii] = ionames[nInputs+ii]; 
    }
    free(ionames);
    (*nInps) = nInputs;
    (*nOuts) = nOutputs;
  }
  else (*inames) = ionames;
  return 0;
}

