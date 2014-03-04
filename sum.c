#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void print_matrix(double *matrix, int m) {
  int i;
  for (i=0; i<m; i++)
    printf("%g\n",matrix[i]);
}

main(int argc, char*argv[])
{
  FILE *fp_in;
  int n;
  long int N,i,ni=0;
  double *v,sum=0;

  if(argc!=2) {
    printf("Sintaxe:\n\n%s ficheiro\n\n",argv[0]);
    exit(1);
  }

  if((fp_in=fopen(argv[1],"r"))==NULL) {
    printf("Impossivel abrir o ficheiro %s\n",argv[1]);
    exit(2);
  }

  fseek(fp_in,0,SEEK_END);
  N = ftell(fp_in)/(sizeof(double));
  rewind(fp_in);
  
  v = (double *)calloc(N,sizeof(double));
  if(v==NULL) puts("Problemas na Alocacao da memoria");
  
  for(i=0;i<N;i++) {
    n=fscanf(fp_in,"%lf",&(v[i]));
    if(n>0) ni+=n;
  }
  
  printf("ni=%ld, N=%ld\n",ni,N);

  printf(" --input-- \n");
  print_matrix(v,N);

  if(ni<N) {
    N=ni;
    v = (double *)realloc(v,N*sizeof(double));
    if(v==NULL) puts("Problemas na Re-Alocacao da memoria");
  }

  printf("\n --output-- \n");  
  print_matrix(v,N);
  printf("\n");

  for(i=0;i<N;i++)
    sum+=v[i];
  
  printf("soma total = %g\n",sum);
  
  fclose(fp_in);
  free(v);
}
