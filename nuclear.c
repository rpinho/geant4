#include <stdio.h>
#include <stdlib.h>

main(int argc, char*argv[])
{
  char f_out[]="nuclear_logbin_average.dat";
  FILE *fp_in1,*fp_in2,*fp_out;
  int i,j,m,n,M,N,ni=0,mi=0;
  double **v,**u;

  if(argc!=3) {
    printf("Sintaxe:\n\n%s ficheiro_total ficheiro_electronic\n\n",argv[0]);
    exit(1);
  }

  if((fp_in1=fopen(argv[1],"r"))==NULL) {
    printf("Impossivel abrir o ficheiro %s\n",argv[1]);
    exit(2);
  }

  if((fp_in2=fopen(argv[2],"r"))==NULL) {
    printf("Impossivel abrir o ficheiro %s\n",argv[2]);
    exit(3);
  }

  if((fp_out=fopen(f_out,"w"))==NULL) {
    printf("Impossivel criar o ficheiro %s\n",f_out);
    exit(4);
  }

  fseek(fp_in1,0,SEEK_END);
  N = ftell(fp_in1)/2/(sizeof(double)+1);
  rewind(fp_in1);
  
  v = (double **)calloc(N,sizeof(double*));
  if(v==NULL) puts("Problemas na Alocacao da memoria");
  
  for(i=0;i<N;i++) {
    v[i] = (double *)calloc(2,sizeof(double));
    for(j=0;j<2;j++) {
      n=fscanf(fp_in1,"%lf",&(v[i][j]));
      if(n>0) ni+=n;
    }
  }

  ni/=2;
  printf("ni=%ld, N=%ld\n",ni,N);

  //printf(" --input-- \n");
  //print_matrix(v,N,2);

  if(ni<N) {
    N=ni;
    v = (double **)realloc(v,N*sizeof(double*));
    if(v==NULL) puts("Problemas na Re-Alocacao da memoria");
  }

  fseek(fp_in2,0,SEEK_END);
  M = ftell(fp_in2)/2/(sizeof(double)+1);
  rewind(fp_in2);
  
  u = (double **)calloc(M,sizeof(double*));
  if(u==NULL) puts("Problemas na Alocacao da memoria");
  
  for(i=0;i<M;i++) {
    u[i] = (double *)calloc(2,sizeof(double));
    for(j=0;j<2;j++) {
      m=fscanf(fp_in2,"%lf",&(u[i][j]));
      if(m>0) mi+=m;
    }
  }

  mi/=2;
  printf("mi=%ld, M=%ld\n",mi,M);

  //printf(" --input-- \n");
  //print_matrix(v,N,2);

  if(mi<M) {
    M=mi;
    u = (double **)realloc(u,M*sizeof(double*));
    if(u==NULL) puts("Problemas na Re-Alocacao da memoria");
  }
  /*printf("\n --output-- \n");  
  print_matrix(v,N,2);
  printf("\n");
  */

  if(N==M)
    for(i=0;i<N;i++)
      if(v[i][0]==u[i][0])
	fprintf(fp_out,"%lf\t%lf\n",v[i][0],v[i][1]-u[i][1]);
      else {
	printf("x diferentes na linha %d\n",i);
	exit(6);
      }
  else {
    printf("O numero de pontos dos dois ficheiros e' diferente: N=%d, M=%d\n",N,M);
    exit(5);
  }
  
  fcloseall();
  
  for (i=0; i<N; i++)
    free(v[i]);
  free(v);

  for (i=0; i<M; i++)
    free(u[i]);
  free(u);
}
