#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define M 500    /* numero de pontos de output */
#define MAX 3    /* limite maximo da energia em log10 (o max pode ser menor) */
#define MIN 3    /* energia minima apenas como valor de referencia inicial em log10 (o min devera ser menor) */

/* massas atomicas em unidades de g/mole ou amu */
#define M_p 1.01 //1.007276
#define M_a 4.00 //4.0026
#define M_N 14.01
#define M_Al 26.98
#define M_Si 28.09
#define M_Ar 39.95
#define M_Cu 63.55
#define M_Ge 72.64
#define M_Xe 131.29 // ou .30
#define M_Pb 207.19

/* densidades em unidades de g/cm3 */
#define d_Al 2.70
#define d_Cu 8.96
#define d_Pb 11.35
#define d_Si 2.33
#define d_LAr 1.39
#define d_Ge 5.32
#define d_LXe 2.9
#define d_XeGas 0.005458

/* particula incidente */
#define AI 0.0  /* massa atomica da particula incidente */

/* normalizacao do output */
#define AN 0.0     /* massa atomica do iao de normalizacao do output */
#define Z 0    /* numero atomico do iao de normalizacao do output: protao=1, alfa=2, Al=13 */

/* alvo */
#define AA 0.0  /* massa atomica do material do alvo */
#define RO 0.0  /* densidade do material do alvo */

/* constantes */
#define NA 6.022e23 /* numero de Avogrado em unidades /mole */

void joao_viana(FILE *fp,double **v,long int N,double *min,double *max);

void print_matrix(double **matrix, int m, int n) {
  int i,j;
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++)
      printf("%lf ",matrix[i][j]);
    printf("\n");
  }
}

main(int argc, char*argv[])
{
  char f_out[50],*name="_logbin_average.dat";
  FILE *fp_in,*fp_out;
  int n,j=0,l;
  long int N,i,ni=0;
  double **v,min,max;

  if(argc!=2) {
    printf("Sintaxe:\n\n%s ficheiro\n\n",argv[0]);
    exit(1);
  }

  if((fp_in=fopen(argv[1],"r"))==NULL) {
    printf("Impossivel abrir o ficheiro %s\n",argv[1]);
    exit(2);
  }

  strcpy(f_out,argv[1]);
  l=strlen(f_out)-4;
  while(f_out[l]=name[j]) {
    l++;j++;
  }

  if((fp_out=fopen(f_out,"w"))==NULL) {
    printf("Impossivel criar o ficheiro %s\n",f_out);
    exit(3);
  }

  puts("\nAtencao a massa da particula incidente, massa e carga do iao de normalizacao e massa e densidade do alvo definidas no codigo ");
  printf("Mi = %.2lf amu\t\tMn = %.2lf amu\t\tZ = %d\t\tMa = %.2lf amu\t\tdensidade = %.2lf g/cm3\n\n",AI,AN,Z,AA,RO);

  fseek(fp_in,0,SEEK_END);
  N = ftell(fp_in)/2/(sizeof(double)+1);
  rewind(fp_in);
  
  v = (double **)calloc(N,sizeof(double*));
  if(v==NULL) puts("Problemas na Alocacao da memoria");
  
  for(i=0;i<N;i++) {
    v[i] = (double *)calloc(2,sizeof(double));
    for(j=0;j<2;j++) {
      n=fscanf(fp_in,"%lf",&(v[i][j]));
      if(n>0) ni+=n;
    }
    v[i][0]=log10(v[i][0]);
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

  /*printf("\n --output-- \n");  
  print_matrix(v,N,2);
  printf("\n");
  */

  joao_viana(fp_out,v,N,&min,&max);
  
  printf("\nA energia (x) foi dividida por %lu bins logaritmicamente espacados de %g a %g\n\nImpresso o ficheiro %s\n\n",M,pow(10,min),pow(10,max),f_out);
  
  fcloseall();
  for (i=0; i<N; i++)
    free(v[i]);
  free(v);
}


void joao_viana(FILE *fp,double **v,long int N,double *min,double *max)
{
  long int i,k,n[M]={0};
  double x,y,r,energy,media[M]={0};

  *min=MIN;*max=0; /* valores iniciais apenas */

  for(i=0;i<N;i++) {
    if(v[i][0] && v[i][0]<*min)
      *min=v[i][0];//printf("min = %lf\n",min);}
    if(v[i][0]<MAX && v[i][0]>*max)
      *max=v[i][0];//printf("max = %lf\n",max);}
  }
   
  for(i=0;i<N;i++)
    if(v[i][0]<MAX && v[i][0]*v[i][1]) {
      r=(v[i][0]-*min)/(*max-*min)*(M-1);//printf("r = %lf\t",r);
      k=(long int)rint(r);//printf("k = %lu\n",k);
      media[k]+=v[i][1];
      n[k]++;//printf("n[%lu] = %lu\n",k,n[k]);
    }

  for(k=0;k<M;k++)
    if(media[k]){
      energy=pow(10,*min+(double)k*(*max-*min)/(M-1));
      if(energy) {

	/***************** unidades de output ****************************/
	/* a energia cinetica T de imput, energy, vem em unidades de MeV */
	x=energy;        /* protao: MeV ou ja vem nas unidades que quero */
	//	x*=1000;                                         /* Alessio: keV */
	//	x/=AI;                                          /* alfa: KeV/amu */
	//	x*=AN;                      /*  outros ioes: KeV/amu normalizado */
	/* o dE/dx de imput, media[], vem em unidades de MeV/mm **********/
	y=media[k]/n[k];      /* MeV/mm ou ja vem nas unidades que quero */
	//	y*=10;                                                 /* MeV/cm */
	//	y/=RO;                                   /* protao: MeV cm^2 / g */
	//	y/=1000;                               /* Alessio: MeV cm^2 / mg */
	//	y*=AA/NA*pow(10,24);                    /* alfa: eV cm^2 / 10^15 */
	//	y*=pow(Z,2);         /* outros ioes: eV cm^2 / 10^15 normalizado */
	/*****************************************************************/
	
	fprintf(fp,"%lf\t%lf\n",x,y);
      }
    }
}
