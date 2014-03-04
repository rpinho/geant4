#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* constantes */
#define NA 6.022e23           /* numero de Avogrado em unidades /mole */
#define factorAMU2PDG 931.495 /* Factor to convert AMU (g/mole) into PDG mass unit (MeV/c^2) */

/* numeros atomicos em unidades de e (carga elementar) */
#define Z_p 1
#define Z_a 2
#define Z_N 7
#define Z_Al 13
#define Z_Si 14
#define Z_Ar 18
#define Z_Ge 32
#define Z_Xe 54

/* massas atomicas em unidades de g/mole ou amu */
#define M_p 1.01 //1.007276
#define M_a 4.00 //4.0026
#define M_N 14.01
#define M_Al 26.98
#define M_Si 28.09
#define M_Ar 39.95
#define M_Ge 72.64
#define M_Xe 131.29 // ou .30
#define M_Pb 207.19

/* densidades em unidades de g/cm3 */
#define d_Al 2.70
#define d_Si 2.33
#define d_LAr 1.39
#define d_Ge 5.32
#define d_LXe 2.9
#define d_XeGas 0.005458
#define d_Pb 11.35

/* particula incidente */
#define Z1 Z_Xe  /* numero atomico da particula incidente */
#define M1 M_Xe  /* massa atomica da particula incidente */

/* alvo */
#define Z2 Z_Xe  /* numero atomico do material do alvo */
#define M2 M_Xe  /* massa atomica do material do alvo */
#define RO d_LXe /* densidade do material do alvo */

/* normalizacao do output */
#define MN 0.0   /* massa atomica do iao de normalizacao do output */
#define ZN 0     /* numero atomico do iao de normalizacao do output: protao=1, alfa=2, Al=13 */

void output(FILE *fp,double **v,long int N);

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
  char f_out[50],*name="_converted_units.dat";
  FILE *fp_in,*fp_out;
  int n,j=0,l;
  long int N,i,ni=0;
  double **v;

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

  puts("\nAtencao a massa e carga da particula incidente, massa, carga e densidade do alvo e massa e carga do iao de normalizacao definidas no codigo");
  printf("Z1 = %d\t\tM1 = %.2lf amu\t\tZ2 = %d\t\tM2 = %.2lf amu\t\tdensidade = %.2lf g/cm3\t\tMn = %.2lf amu\t\tZn = %d\n\n",Z1,M1,Z2,M2,RO,MN,ZN);

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

  output(fp_out,v,N);
  
  printf("\nImpresso o ficheiro %s\n\n",f_out);
  
  fclose(fp_in);fclose(fp_out);
  for (i=0; i<N; i++)
    free(v[i]);
  free(v);
}


void output(FILE *fp,double **v,long int N)
{
  int A1, A2;
  long int i;
  double x,y,rm,n,theZieglerFactorinverse,Linhard;

  A1 = (int) rint(M1);
  A2 = (int) rint(M2); printf("\nA1=%d\tA2=%d\n",A1,A2);
  /* Geant4's reduced mass */
  rm = (A1 + A2) * ( pow(Z1, .23) + pow(Z2, .23) ) ; printf("Geant4's reduced mass = %g\n",rm);
  /* Total Number Of Atoms Per Volume of material */
  n = RO*NA/M2; printf("Total Number Of Atoms Per Volume of material = %g nm^-3\n",n/pow(10,21));
  /* Factor to convert the Geant4 dE/dx [MeV/mm] unit into the Stopping Power unit [ev/(10^15 atoms/cm^2] */
  theZieglerFactorinverse = (pow(10,22)*M2)/(RO*NA); printf("theZieglerFactorinverse=%lf\n",theZieglerFactorinverse);
  /* Factor to return to Linhard's units */
  Linhard = 8.462 * Z1 * Z2 * A1 / rm ; printf("Linhard = %lf\n",Linhard);
  
  for(i=0;i<N;i++) {
    
    /***************** unidades de output ****************************/
    /* a energia cinetica T de imput, energy, vem em unidades de MeV */
    x=v[i][0];                                        /* protao: MeV */
    x*=1000;                                         /* Alessio: keV */
    //    x/=M1;                                          /* alfa: KeV/amu */
    //    x*=MN;                      /*  outros ioes: KeV/amu normalizado */
    x=sqrt(32.536*A2/(Z1*Z2*rm)*x); /* Linhard's reduced energy SQRT */
    /* o dE/dx de imput vem, por defeito, em unidades de MeV/mm ******/
    y=v[i][1];                                             /* MeV/mm */
    //    y*=10;                                                 /* MeV/cm */
    //    y/=RO;                                   /* protao: MeV cm^2 / g */
    //    y/=1000;                               /* Alessio: MeV cm^2 / mg */
    //    y*=M2/NA*pow(10,24);   /* alfa: eV cm^2 / 10^15, Ziegler's units */
    //    y*=pow(ZN,2);         /* outros ioes: eV cm^2 / 10^15 normalizado */
    //   y*=pow(10,3);  /* se o dE/dx de imput vier em unidades de MeV/um */
    y*=theZieglerFactorinverse/Linhard;           /* Linhard's units */
    /*****************************************************************/
    
    fprintf(fp,"%lf\t%lf\n",x,y);
  }
}
