#include <stdio.h>
#include <math.h>

#define MAX 3e11
#define MIN 1e7
#define N 100

#define M_Xe 122298642 // massa do Xe em keV/c^2
#define M_X 50e6 // massa do neutralino em keV/c^2

#define AMU2PDG 931494 // para converter de amu para keV

main()
{
  int i;
  double x, y, min, max;

  min=log10(MIN);
  max=log10(MAX);

  for(i=0;i<=N;i++) {
    x=pow(10,min+(max-min)/N*i);
    y=2e-2*x*x*pow(1/137.,2)*M_Xe/pow(M_Xe+x,2);

    //    printf("%g\t%g\n",x*pow(10,-6),y);
  }

  for(i=1;i<=226;i++) {
    x=i*AMU2PDG;
    y=2e-2*pow(M_X,2)*pow(1/137.,2)*x/pow(M_X+x,2);

    printf("%d\t%g\n",i,y);
  }
}
