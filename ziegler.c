#include <stdio.h>
#include <math.h>

#define MIN -9
#define MAX 6
#define N 100

#define a 1.1383
#define b 0.01321
#define c 0.21226
#define d 0.19593

main()
{
  double STEP,i,A,B,x,f,g;
  
  STEP=(double)(MAX-MIN)/N;

  //  printf("%lf\n",log(exp(1)));

  for(i=MIN;i<MAX;i+=STEP)
    {
      x=pow(10,i);
      A=1+a*x;
      B=x+b*pow(x,c)+d*sqrt(x);
      f=log10(A)/2/B + a*x/2/A/B - x*log10(A)*(1+b*c*pow(x,c-1)+d/2/sqrt(x))/2/B/B;
      //      printf("%g\t%g\n",x,f);

      g=log10(A)/B/x/x + a/A/B/x - log10(A)*(1+b*c*pow(x,c-1)+d/2/sqrt(x))/B/B/x;
      printf("%g\t%g\n",x,g);
    }
}
