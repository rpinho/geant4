#include <stdio.h>
#include <math.h>

#define MinKineticEnergy 0.01 // energy in keV

#define MIN MinKineticEnergy
#define MAX 100 // energy in keV
#define lMIN -3
#define lMAX -1.15
#define N 100 // num. de pontos

// Projectile nucleus
#define z1 54
#define m1 130.876

// Target nucleus
#define z2 54
#define m2 131.303
#define TotNbOfAtomsPerVolume 13.30 // nm^-3

// Tmax = gamma*E; approx 1 (resultado e' igual)
#define gamma 4*m1*m2/(m1+m2)/(m1+m2)

// Ziegler's constants
#define a 1.1383
#define b 0.01321
#define c 0.21226
#define d 0.19593

// Screening length choice
#define screeningLength 'U' // L=Lindhard, F=Firsov, U=Ziegler

// Screening function choice
#define S 4

// Eckstein constants
double params[6][3] = {
  {1.309, 0.333, 0.667}, // S=0 Thomas-Fermi
  { 1.70, 0.311, 0.588}, // S=1 Thomas-Fermi-Sommerfeld
  { 2.37, 0.103, 0.570}, // S=2 Bohr
  { 2.92, 0.191, 0.512}, // S=3 Lenz-Jensen
  { 3.07, 0.216, 0.530}, // S=4 Moliere
  { 3.35, 0.233, 0.445}  // S=5 Kr-C
};

// Eckstein parameterisation of f(x)
double fEckstein(double x) {
  double l,m,q;
  l = params[S][0];
  m = params[S][1];
  q = params[S][2];
  return l*pow(x,1-2*m)*pow(1+pow(2*l*pow(x,2-2*m),q),-1./q);
};

// Eckstein parameterisation of f(t^1/2)
double ftEckstein(double t) {
  double l,m,q;
  l = params[S][0];
  m = params[S][1];
  q = params[S][2];
  return l*pow(t,1/2-m)*pow(1+pow(2*l*pow(t,1-m),q),-1./q);
};

// Ziegler
double fZiegler(double t, void * params) {

  // Ziegler's constants
  double a = 1.1383,b = 0.01321,c = 0.21226,d = 0.19593;
  double x,A,B,f,g;

  x = sqrt(t);
  A = 1+a*x;
  B = x+b*pow(x,c)+d*sqrt(x);
  f = log10(A)/2/B + a*x/2/A/B - x*log10(A)*(1+b*c*pow(x,c-1)+d/2/sqrt(x))/2/B/B;
  g = log10(A)/B/x/x + a/A/B/x - log10(A)*(1+b*c*pow(x,c-1)+d/2/sqrt(x))/B/B/x;

  return f/(2*pow(t,3./2)); // F(t^1/2)
}

void splineFitParametersInitialization();
void littmarkSpline();
void zieglerSpline();

main()
{
  double min=lMIN,max=lMAX,step,i,x;
    
  step=(double)(max-min)/N;
  for(i=min;i<max;i+=step)
    {
      x=pow(10,i);
      //      printf("%g\t%g\n",x,fZiegler(x));
      //      printf("%g\t%g\n",x,fEckstein(x));
      printf("%g\t%g\t%g\n",x,ftEckstein(x),ftEckstein(x)/(2*pow(x,3./2)));
    }

  //  puts("x\t\t\tm\t\t\tm-mOld\t\tl\t\tfspline");
  //  littmarkSpline();
  //  zieglerSpline();
}

void zieglerSpline()
{
  int j;
  double min=lMIN,max=lMAX,step,i,x;
  double sum,m,l,lOld,mOld=0,fspline;
  // Spline fit constants
  double C[] = {-2.432,-0.1509,2.648,-2.742,1.215,-0.1665};
  
  step=(double)(max-min)/N;
  for(i=min;i<max;i+=step)
    {
      x=pow(10,i);
      sum=0;
      for(j=0;j<6;j++)
	sum += C[j]*pow(0.1*(i-min),j);
      m = 1-exp(-exp(sum));
      if(i==min) l = 124.6;
      else l = lOld*pow(x,2*(m-mOld));
      fspline=l*pow(x,1-2*m);
      //      printf("%g\t%g\n",x,fspline);
      printf("%g\t\t%g\t\t%g\t\t%g\t\t%g\t\t%g\n",x,m,m-mOld,l,fspline,sum);
      mOld=m;lOld=l;
    }
}

void littmarkSpline()
{
  double min=lMIN,max=lMAX,step,i,x;
  double m,l,lOld,mOld=0,fspline;
  
  step=(double)(max-min)/N;
  for(i=min;i<max;i+=step)
    {
      x=pow(10,i);
      m = 1-exp(-0.2089-3.235e-9*pow(10+log10(x),8.58));
      if(i==min) l = 7.575;
      else l = lOld*pow(x,2*(m-mOld));
      fspline=l*pow(x,1-2*m);
      //      printf("%g\t%g\n",x,fspline);
      printf("%g\t\t%g\t\t%g\t\t%g\t\t%g\n",x,m,m-mOld,l,fspline);
      mOld=m;lOld=l;
    }
}
