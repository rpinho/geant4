#include <stdio.h> 
#include <math.h> 
#include <gsl/gsl_integration.h> 

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

// Screening length choice
#define screeningLength 'U' // L=Lindhard, F=Firsov, U=Ziegler

// Screening function choice
#define S 0

// Eckstein constants
double params[6][3] = {
  {1.309, 0.333, 0.667}, // S=0 Thomas-Fermi
  { 1.70, 0.311, 0.588}, // S=1 Thomas-Fermi-Sommerfeld
  { 2.37, 0.103, 0.570}, // S=2 Bohr
  { 2.92, 0.191, 0.512}, // S=3 Lenz-Jensen
  { 3.07, 0.216, 0.530}, // S=4 Moliere
  { 3.35, 0.233, 0.445}  // S=5 Kr-C
};

struct function_params{double a,b,c;};

double f (double t, void * params) { 
  struct function_params *p 
    = (struct function_params *) params;
  
  double l = p->a;
  double m = p->b;
  double q = p->c;
  
  double f = l*pow(t,1./2-m)*pow(1+pow(2*l*pow(t,1-m),q),-1./q);

  printf("%g\n",f/(2*pow(t,3./2)));
  return f/(2*pow(t,3./2)); // F(t^1/2)
}

// Ziegler
double fZiegler(double t, void * params) {

  // Ziegler's constants
  double a = 1.1383,b = 0.01321,c = 0.21226,d = 0.19593;
  double x,A,B,f;

  x = sqrt(t);
  A = 1+a*x;
  B = x+b*pow(x,c)+d*sqrt(x);
  f = log10(A)/2/B + a*x/2/A/B - x*log10(A)*(1+b*c*pow(x,c-1)+d/2/sqrt(x))/2/B/B;
  
  printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	 t,x,A,B,log10(A),log10(A)/2,log10(A)/2/B,x/2,a*x/2/A/B,d/2,
	 1+b*c*pow(x,c-1),x*log10(A)*(1+b*c*pow(x,c-1)+d/2/sqrt(x))/2,f,2*pow(t,3./2));
//  printf("%g\n",f/(2*pow(t,3./2)));

  return f/(2*pow(t,3./2)); // F(t^1/2)
}

int main (void) 
{ 
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000); 
 
  double result, error;
  double x_lo = 5e-6, x_hi = 1e3;
  double scale,a,energy,er,meanFreePath;
  double i,min=log10(MIN),max=log10(MAX),step=(double)(max-min)/N;
  struct function_params p 
    = {params[S][0], params[S][1], params[S][2]};
  gsl_function F; 

  // Normalization check
  /*{
    x_hi = 1e4;
    x_lo = 1e-3;
    
    F.function = &f;//Ziegler; 
    F.params = &p;
    
    gsl_integration_qags (&F, x_lo, x_hi, 0, 1e-7, 1000, 
			  w, &result, &error); 
    
    printf ("result = %g\n", result); 
    printf ("estimated error = %g\n", error); 
    printf ("intervals = %d\n", w->size);
  }*/
  
  switch(screeningLength) {
  case 'U': // Ziegler: Geant4's choice for ICRU_R49 and Ziegler1985
    scale = ( pow(z1, .23) + pow(z2, .23) ); break;
  case 'L': // Lindhard: Geant4's choice for Ziegler1977
    scale = sqrt( pow(z1, .667) + pow(z2, .667) ); break;
  case 'F': // Firsov: a popular choice
    scale = pow ( sqrt(z1) + sqrt(z2), .667); break;
  default : // Ziegler
    scale = ( pow(z1, .23) + pow(z2, .23) );
  }
  a = 0.8853 * 0.0529 / scale; // the screeningLength in nm

  for(i=min;i<=max;i+=step)
    {  
      energy = pow(10,i);
      // reduced energy 
      er = 32.536 * m2/(m1 + m2) / (z1 * z2 * scale) * energy ;
 
      x_hi = er*er;
      x_lo = x_hi*MinKineticEnergy/(gamma*energy);
      
      F.function = &f;//Ziegler; 
      F.params = &p;
      
      gsl_integration_qags (&F, x_lo, x_hi, 0, 1e-7, 1000, 
			    w, &result, &error); 
      /*
      printf ("result = %g\n", result); 
      printf ("estimated error = %g\n", error); 
      printf ("intervals = %d\n", w->size);
      */
     
      result *= M_PI*a*a; // scale back to physical units: nm^2
      //      printf("%g\t%g\n",energy,result);
      
      meanFreePath = 1. / (TotNbOfAtomsPerVolume*result); // nm
      //      printf("%g\t%g\n",energy,meanFreePath);
    }
  
  gsl_integration_workspace_free (w); 
  
  return 0; 
} 
