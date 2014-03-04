#include <stdio.h> 
#include <math.h> 
#include <gsl/gsl_integration.h> 

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

struct quadratic_params{double a,b,c;};

double f (double t, void * params) { 
  struct quadratic_params *p 
    = (struct quadratic_params *) params;
  
  double l = p->a;
  double m = p->b;
  double q = p->c;
  
  double f = l*pow(t,1/2-m)*pow(1+pow(2*l*pow(t,1-m),q),-1./q);
  return f/(2*pow(t,3/2)); 
} 

int main (void) 
{ 
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000); 
 
  double result, error; 
  double x_lo = 5e-6, x_hi = 1e3;
  struct quadratic_params param 
    = {params[S][0], params[S][1], params[S][2]};
  gsl_function F; 
  
  F.function = &f; 
  F.params = &param;
  
  gsl_integration_qags (&F, x_lo, x_hi, 0, 1e-7, 1000, 
			w, &result, &error); 
  
  printf ("result = % .18f\n", result); 
  printf ("estimated error = % .18f\n", error); 
  printf ("intervals = %d\n", w->size); 
  
  gsl_integration_workspace_free (w); 
  
  return 0; 
} 
