#include <stdio.h>

#define MAX 150
#define MIN 5
#define STEP 5
#define N 1

#define Z 54
#define A 131

main()
{
  int i;
  FILE *fp;

  fp=fopen("recoils-gps.mac","w");

  // Set the verbosity of the output  
  fprintf(fp,"/run/verbose 0\n/event/verbose 0\n/tracking/verbose 2\n");
  // Set particle type
  fprintf(fp,"/gps/particle ion\n/gps/ion %d %d %d\n",Z,A,Z);
  // Set source geometry and starting position of the particle
  fprintf(fp,"/gps/pos/type Point\n/gps/pos/centre 0 0 0\n");
  // Set isotropic direction distribution
  fprintf(fp,"/gps/ang/type iso\n");

  // Set kinetic energy and Start run defining the number of events
  for(i=MAX;i>=MIN;i-=STEP)
    fprintf(fp,"/gps/energy %d keV\n/run/beamOn %d\n",i,N);

  fclose(fp);
}
