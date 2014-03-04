#include <stdio.h>

// number of particles to shoot at each run (== #events)
#define N 1

#define MAX 150
#define MIN 5
#define STEP 5

// (don't) print very small recoil energies:
// 0 or 1 (boolean flag)
#define SMALL 1 

// material of the target (world):
// Silicon, Germanium, Argon, Xenon
#define MAT "Silicon"

// particle type:
// Xe: Z=54, A=131; Ge: Z=32, A=73; Ar: Z=18, A=40; Si: Z=14, A=28
#define Z 14
#define A 28

main()
{
  int i;
  char f[]="nuclear.mac",f1[]="recoils1.mac",f2[]="recoils2.mac";
  FILE *fp;

  if(SMALL)
    fp=fopen(f2,"w");
  else
    fp=fopen(f1,"w");

  // Set the verbosity of the output
  fprintf(fp,"/run/verbose 0\n/event/verbose 0\n/tracking/verbose 2\n");
  // Set material of the target (world)
  fprintf(fp,"/N02/det/setTargetMaterial %s\n",MAT);
  // Set particle type
  fprintf(fp,"/gun/particle ion\n/gun/ion %d %d %d\n",Z,A,Z);
  
  // Set kinetic energy and Start run defining the number of events
  for(i=MAX;i>=MIN;i-=STEP)
    fprintf(fp,"/gun/energy %d keV\n/run/beamOn %d\n",i,N);

  //print very small recoil energies
  if(SMALL)
    for(i=MIN-1;i>0;i--)
      fprintf(fp,"/gun/energy %d keV\n/run/beamOn %d\n",i,N);
  
  fclose(fp);
}
