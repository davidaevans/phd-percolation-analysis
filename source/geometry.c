#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

/*
Returns the area of a 2D disc of given diameter and shell
thickness "shell" in units of D^2
*/

double disc_area(double diameter)
{
   return SQ(diameter/2) * M_PI;
}

/*..............................................................................*/

/*
Randomly assign species labels to the discs based on the fraction of each
specified by the user.  Returns the number of species 1.
*/

long assign_diameters(long npart, struct disc *particle, double fraction,
   double *diameter)
{
   long i;
   long n1;

   n1 = 0;
   for (i=0; i<npart; i++) {
      if (ran2(&seed) < fraction) {
         particle[i].species = 0;
         particle[i].diameter = diameter[0];
         n1++;
      } else {
         particle[i].species = 1;
         particle[i].diameter = diameter[1];
      }
   }

   return (n1);
}

/*..............................................................................*/

/*
Returns the cell number for a position in reduced coordinated
(on the unit square -0.5<=x<0.5, -0.5<=y<0.5).
*/

long getcell(struct vector pos, long ncellx, long ncelly)
{
   return (long)( (pos.y+0.5) * ncelly) * ncellx + (long)( (pos.x+0.5) * ncellx);
}

/*..............................................................................*/
