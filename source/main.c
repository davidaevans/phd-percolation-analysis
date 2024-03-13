#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

int main()
{
   char filename[80];          /*Name of file where LAMMPS configs are */
   char blockdir[60];          /* Name of the directory for block sizes */
   char fulldir[80];           /* Name of the directory for block and density */
   double fraction;            /* Fraction of rod species 1 */
   double *diameter;             /* Length of the cylinders */
   //double meanl;               /* Mean L/D */
   //double moment, moment0;     /* Moments of the L/D distribution */
   double smin;                /* minimum shell diameter */
   double smax;                /* maximum shell diameter */
   double sstep;               /* step in shell diameter */
   double dr;                  /* step in r for binning RDF */
   int fmethod;                /* Method for counting fraction of rod species */
   int periodic;               /* 1=wrapping criterion, 0=spanning */
   int rdfflag;                /* boolean flag for whether to calculate RDF or not */
   long sourcetype;             /* 0 = LAMMPS, 1 = MC sim */
   long typetest;               /* 0 = test all structure types */
   long i;
   long movie;                 /* Number of sweeps between movie frames */
   long configs;               /* Number of sweeps between writing configs for debugging */
   long npart;                 /* Number of particles */
   long nsweeps;               /* Number of production sweeps */
   long nblocks;               /* Number of blocks for averaging over */
   long blocksize;             /* Number of sweeps per block */
   long blockstart, blockend;  /* Sweep to start and end blocks */
   long report;                /* Number of sweeps between statistics reports */
   long equilibrate;           /* Number of sweeps to wait before sampling statistics */
   struct disc *particle;  /* Configuration of entire system */
   struct vector box;          /* Simulation cell dimensions */
   //struct mystat st;

   printf ("\nDiscs in 2D");
   printf ("\n--------------------------\n\n");
   

   diameter = (double *)malloc( 2 * sizeof(double) );

   /* Get user parameters */
   read_options(&npart, &box, diameter, &nsweeps, &report, &movie,
                &fraction, &fmethod, &periodic, &configs, 
                &typetest, &smin, &smax, &sstep, &equilibrate, filename, &dr, 
                &rdfflag, &sourcetype, &nblocks);
   
   /* Set aside memory for the configuration */
   particle = (struct disc *)malloc(npart * sizeof(struct disc));
   blocksize = (long) (nsweeps-equilibrate)/(nblocks);
   for (i=0; i<npart; i++) { particle[i].idx=i; }

   /* Display moments of the length distribution 
   printf ("Moments about zero and about mean of target L/D distribution:\n");
   meanl = fraction * diameter[0] + (1.0 - fraction) * diameter[1];
   for (i=0; i<=6; i++) {
      moment0 = pow(diameter[0], (double)i) * fraction + pow(diameter[1], (double)i) * (1.0 - fraction);
      moment = pow(diameter[0] - meanl, (double)i) * fraction
             + pow(diameter[1] - meanl, (double)i) * (1.0 - fraction);
      printf ("moment %ld:  %.8le  %.8le\n", i, moment0, moment);
   }
   printf ("\n");
   fflush (stdout);
   */
   // Calculate properties per block
   printf ("Starting simulation\n\n");
   fflush (stdout);
   for (i=0; i<nblocks; i++){
      printf("Beginning block %ld\n", i);
      blockstart = equilibrate + i*blocksize;
      blockend = blockstart + blocksize;
      printf("Start sweep: %ld\nEnd sweep: %ld\n", blockstart, blockend);
      sprintf(blockdir, "block%ld-%.0ld-%.0ld", i, blockstart, blockend);
      // check if dir has already been created and create it if not
      DIR* dir = opendir(blockdir);
      if (dir) {
         /* Directory exists. */
         closedir(dir);
      } else if (ENOENT == errno) {
         mkdir(blockdir,0777);
      } else {
         die ("Could not determine whether directory exists or not");
      }
      sprintf(fulldir, "%s/rho%.3lf", blockdir, (npart/(box.x*box.y)));
      DIR* dendir = opendir(fulldir);
      if (dendir) {
         /* Directory exists. */
         closedir(dendir);
      } else if (ENOENT == errno) {
         mkdir(fulldir,0777);
      } else {
         die ("Could not determine whether density directory exists or not");
      }

      simulate(npart, fraction, &box, diameter, blockend, report, movie,
         periodic, particle, configs, typetest, smin, smax, sstep, blockstart,
         filename, dr, rdfflag, sourcetype, fulldir);


   }


   
   
   // simulate(npart, fraction, &box, diameter, nsweeps, report, movie,
   //    periodic, particle, configs, typetest, smin, smax, sstep, equilibrate,
   //    filename, dr, rdfflag, sourcetype);
   
   printf ("\nDone\n\n");

   return 0;
}

/*..............................................................................*/

/*
Print error message and exit.
*/

void die(char string[])
{
   fprintf (stderr, "\nERROR: %s\n\n", string);
   exit (1);
}

/*................................................................................*/
