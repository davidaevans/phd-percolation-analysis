#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

/*
Dumps a configuration to the supplied file handle.  The director is normalised
to the length of the rod.
*/

void draw(FILE *outfile, struct vector box, long npart,
          struct disc *particle)
{
   long i;

   for (i=0; i<npart; i++) {
      fprintf (outfile,
         "%15.8le %15.8le 0.0   0.0 0.0 0.0\n",
         box.x * (particle[i].pos.x - anint(particle[i].pos.x)),
         box.y * (particle[i].pos.y - anint(particle[i].pos.y))
         );
   }
}

/*..............................................................................*/


/*
Writes the configuration as above, but in a different format for debugging whether the
percolating cluster detection is working correctly

Format:

sweep# xpos ypos radius percolating

*/

void write_config(FILE *testfile, struct vector box, long npart,
                  struct disc *particle, long sweep, int percolating)
{
   long i;

   for (i=0; i<npart; i++) {
      fprintf (testfile,
         "%ld %15.8le %15.8le %lf %d %d\n",
         sweep, 
         box.x * (particle[i].pos.x - anint(particle[i].pos.x)),
         box.y * (particle[i].pos.y - anint(particle[i].pos.y)),
         particle[i].diameter,
         percolating,
         particle[i].structure
      );
   }
}

/*..............................................................................*/


/*
Writes un-normalised and normalised cluster size histograms for a given shell thickness

Format:

clustersize unnorm norm 

*/

void write_distribution(double shell, long *clustersize, long *percclusters,
                        double *normclustersize_ns,
                        double *normclustersize_sns, 
                        double *normclustersize_ns_expc,
                        double *normclustersize_sns_expc,
                        long npart,
                        char *fulldir) {
   
   long i;
   char filename[18];
   FILE *wfile;
   char fullname[100];
   sprintf(filename, "cs-dist-%.3lf.dat", shell);
   sprintf(fullname, "%s/%s", fulldir, filename);
   printf("Writing histogram to file: %s\n", filename);

   wfile = fopen(fullname, "w");
   if (wfile == NULL) die ("Could not open cluster distribution histogram file for writing.");
   fprintf(wfile, "#clustersize numclusters numclusters_ex_perc ns sns ns_experc sns_experc\n");
   for (i=0;i<npart;i++) {
      fprintf(wfile,
         "%ld %ld %ld %lf %lf %lf %lf\n",
         i + 1,
         clustersize[i],
         clustersize[i] - percclusters[i],
         normclustersize_ns[i],
         normclustersize_sns[i],
         normclustersize_ns_expc[i],
         normclustersize_sns_expc[i]);
   }
}

void write_coordinationnumbers(long *coordnum, double shell, long npart, char *fulldir) {
   
   long i;
   char filename[20];
   char fullfilename[100];
   FILE *wfile;

   sprintf(filename, "coord-num-%.3lf.dat", shell);
   sprintf(fullfilename, "%s/%s", fulldir, filename);
   printf("Writing coordination numbers to file: %s\n\n", fullfilename);
   
   wfile = fopen(fullfilename, "w");
   if (wfile == NULL) die ("Could not open coordination number file for writing.");
   fprintf(wfile, "#coordination_number number_of_particles\n");

   for (i=0;i<npart;i++){
      fprintf(wfile,
         "%ld %ld\n",
         i+1,
         coordnum[i]);
   }
}