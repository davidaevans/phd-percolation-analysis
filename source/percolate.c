#include "global.h"
#include "prototypes.h"

/*................................................................................*/

/*
Tests for percolation in the x and y directions seperately, returning 1 or 0
(as double precision values) in each component of a vector for percolation and no
percolation, respectively.  Also returns the total number of contacts in the
system.
*/

struct vector percolate(long npart, struct disc *particle,
   struct vector box, long *ctot,
   struct mystat *npc1, struct mystat *npc2,
   struct disc **cfirst, long **neighbour, long typetest, double shell,
   long *clustersize, long *percclusters, long *coordnum, char *fulldir,
   struct mystat *avecoordnum)
{

   //int analyse;
   long build;
   long explore;
   long n[2];
   long placed;
   long first;
   long i, j, k;
   long ind;
   struct vector link;
   struct vector psep;
   //struct vector pimg;
   struct vector pcomp;
   static int *done;
   static long **conn=NULL; // array of particles - [particle number][connection number]
   static long *members=NULL; /* list of particles in the cluster */
   static long *nc=NULL;
   static struct vector *phys=NULL;
   FILE *clusterfile = NULL;
   char filename[22];
   char fullfilename[100];

   sprintf(filename, "cluster_ids_%.3lf.dat", shell);
   sprintf(fullfilename, "%s/%s", fulldir, filename);
   clusterfile = fopen(fullfilename, "a");
   /* Allocate memory on first call, then keep it */
   if (!nc) {
      conn = (long **)malloc(sizeof(long *) * npart); //2D array - npart x maxconnections
      for (i=0; i<npart; i++) conn[i]=(long *)malloc(sizeof(long) * MAXO);
      done = (int *)malloc(sizeof(int) * npart);
      members = (long *)malloc(sizeof(long) * npart); /* list of particles in the cluster */
      nc = (long *)malloc(sizeof(long) * npart);
      phys = (struct vector *)malloc(sizeof(struct vector) * npart);
   }

   /* Nonpercolating until proven percolating */
   pcomp.x = pcomp.y = 0.0;
   
   /* Compile lists of contacts for each spherocylinder */
   touching(npart, particle, box, nc, conn, ctot, cfirst, neighbour, shell);

   /* Test the cluster to which each particle belongs, unless it has already
      been tested. */

   for (i=0; i<npart; i++) {
      //reset counter
      done[i]=0;
      //add nc data to coord number array
      coordnum[nc[i]]++;
      //add nc data to average coordination number
      accumulate(avecoordnum, nc[i]);
   }

   for (first=0; first<npart; first++) {
         if (done[first]) continue;

         /* Construct a unified instance of the cluster */
         /* Put the first particle at the origin */
         members[0] = first;  //assign this particle to be the start of the cluster
         phys[first].x = phys[first].y = 0.0; // set the location of this particle to be origin
         done[first] = 1; // set this particle as done/visited
         placed = 1; // have set the particle at the origin
         n[0] = n[1] = 0; // set number of each species to 0
         n[particle[first].species] = 1; // add 1 to whichever species this particle belongs to
         explore = 0;
         //analyse = 0;
         /* Place touching particles relative to those already placed */
         // continue only if we have explored fewer than we've placed
         // i.e. when we have explored more particles than have been placed in the cluster, stop
         while (explore < placed) {
            ind = members[explore]; // index is the particle to explore around
            for (j=0; j<nc[ind]; j++) { //for each particle connected to this particle (nc[ind] is number of connections for particle ind)
               build = conn[ind][j]; //build index is connection j for particle ind
               if (!done[build]) { //if this particle hasn't been placed
                  members[placed] = build; //add this particle to list for cluster
                  done[build] = 1; //build particle is then done 
                  link = image(particle[build].pos, particle[ind].pos, box); // get vector between start of cluster particle 
                  /* vector of the particle that has just been placed/added to cluster, vector of particle 
                  we're exploring around + vector between that particle and the one being added
                  */
                  phys[build].x = phys[ind].x + link.x; 
                  phys[build].y = phys[ind].y + link.y;
                  placed++;
                  n[particle[build].species]++;
               }
               //go to next particle that is connected to particle being explore around
            }
            // go to next particle to explore around if we have placed another particle into the cluster
            explore++;
         }

         clustersize[placed-1]++; //add one to histogram of cluster size placed, remembering index starts from 0, so index 0 is cluster size 1


         /* Now check distances between particles that are bonded when the periodic
            boundary conditions are switched off */

         /* x percolation (turn off just x boundary conditions) */
         
         if (pcomp.x == 0.0) {   //Skip test if a previous cluster percolated 
            for (j=0; j<placed; j++) { // loop over all particles placed so far
               ind = members[j];
               for (i=0; i<nc[ind]; i++) { // loop over all connections of current particle
                  
                  k = conn[ind][i]; // k is index of particle connected to current particle 
                  
                  //Separation vector for current placed particle and particles connected to it
                  //In y direction, can be closest image, in x direction it is *not* closest image
                  
                  psep.x = phys[ind].x - phys[k].x; 
                  psep.y = phys[ind].y - phys[k].y;
                  psep.y = psep.y - box.y * anint(psep.y/box.y);
                   
                  //If this separation vector means the two particles are not overlapped
                  //then the cluster percolates in the x direcion
                  
                
                  if (!overlap(psep, particle[ind].diameter, particle[k].diameter, shell)) {
                     pcomp.x = 1.0;
                     goto ytest;   // Skip to test for y-percolation 
                  }
               }
            }
         }
         
         /* y percolation (turn off just y boundary conditions)
         Same as above but in y direction
          */
         ytest:;
         
         if (pcomp.y == 0.0) {  // Skip test if a previous cluster percolated 
            for (j=0; j<placed; j++) {
               ind = members[j];
               for (i=0; i<nc[ind]; i++) {
                  k = conn[ind][i];
                  psep.x = phys[ind].x - phys[k].x;
                  psep.y = phys[ind].y - phys[k].y;
                  psep.x = psep.x - box.x * anint(psep.x/box.x);

                  //separation vector for closest image of the neighbouring particles
                  
                  if (!overlap(psep, particle[ind].diameter, particle[k].diameter, shell)) {
                     pcomp.y = 1.0;
                     goto stoptest;
                  }
               }
            }
         }
         
         stoptest:;

    
         //if there is a percolating cluster, then add to percclusters array
         if (pcomp.x > 0.5 || pcomp.y > 0.5) {
            percclusters[placed-1]++;
         }

         //write all members of this cluster to the same line on a file
         for (i=0;i<placed;i++){
            if((members[i]+1)>npart) die("Particle ID greater than npart in percolate.c (cluster writing)");
	    fprintf(clusterfile,
                  "%ld ", members[i]+1); // the plus one for the id is to match with the LAMMPS input files which start at 1 instead of 0
         }
         fprintf(clusterfile, "\n");

         /*
         Analyses nematic order parameters for percolating clusters
         Doesn't make sense for 2D discs with no direction

         analyse is only 1 if a percolating cluster has been detected (see above code)

         Mark's comment: If a percolating cluster has been found then gather some statistics on it
         */

        
         /*
         if (analyse == 1) {
            
            accumulate(npc1, (double)n[0]);
            accumulate(npc2, (double)n[1]);
            //if (pcomp.x + pcomp.y > 1.99)  
            // if percolating cluster found in each direction, can exit early
               
            
         }
         */
         

         

   }  /* End of loop over particles as potential starting points for a chain */

   //escape:;
   fprintf(clusterfile, "\n");
   fclose(clusterfile);
   return pcomp;
}

/*................................................................................*/

/*
For each particle, makes a list of other spherocylinders connected to it
(i.e., with axes within a distance of "shell").  The list is returned in "conn"
and the number of connections for each spherocylinder is returned in "nc".
The total number of connections is returned as ctot.
*/

void touching(long npart, struct disc *particle, struct vector box,
              long *nc, long **conn, long *ctot,
              struct disc **cfirst, long **neighbour, double shell)
{
   long i;
   long *cell;
   struct disc *test;
   struct vector r_cm;

   *ctot = 0;
   for (i=0; i<npart; i++) nc[i]=0;

   /* Loop over all particles */
   for (i=0; i<npart; i++) {
      /* Loop over all cells adjacent to particle */
      cell = &neighbour[particle[i].cell][0];
      while (*cell >= 0) {
         /* Loop over all particles in cell */
         test = cfirst[*cell];
         while (test) {
            
            if (i != test->idx) {
               r_cm = image(particle[i].pos, test->pos, box);
               if ( overlap(r_cm, particle[i].diameter, test->diameter, shell) ) {
                  
                  if (nc[i] >= MAXO || nc[test->idx] >= MAXO) {
                     fprintf (stderr,
                        "ERROR: Compiled maximum number of overlaps per rod exceeded\n");
                     fprintf (stderr, "MAXO = %d\n", MAXO);
                     exit (99);
                  } 
                  conn[i][nc[i]] = test->idx; // Add test particle index to connection array for particle i at connection number
                  nc[i]++;
                  (*ctot)++; //total number of connections
               }  /* End of pair overlap test */
            }

         test = test->next; // go to next particle in this cell
         }  /* End of loop over particles in adjacent cell */

         cell++;
      }  /* End of loop of adjacent cells */
   }  /* End of loop over all particles */

   /* Correct for fact that each pair is treated twice. */
   *ctot /= 2;

}

/*................................................................................*/

/*
For each particle, makes a list of other particles connected to it.
The list is returned in "conn" and the number of connections 
for each particle is returned in "nc".
The total number of connections is returned as ctot.
*/

void touching1(long npart, struct disc *particle, struct vector box,
              long *nc, long **conn, long *ctot,
              struct disc **cfirst, long **neighbour, long typetest,
              double shell)
{
   long i,j;
   struct disc *test;
   struct vector r_cm;

   *ctot = 0;
   typetest = 0;
   for (i=0; i<npart; i++) nc[i]=0;

   /* Loop over all particles */
   for (i=0; i<npart; i++) {


      /* Loop over all Voronoi neighbours of particle i up to maximum of 10*/
      for (j=0; j<10; j++) {
         //break loop if neighbour doesn't exist
         if (particle[i].neighbours[j] == -1) break;

         long test_idx = particle[i].neighbours[j];
         test = &particle[test_idx];
         r_cm = image(particle[i].pos, test->pos, box);
         if (overlap(r_cm, particle[i].diameter, test->diameter, shell)){
            conn[i][nc[i]] = test->idx;
            nc[i]++;
            (*ctot)++;
         }
         /*
         //no structure or distance tests
         if (typetest == 0) {
            
            conn[i][nc[i]] = test->idx;
            nc[i]++;
            (*ctot)++; //increase total number of connections

         } else {
            
            //test for structure type
            if (test->structure == typetest && test->structure == particle[i].structure) {
               conn[i][nc[i]] = test->idx;
               nc[i]++;
               (*ctot)++; //increase total number of connections

            }
         }
         */


      } /* End loop over neighbours */


   }  /* End of loop over all particles */

   /* Correct for fact that each pair is treated twice. */
   *ctot /= 2;

}
