#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

void simulate(long npart, double fraction, struct vector *box, double *diameter,
              long nsweeps, long report, long movie, int periodic,
              struct disc *particle, long configs, long typetest,
              double smin, double smax, double sstep, long equilibrate,
              char *filename, double dr, int rdfflag, long sourcetype,
              char *fulldir)
{
   double longer;      /* Length of the longer of the two rod species */
   double pf;          /* Instantaneous packing fraction */
   double shell;       /* diameter ratio of the connectivity shell */
   double maxsep;      /* Maximum separation for binning of RDF */
   double shorter;     /* shorter box length to determine max cutoff for RDF */
   double meanclustersize;
   double meanclustersize_expc; /* Mean cluster size excluding any percolating clusters */
   long i, j;
   long ncellx, ncelly;      /* Number of cell-list cells in each direction */
   long ncells;        /* Total number of cell-list cells */
   long **neighbour;   /* List of neighbouring cells for each cell */
   long n1;            /* Number of species 1 in present configuration */
   //long next_dump;     /* Next sweep number for reporting statistics */
   //long next_frame;    /* Next sweep number for dumping a movie frame */
   //long next_debug;     /* Next sweep number for dumping configs for debugging */
   long nolap;         /* Current number of pairwise overlaps */
   long sweep;         /* Current sweep number */
   long numbins;        /* Number of bins in RDF histogram */
   long *rdfhist;      /* Array for RDF histogram */
   long *clustersize;   /* Array for histogram of cluster sizes */
   long *percclusters;  /* Array for storing sizes of percolating clusters */
   long *coordnum;      /* Array for coordination numbers of each particle */
   struct disc **cfirst;    /* Array of pointers to first particle in cell list */
   struct mystat contacts;  /* Shell contacts per particle statistics */
   struct mystat avecoordnum; /* Stats struct for average coordination number */
   struct mystat frac1;  /* Number fraction of species 1 */
   struct mystat npc1, npc2;   /* Number of each species in percolating clusters */
   struct mystat packfrac;     /* Packing fraction statistics */
   struct mystat percx;  /* Percolation statistics specifically in the x direction */
   struct mystat percy;  /* Percolation statistics specifically in the y direction */
   struct mystat percxyboth; /* Percolation statistics in the xy plane (x or y or both) */
   struct mystat percxandy; /* Percolation statistics for x and y directions */
   struct vector pcomp; /* Flags (1.0 or 0.0) for percolation in each of the Cartesian directions */
   FILE *mf = NULL;           /* Handle for movie file */
   FILE *testfile = NULL; /* debugging file */
   FILE *shellseries = NULL; /* file where time series of critical shell variables are written to */
   FILE *timeseries = NULL;  /* time series of whether percolating cluster is present in each frame */
   FILE *configfile = NULL;
   char shellseriesname[80];
   char timeseriesname[80];

   const struct mystat nullstat = {0.0, 0.0, 0, 0.0, 0.0};
   
   /*=== Initialise counters etc. ===*/
   //next_dump = report + equilibrate;
   //next_frame = movie + equilibrate;
   //next_debug = configs + equilibrate;

   packfrac = frac1 = contacts = avecoordnum = nullstat;
   npc1 = npc2 = nullstat;
   percx = percy = percxyboth = percxandy  = nullstat;

   numbins = 0;
   maxsep = 0;
   rdfhist = NULL;
   if (rdfflag) {
      shorter = box->y > box->x ? box->x : box->y;
      maxsep = shorter/2;
      numbins = maxsep/dr + 1;
      rdfhist = (long *)malloc(numbins * sizeof(long));
      //set all elements to zero as no guarantee they will all be accessed
      for (i=0;i<numbins;i++) rdfhist[i] = 0;
   }

   if (movie > 0) {
      mf = fopen("movie", "w"); 
   } else {
      mf = NULL;
   }

   if (configs > 0){
      testfile = fopen("test-data", "a");
   } else {
      testfile = NULL;
   }
   sprintf(shellseriesname, "%s/%s", fulldir, "shellseries.dat");

   shellseries = fopen(shellseriesname, "w");
   if (shellseries == NULL) die("Could not open shell series file");
   fprintf(shellseries, "#shell x xrms y yrms xyboth xybothrms xory xoryrms xandy xandyrms meanclustersize meanclustersize_exc_perc_clusters mean_coord_num\n");

   shell = smin;

   sprintf(timeseriesname, "%s/%s", fulldir, "timeseries.dat");

   timeseries = fopen(timeseriesname, "w");
   if (timeseries == NULL) die ("Could not open time series file");
   fprintf(timeseries, "#xyorboth\n");

   //increase shell thickness until reach maximum
   while (shell <= smax) {
      //This does not generalise if two diameters are not the same
      printf("Shell: %lf\n", shell);

      /*=== Initialise cell list ===*/
      longer = diameter[0];
      if (diameter[1] > diameter[0]) { longer = diameter[1]; }
      ncellx = (long)(box->x / (shell*longer)); 
      ncelly = (long)(box->y / (shell*longer));
      ncells = ncellx * ncelly;
      cfirst = (struct disc **)malloc(sizeof(struct disc *) * ncells);
      neighbour = (long **)malloc(sizeof(long *) * ncells);
      for (i=0; i<ncells; i++) {
         neighbour[i] = (long *)malloc(sizeof(long) * 10);
      }
      /* Work out neighbouring cells for each cell by pointer */
      /* Interior of box */
      for (i=0; i<ncellx; i++) {
         for (j=0; j<ncelly; j++) {
            neighbour[j*ncellx+i][0] = ((j-1) + (j==0?ncelly:0))*ncellx + ((i-1) + (i==0?ncellx:0));
            neighbour[j*ncellx+i][1] = ((j-1) + (j==0?ncelly:0))*ncellx + i;
            neighbour[j*ncellx+i][2] = ((j-1) + (j==0?ncelly:0))*ncellx + ((i+1) - (i==ncellx-1?ncellx:0));
            neighbour[j*ncellx+i][3] = j*ncellx + ((i-1) + (i==0?ncellx:0));
            neighbour[j*ncellx+i][4] = j*ncellx + i;
            neighbour[j*ncellx+i][5] = j*ncellx + ((i+1) - (i==ncellx-1?ncellx:0));
            neighbour[j*ncellx+i][6] = ((j+1) - (j==ncelly-1?ncelly:0))*ncellx + ((i-1) + (i==0?ncellx:0));
            neighbour[j*ncellx+i][7] = ((j+1) - (j==ncelly-1?ncelly:0))*ncellx + i;
            neighbour[j*ncellx+i][8] = ((j+1) - (j==ncelly-1?ncelly:0))*ncellx + ((i+1) - (i==ncellx-1?ncellx:0));
            neighbour[j*ncellx+i][9] = -1;  /* end token */
         }
      }
      if (!periodic) {
         /* Overwrite periodic results along the boundaries */
         /* Edges */
         for (i=1; i<ncellx-1; i++) {
            /* top */
            neighbour[i][0] = i-1;
            neighbour[i][1] = i;
            neighbour[i][2] = i+1;
            neighbour[i][3] = ncellx + (i-1);
            neighbour[i][4] = ncellx + i;
            neighbour[i][5] = ncellx + (i+1);
            neighbour[i][6] = -1;
            /* bottom */
            neighbour[(ncelly-1)*ncellx + i][0] = (ncelly-2)*ncellx + (i-1);
            neighbour[(ncelly-1)*ncellx + i][1] = (ncelly-2)*ncellx + i;
            neighbour[(ncelly-1)*ncellx + i][2] = (ncelly-2)*ncellx + (i+1);
            neighbour[(ncelly-1)*ncellx + i][3] = (ncelly-1)*ncellx + (i-1);
            neighbour[(ncelly-1)*ncellx + i][4] = (ncelly-1)*ncellx + i;
            neighbour[(ncelly-1)*ncellx + i][5] = (ncelly-1)*ncellx + (i+1);
            neighbour[(ncelly-1)*ncellx + i][6] = -1;
         }
         for (j=1; j<ncelly-1; j++) {
            /* left */
            neighbour[j*ncellx][0] = (j-1)*ncellx;
            neighbour[j*ncellx][1] = (j-1)*ncellx + 1;
            neighbour[j*ncellx][2] = j*ncellx;
            neighbour[j*ncellx][3] = j*ncellx + 1;
            neighbour[j*ncellx][4] = (j+1)*ncellx;
            neighbour[j*ncellx][5] = (j+1)*ncellx + 1;
            neighbour[j*ncellx][6] = -1;
            /* right */
            neighbour[(j+1)*ncellx-1][0] = j*ncellx-2;
            neighbour[(j+1)*ncellx-1][1] = j*ncellx-1;
            neighbour[(j+1)*ncellx-1][2] = (j+1)*ncellx-2;
            neighbour[(j+1)*ncellx-1][3] = (j+1)*ncellx-1;
            neighbour[(j+1)*ncellx-1][4] = (j+2)*ncellx-2;
            neighbour[(j+1)*ncellx-1][5] = (j+2)*ncellx-1;
            neighbour[(j+1)*ncellx-1][6] = -1;
         }
         /* Corners */
         /* Top left */
         neighbour[0][0] = 0;
         neighbour[0][1] = 1;
         neighbour[0][2] = ncellx;
         neighbour[0][3] = ncellx+1;
         neighbour[0][4] = -1;
         /* Top right */
         neighbour[ncellx-1][0] = ncellx-2;
         neighbour[ncellx-1][1] = ncellx-1;
         neighbour[ncellx-1][2] = 2*ncellx-2;
         neighbour[ncellx-1][3] = 2*ncellx-1;
         neighbour[ncellx-1][4] = -1;
         /* Bottom left */
         neighbour[ncellx*(ncelly-1)][0] = ncellx*(ncelly-2);
         neighbour[ncellx*(ncelly-1)][1] = ncellx*(ncelly-2) + 1;
         neighbour[ncellx*(ncelly-1)][2] = ncellx*(ncelly-1);
         neighbour[ncellx*(ncelly-1)][3] = ncellx*(ncelly-1) + 1;
         neighbour[ncellx*(ncelly-1)][4] = -1;
         /* Bottom right */
         neighbour[ncellx*ncelly-1][0] = ncellx*(ncelly-1) - 2;
         neighbour[ncellx*ncelly-1][1] = ncellx*(ncelly-1) - 1;
         neighbour[ncellx*ncelly-1][2] = ncellx*ncelly - 2;
         neighbour[ncellx*ncelly-1][3] = ncellx*ncelly - 1;
         neighbour[ncellx*ncelly-1][4] = - 1;
      }
      printf ("Cell list grid: %ld x %ld\n", ncellx, ncelly);
      printf ("Cell size x:    %lf\n", box->x / ncellx);
      printf ("Cell size y:    %lf\n", box->y / ncelly);

      //create clustersize arrays and variables and set to zero
      clustersize = (long *)malloc(npart * sizeof(long));
      percclusters = (long *)malloc(npart * sizeof(long));
      meanclustersize = meanclustersize_expc = 0;

      //want to make sure they are all zero initially as not guaranteed every element will be accessed in program
      for (i = 0;i < npart; i++){
         clustersize[i] = 0;
         percclusters[i] = 0;
      }

      //set up coordination number array
      coordnum = (long *)malloc(npart * sizeof(long));

      //set percolation statistcs back to zero
      percx = percy = percxyboth = percxandy  = nullstat;

      //open LAMMPS config file for each shell thickness
      configfile = fopen(filename, "r");
      
      //loop through configurations
      for (sweep = 1; sweep <= nsweeps; sweep++){
         //printf("sweep: %ld\n", sweep);
         //skip equilibrate steps
         if (sweep <= equilibrate){
            load(configfile, particle, npart, ncellx, ncelly, cfirst, *box, sourcetype);
            continue;
         }

         //gather stats on particles
         n1 = assign_diameters(npart, particle, fraction, diameter);
         pf = (n1 * disc_area(diameter[0]) + (npart - n1) * disc_area(diameter[1]))
            / ( (box->x) * (box->y) );
         accumulate(&frac1, (double)n1/(double)npart);
         accumulate(&packfrac, pf);

         //reset pointers for cell lists and histogram data
         for (i=0; i<npart; i++) { particle[i].next = NULL; }
         for (i=0; i<ncells; i++) { cfirst[i] = NULL; }

         //printf("Sweep: %ld\n", sweep);      
         load(configfile, particle, npart, ncellx, ncelly, cfirst, *box, sourcetype);

         /* Percolation*/  
         
         //printf("Shell: %lf", shell);
         if (periodic) {
            pcomp = percolate(npart, particle, *box, &nolap, &npc1, &npc2,
               cfirst, neighbour, typetest, shell, clustersize, percclusters, coordnum, fulldir,
               &avecoordnum);
         } else {
            pcomp = spantest(npart, particle, *box, &nolap, &npc1, &npc2,
               cfirst, neighbour, shell, clustersize, percclusters, coordnum, fulldir, 
               &avecoordnum);
         }
         
         //gather statistics on percolation
         accumulate (&contacts, 2.0*nolap/npart);
      
         accumulate (&percx,                pcomp.x); // x only
         accumulate (&percy,                pcomp.y); // y only
         accumulate (&percxyboth, (double)( pcomp.x>0.0 || pcomp.y>0.0)); // x, y or both
         accumulate (&percxandy,  (double)( pcomp.x>0.0 && pcomp.y>0.0)); // x and y directions simultaneously
      
         //write time series data
         fprintf(timeseries, "%d\n", pcomp.x > 0.0 || pcomp.y > 0.0 ? 1 : 0);

         /* Calculate RDF if required */
         if (rdfflag) {
            rdf(particle, rdfhist, npart, nsweeps, equilibrate, numbins, dr, maxsep, *box);
         }

      } //end of configurations loop

      rdfflag = 0;

      cluster_distribution(clustersize, &meanclustersize, &meanclustersize_expc, percclusters, npart, nsweeps, equilibrate, shell, fulldir);
      write_coordinationnumbers(coordnum, shell, npart, fulldir);
      //close LAMMPs file so it can be reset for next shell thickness
      fclose(configfile);

      //write percolation probability for this shell thickness
      fprintf(shellseries, 
         "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
         shell,
         percx.mean,
         percx.rms,
         percy.mean,
         percy.rms,
         percxyboth.mean,
         percxyboth.rms,
         percxandy.mean,
         percxandy.rms,
         meanclustersize,
         meanclustersize_expc,
         avecoordnum.mean);
       
      printf("Average coordination number: %lf\n", avecoordnum.mean);
      printf("Average cluster size: %lf\n", meanclustersize);
      printf("Average cluster size (excluding percolating clusters): %lf\n\n", meanclustersize_expc);

      shell += sstep;
      //free up arrays
      free(cfirst);
      for (i=0; i<ncells; i++) {
         free(neighbour[i]);
      }
      free(neighbour);
      free(clustersize);
   }

   //end of shell increases
   if (rdfhist) {
      normalise_and_write_rdf(rdfhist, numbins, npart, nsweeps, equilibrate, dr, *box, fulldir);
   }

   if (movie > 0) fclose (mf);
   if (movie > 0) fclose (testfile);
   fclose(shellseries);
   fclose(timeseries);   
   
}


/*..............................................................................*/

/*
Accumulate a value into the statistics and update the mean and rms values.
*/
void accumulate(struct mystat *q, double x)
{
   (*q).sum += x;
   (*q).sum2 += x*x;
   (*q).samples++;
   (*q).mean = (*q).sum / (*q).samples;
   (*q).rms = sqrt(fabs((*q).sum2 / (*q).samples -
                        (*q).sum * (*q).sum / (*q).samples / (*q).samples));
}

/*................................................................................*/

/*
Accumulate a batch of values into the statistics and update the mean and rms values.
*/
void baccumulate(struct mystat *q, double x, double x2, long n)
{
   if (LONG_MAX - (*q).samples < n) {
      fprintf (stderr, "ERROR: Overflow of sample number in baccumulate\n");
      exit (1);
   }
   (*q).sum += x;
   (*q).sum2 += x2;
   (*q).samples += n;
   if ((*q).samples > 0) {
      (*q).mean = (*q).sum / (*q).samples;
      (*q).rms = sqrt(fabs((*q).sum2 / (*q).samples -
                        (*q).sum * (*q).sum / (*q).samples / (*q).samples));
   }
}

/*................................................................................*/
