#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

/*
Reads the run parameters from the external file "options".  See INSTRUCTIONS
for a list of keywords.
*/

void read_options(long *npart, struct vector *box, double *diameter, long *nsweeps,
   long *report, long *movie,
   double *fraction, int *fmethod, int *periodic,
   long *configs, long *typetest, 
   double *smin, double *smax, double *sstep,
   long *equilibrate, char *filename,
   double *dr, int *rdfflag, long *sourcetype, long *nblocks)
{
   int i;
   char command[20];
   char option[20];
   char error[200];
   char fname[80]; //name of file where configurations are read from
   double cellarea;   /* Area of simulation cell */
   double minbox;
   double fractionl;  /* Fraction of species 1 by cylinder length */
   double fractionv;  /* Fraction of species 1 by area */
   double pf;         /* Target packing fraction */
   FILE *infile = NULL;
   FILE *datafile = NULL; //LAMMPS file


   /* Prototypes for keywords library */
   int  read_line(FILE *);
   int  get_string(char [], int);
   void upper_case(char []);
   int  get_int(long int *);
   int  get_double(double *);


   /*--- 0. Defaults ---*/
   box->x = box->y -1.0;              /* Box dimensions */
   *fmethod = 0;                      /* 0=fraction of rods by number, 1=fraction by area */
   *fraction = 1.0;                   /* Fraction of species 1 rods */
   diameter[0] = 0.0;                   /* L/D aspect ratio of cylindrical part of rod species 1 */
   diameter[1] = 0.0;                   /* L/D aspect ratio of cylindrical part of rod species 2 */
   *npart = 10;                       /* Number of rods */
   *nsweeps = 100;                    /* Number of accumulation sweeps */
   *periodic = 1;                     /* Wrapping rather than spanning percolation criterion */
   *report = 100;                     /* Number of sweeps between statistics reports */
   seed = -1;                         /* Random number seed */
   *configs = 0;                       /* Number of sweeps before writing configs for debugging */
   *smin = 0.0;                       /* minimum connectivity shell thickness */
   *smax = 0.5;                       /* maximum connectivity shell thickness */
   *sstep = 0.05;                     /* step in connectivity shell thickness */
   *typetest = 0;                     /* Default to test for all particles of all type */
   *equilibrate = 0;                  /* Default equilibration sweeps is 0 */
   *dr = 0;                           /* Width of bins for RDF histogram */
   *rdfflag = 0;                      /* Boolean flag whether to perform RDF calc or not - 0 = no, 1 = yes */
   *sourcetype = 0;                   /* 0 = LAMMPS, 1 = MC */
   *nblocks = 1;                      /* Number of blocks for average */

   /*--- 1. Read in values ---*/

   infile = fopen("options", "r");
   if (infile == NULL) die ("Could not open \"options\" file");
   printf ("Reading run options from the \"options\" file\n\n");

   while ( read_line(infile) ) {
      get_string(command, sizeof(command));
      upper_case(command);

      if (*command == '#') {
         continue;

      //} else if (strcmp(command, "BOX") == 0) {
         //if (!get_double(&(box->x))) die ("Could not read x box dimension after BOX");
         //if (!get_double(&(box->y))) die ("Could not read y box dimension after BOX");
      
      } else if (strcmp(command, "FILENAME") == 0) {
         if(!get_string(fname, sizeof(fname))) die ("Could not read filename after FILENAME");
         strcpy(filename, fname);
         //sprintf(error, "fname: %s, filename: %s", fname, filename);
         //die(error);
      } else if (strcmp(command, "DEFINITION") == 0) {
         get_string(option, sizeof(option));
         upper_case(option);
         if (strcmp(option, "SPANNING") == 0) {
            *periodic = 0;
         } else if (strcmp(option, "WRAPPING") == 0) {
            *periodic = 1;
         } else {
            sprintf (error, "Unrecognised option after DEFINITION keyword: %s", option);
            die (error);
         }

      } else if (strcmp(command, "DIAMETER") == 0) {
         if (!get_double(&diameter[0])) die ("Could not read diameter of species 1 after DIAMETER");
         if (!get_double(&diameter[1])) die ("Could not read diameter of species 2 after DIAMETER");
      } else if (strcmp(command, "EQUILIBRATE") == 0){
         if (!get_int(equilibrate)) die ("Could not read equilibration sweeps after EQUILIBRATE");
      } else if (strcmp(command, "MOVIE") == 0) {
         if (!get_int(movie)) die ("Could not read number sweeps between movie frames after MOVIE");

      //} else if (strcmp(command, "DISCS") == 0) {
         //if (!get_int(npart)) die ("Could not read number of discs after DISCS");

      } else if (strcmp(command, "SEED") == 0) {
         if (!get_int(&seed)) die ("Could not read random number seed after SEED");
         seed = -abs(seed);
      } else if (strcmp(command, "SHELL") == 0){
         if (!get_double(smin)) die ("Could not read minimum shell thickness after SHELL");
         if (!get_double(smax)) die ("Could not read maximum shell thickness after SHELL");
         if (!get_double(sstep)) die ("Could not read step in shell thickness after SHELL");      
      } else if (strcmp(command, "SPECIES1") == 0) {
         get_string(option, sizeof(option));
         upper_case(option);
         if (strcmp(option, "NUMBER") == 0) {
            *fmethod = 0;
            if (!get_double(fraction)) die ("Could not read fraction after SPECIES1 NUMBER");
         } else if (strcmp(option, "AREA") == 0) {
            *fmethod = 1;
            if (!get_double(&fractionv)) die ("Could not read fraction after SPECIES1 AREA");
         } else if (strcmp(option, "LENGTH") == 0) {
            *fmethod = 2;
            if (!get_double(&fractionl)) die ("Could not read fraction after SPECIES1 LENGTH");
         } else {
            sprintf (error, "Unrecognised option after SPECIES1 keyword: %s", option);
            die (error);
         }

      } else if (strcmp(command, "STATISTICS") == 0) {
         if (!get_int(report)) die ("Could not read number sweeps between statistics reports after STATISTICS");
      } else if (strcmp(command, "BLOCKS") == 0) {
         if (!get_int(nblocks)) die ("Could not read number of blocks after BLOCKS");

      } else if (strcmp(command, "SWEEPS") == 0) {
         if (!get_int(nsweeps)) die ("Could not read number of sweeps after SWEEPS");

      } else if (strcmp(command, "DEBUG") == 0) {
         if (!get_int(configs)) die ("Could not read number of sweeps after SWEEPS");
      
      } else if (strcmp(command, "TYPE") == 0) {
         if (!get_int(typetest)) die ("Could not read structure type after TYPE");
      } else if (strcmp(command, "SOURCETYPE") == 0) {
         if (!get_int(sourcetype)) die("Could not read source file type (LAMMPS or MC).");
      } else if (strcmp(command, "BINSIZE") == 0) {
         if (!get_double(dr)) die ("Could not read bin size after BINSIZE");
         *rdfflag = 1;
      } else {
         sprintf (error, "Unrecognised keyword: %s", command);
         die (error);
      }
   }

   fclose(infile);

   datafile = fopen(filename, "r");
   if (datafile == NULL) die ("Could not open configurations file");
   
   for (i = 0; i < 9; i++) {
      read_line(datafile);

      if (i==3) {
         if(!get_int(npart)) die ("Could not read number of particles from LAMMPS configuration file");
      }

      if (i==5) {
         double tmp_xa, tmp_xb;
         if (!get_double(&tmp_xa)) die ("Could not read first box size in x direction from LAMMPS configuration file");
         if (!get_double(&tmp_xb)) die ("Could not read second box size in x direction from LAMMPS configuration file");
         box->x = tmp_xb - tmp_xa;
      }

      if (i==6) {
         double tmp_ya, tmp_yb;
         if (!get_double(&tmp_ya)) die ("Could not read first box size in y direction from LAMMPS configuration file");
         if (!get_double(&tmp_yb)) die ("Could not read second box size in y direction from LAMMPS configuration file");
         box->y = tmp_yb - tmp_ya;
      }
   }

   fclose(datafile);
   
   

   /*--- 2. Validity checks ---*/

   if (*npart < 1) {
      die ("The number of discs must be at least 1.");
   }

   if (*equilibrate < 0) {
      die ("Equilibration time must be greater than or equal to 0.");
   }

   if (diameter[0] < 0.0 || diameter[1] < 0.0) {
      die ("The diameter of discs cannot be negative.");
   }

   if (*smin < 1 ||  *smax < 1) {
      die ("Minimum/maximum shell diameter must be greater than hard core diameter.");
   }

   if (*sstep <= 0) {
      die ("Shell diameter step size must be greater than 0.");
   }

   if (*nblocks <= 0) {
      die ("Nblocks must be greater than 0");
   }

   if (*typetest < 0 || *typetest > 6) {
      die ("Structure type must be an integer between 0 and 6 inclusive.");
   }

   if (*fraction < 0.0 || *fraction > 1.0) {
      die ("The fraction of discs of species 1 must lie in the range 0 to 1.");
   }

   if (*dr < 0) {
		die ("Bin size must be greater than 0.");
	}

   if (seed == 0) {
      die ("The random seed must be a negative integer (not zero).");
   }

   minbox = diameter[0] * 2.0;
   if (diameter[1] > diameter[0]) {
      minbox = diameter[1] * 2.0;
   }
   if ( (box->x < minbox) || (box->y < minbox) ) {
      die ("Both box lengths must be at least two full diameters long for the larger species.");
   }

   if (*report > *nsweeps) *report=*nsweeps;

//FINITE THICKNESS
   if (*fmethod == 0) {
      //fraction of species 1 by area
      fractionv = *fraction * disc_area(diameter[0]) /
         (*fraction * disc_area(diameter[0]) + (1.0 - *fraction) * disc_area(diameter[1]));
      //fraction of species 1 by diameter
      fractionl = *fraction * diameter[0] /
         (*fraction * diameter[0] + (1.0 - *fraction) * diameter[1]);
   } else if (*fmethod == 1) {
      *fraction = (disc_area(diameter[1]) * fractionv) /
         (disc_area(diameter[0]) - fractionv * disc_area(diameter[0]) +
         fractionv * disc_area(diameter[1]));
      fractionl = *fraction * diameter[0] /
         (*fraction * diameter[0] + (1.0 - *fraction) * diameter[1]);
   } else {
      *fraction = (diameter[1] * fractionl) /
         (diameter[0] - fractionl * diameter[0] + fractionl * diameter[1]);
      fractionv = *fraction * disc_area(diameter[0]) /
         (*fraction * disc_area(diameter[0]) + (1.0 - *fraction) * disc_area(diameter[1]));
   }


   cellarea = (box->x) * (box->y) ;
   pf = (*fraction * disc_area(diameter[0]) + (1.0 - *fraction) * disc_area(diameter[1]))
        * *npart/cellarea;

   /*--- 3. Summarize results on standard output ---*/
   printf (" Configurations file:                      %s\n", filename);
   if (*typetest == 0) {
      printf (" Analysing clusters of type:               all\n");
   } else {
      printf (" Analysing clusters of type:               %ld\n", *typetest);
   }
   printf (" Simulation cell dimensions:               %.8lf, %.8lf\n",
           box->x, box->y);
   printf (" Total number of discs:                    %ld\n", *npart);
   printf (" Diameters of each species:                %.8lf, %.8lf\n", diameter[0], diameter[1]);

   printf (" Target fraction of species 1\n");
   printf ("    by number:                             %.10lf\n", *fraction);
   printf ("    weighted by particle area:             %.10lf\n", fractionv);
   printf ("    weighted by disc diameter:             %.10lf\n", fractionl);
   printf (" Overall number density:                   %.10le\n", *npart/cellarea);
   printf (" Target number density of species:         %.10le, %.10le\n",
      *fraction * *npart/cellarea, (1 - *fraction) * *npart/cellarea);
   printf (" Target overall packing fraction:          %.10le\n", pf);
   printf (" Target packing fraction of species:       %.10le, %.10le\n",
      *fraction * disc_area(diameter[0]) * *npart/cellarea,
      (1.0 - *fraction) * disc_area(diameter[1]) * *npart/cellarea);
   //sweeps
   printf (" Total sweeps in data:                     %ld\n", *nsweeps);
   printf (" Equilibration sweeps:                     %ld\n", *equilibrate);
   printf (" Total sweeps for statistics:              %ld\n", *nsweeps - *equilibrate);
   printf (" Sweeps between statistics reports:        %ld\n", *report);
   //shell
   printf (" Minimum shell thickness:                  %lf\n", *smin);
   printf (" Maximum shell thickness:                  %lf\n", *smax);
   printf (" Step in shell thickness:                  %lf\n", *sstep);
   printf (" Number of blocks for averaging:           %ld\n", *nblocks);
   printf (" Random number seed:                       %ld\n", seed);
   if (*movie > 0) {
      printf (" Sweeps between movie frames:              %ld\n", *movie);
   } else {
      printf (" No movie\n");
   }
   printf (" Percolation definition:                   %s\n", *periodic?"WRAPPING":"SPANNING");
   printf ("\n");

}

/*..............................................................................*/
