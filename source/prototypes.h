
void accumulate(struct mystat *q, double x);

double anint(double arg);

long assign_diameters(long npart, struct disc *particle, double fraction,
   double *diameter);

void baccumulate(struct mystat *q, double x, double x2, long n);

void die(char string[]);

void draw(FILE *outfile, struct vector box, long npart,
          struct disc *particle);

long getcell(struct vector pos, long ncellx, long ncelly);

struct vector image(struct vector r1, struct vector r2, struct vector box);

double imagesep(struct vector r1, struct vector r2, struct vector box);

int overlap(struct vector r_cm, double diameter1, double diameter2, double shell);

struct vector percolate(long npart, struct disc *particle,
                        struct vector box, long *ctot,
                        struct mystat *npc1, struct mystat *npc2,
                        struct disc **cfirst, long **neighbour,
                        long typetest, double shell, long *clustersize,
                        long *percclusters, long *coordnum, char *fulldir,
                        struct mystat *avecoordnum);

struct vector spantest(long npart, struct disc *particle,
                       struct vector box, long *ctot,
                       struct mystat *npc1, struct mystat *npc2,
                       struct disc **cfirst, long **neighbour,
                       double shell,long *clustersize, long *percclusters,
                       long *coordnum, char *fulldir,
                       struct mystat *avecoordnum);

double ran2(long *idum);

void read_options(long *npart, struct vector *box, double *diameter, long *nsweeps,
   long *report, long *movie,
   double *fraction, int *fmethod, int *periodic,
   long *configs, long *typetest,
   double *smin, double *smax, double *sstep,
   long *equilibrate, char *filename,
   double *dr, int *rdfflag, long *sourcetype, long *nblocks);

double disc_area(double diameter);

void simulate(long npart, double fraction, struct vector *box, double *diameter,
   long nsweeps, long report, long movie,
   int periodic, struct disc *particle, long configs, long typetest,
   double smin, double smax, double sstep, long equilibrate, char *filename, 
   double dr, int rdfflag, long sourcetype, char *fulldir);

void touching(long npart, struct disc *particle, struct vector box,
              long *nc, long **conn, long *ctot,
              struct disc **cfirst, long **neighbour, double shell);

void touching1(long npart, struct disc *particle, struct vector box,
              long *nc, long **conn, long *ctot,
              struct disc **cfirst, long **neighbour, long typetest,
              double shell);

void touching_nopbc(long npart, struct disc *particle, struct vector box,
              long *nc, long **conn, long *ctot,
              struct disc **cfirst, long **neighbour, double shell);

void write_config(FILE *testfile, struct vector box, long npart,
                  struct disc *particle, long sweep, int percolating);

void cluster_distribution(long *clustersize, double *meanclustersize, 
                        double *meanclustersize_expc,
                        long *percclusters,
                        long npart, long nsweeps, long equilibrate, double shell,
                        char *fulldir);

void write_distribution(double shell, long *clustersize, long *percclusters, 
                        double *normclustersize_ns,
                        double *normclustersize_sns,
                        double *normclustersize_ns_expc,
                        double *normclustersize_sns_expc, long npart, char *fulldir);


void rdf(struct disc *particle,
         long *rdfhist,
         long npart,
         long nsweeps,
         long equilibrate,
         long numbins,
         double dr,
         double maxsep,
         struct vector box);

void normalise_and_write_rdf(long *rdfhist,
                             long numbins,
                             long npart,
                             long nsweeps,
                             long equilibrate,
                             double dr,
                             struct vector box,
                             char *fulldir);

void load(FILE *configfile, struct disc *particle, long npart,
         long ncellx, long ncelly, struct disc **cfirst, struct vector box,
         long sourcetype);

void load_LAMMPS(FILE *configfile, struct disc *particle, long npart,
         long ncellx, long ncelly, struct disc **cfirst, struct vector box);

void load_MC(FILE *configfile, struct disc *particle, long npart,
         long ncellx, long ncelly, struct disc **cfirst, struct vector box);

void write_coordinationnumbers(long *coordnum, double shell, long npart, char *fulldir);