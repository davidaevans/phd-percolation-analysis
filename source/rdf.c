#include "global.h"
#include "prototypes.h"


// FOLLOWING ALGORITH 7, pg.86 from Frenkel and Smit 2nd ed.
void rdf(struct disc *particle,
         long *rdfhist,
         long npart,
         long nsweeps,
         long equilibrate,
         long numbins,
         double dr,
         double maxsep,
         struct vector box) {

    
    long i,j;
    double sep, sep2;             /* separation of particles */ 
    int bin;

    for (i=0;i<npart-1;i++) {
        for (j=i+1; j<npart;j++){
            sep2 = imagesep(particle[i].pos, particle[j].pos, box);
            sep = sqrt(sep2);

            if(sep < maxsep) {
                bin = sep/dr;
                rdfhist[bin] += 2;
            }

        }
        //printf("i: %ld\n", i);
    }

}

//normalisation of the rdf called at the end
void normalise_and_write_rdf(long *rdfhist,
                             long numbins,
                             long npart,
                             long nsweeps,
                             long equilibrate,
                             double dr,
                             struct vector box,
                             char *fulldir) {

    long i;
    long sweeps;
    long nparteff;
    double *normrdfhist;
    FILE *outfile;
    char rdfname[80];

    normrdfhist = (double *)malloc(numbins * sizeof(double));

    sweeps = nsweeps - equilibrate;
    nparteff = npart * sweeps;
    
    //normalisation of the histogram by expected number of particles from ideal gas
    double volbin, nideal;
    for (i=0;i<numbins;i++){

        volbin = (SQ(i+1) - SQ(i)) * SQ(dr) * M_PI;
        nideal = volbin * npart/(box.x*box.y);
        normrdfhist[i] = rdfhist[i]/(nparteff*nideal);
    }

    //write the resulting normalised histogram to the rdf.dat file 

    sprintf(rdfname, "%s/%s", fulldir, "rdf.dat");
    outfile = fopen(rdfname, "w");
    if (outfile == NULL) die ("Could not open rdf file for writing");
    for (i = 0; i < numbins; i++) {
        fprintf(outfile, "%lf %lf\n",
                (i+0.5)*dr,
                normrdfhist[i]);

        fflush(outfile);
    }
    fclose(outfile);
}

