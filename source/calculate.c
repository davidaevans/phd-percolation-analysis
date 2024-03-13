#include "global.h"
#include "prototypes.h"

void cluster_distribution(long *clustersize, 
                          double *meanclustersize,
                          double *meanclustersize_expc,
                          long *percclusters,
                          long npart,
                          long nsweeps,
                          long equilibrate,
                          double shell,
                          char *fulldir) {


    double *normclustersize_ns; /* Array for normalised histogram */
    double *normclustersize_sns; /* Array for normalised histogram */

    double *normclustersize_ns_expc; /* Array for normalised histogram excluding percolating clusters*/
    double *normclustersize_sns_expc; /* Array for normalised histogram excluding percolating clusters*/


    long i;
    long statsweeps;  //number of sweeps that count for taking statistics
    long totclusters; //total number of clusters sampled
    long totpercclusters; //total number of percolating clusters 
    long npartperc; //number of particles that belong to a percolating cluster
    long sum, sum_expc;

    normclustersize_ns = (double *)malloc(npart * sizeof(double));
    normclustersize_sns = (double *)malloc(npart * sizeof(double));
    normclustersize_ns_expc = (double *)malloc(npart * sizeof(double));
    normclustersize_sns_expc = (double *)malloc(npart * sizeof(double));
    //set all elements to zero as no guarantee they'll all be visited.
    for (i = 0; i < npart; i++) {
        normclustersize_ns[i] = 0;
        normclustersize_sns[i] = 0;
        normclustersize_ns_expc[i] = 0;
        normclustersize_sns_expc[i] = 0;
    }

    sum = sum_expc = totclusters = totpercclusters = npartperc = 0;
    *meanclustersize = 0;
    statsweeps = nsweeps - equilibrate;

    //get info about cluster distribution for calculating normalisations
    for (i = 0; i < npart; i++) {
        //total number of clusters
        totclusters += clustersize[i];
        totpercclusters += percclusters[i];
        sum += clustersize[i] * SQ(i+1);
        npartperc = percclusters[i] * (i+1);
        sum_expc += (clustersize[i] - percclusters[i]) * SQ(i+1); 

    }
    *meanclustersize = (double) sum/(npart * statsweeps);
    *meanclustersize_expc = (double) sum_expc/((npart*statsweeps)-npartperc);
    printf("Total number of clusters:  %ld\n", totclusters);
    printf("Total number of percolating clusters: %ld\n", totpercclusters);
    

    //normalisation by total number of clusters/particles
    for (i = 0; i < npart; i++){

        //with percolating clusters
        normclustersize_ns[i] = (double) clustersize[i]/totclusters;
        normclustersize_sns[i] = (double) (clustersize[i]*(i+1))/(statsweeps*npart);

        //excluding percolating clusters
        normclustersize_ns_expc[i] = (double) (clustersize[i]-percclusters[i])/(totclusters-totpercclusters);
        normclustersize_sns_expc[i] = (double) ((clustersize[i]-percclusters[i])*(i+1))/((statsweeps*npart)-npartperc);
    }

    write_distribution(shell, clustersize, percclusters, normclustersize_ns, normclustersize_sns, normclustersize_ns_expc, normclustersize_sns_expc, npart, fulldir);
    free(normclustersize_ns);
    free(normclustersize_sns);
}