#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>

#define ZEROTOL 1.0e-12     /* Dot products below ZEROTOL are deemed zero */
#define MAXO 100            /* Maximum number of connections per rod */

/* Square of a variable */
#define SQ(x) ((x) * (x))

/* Dot product in 2D */
#define DOT(a,b) ((a).x * (b).x + (a).y * (b).y)

/* Acceptance ratio */
#define RATIO(a) ( ((a).acc+(a).rej) > 0 ? 1.0*(a).acc/((a).acc+(a).rej) : 0.0 )

struct vector {             /* Define a 2D vector structure */
   double x;
   double y;
};

struct disc {              /* Define a disc */
   long idx;               /* Index of each spherocylinder (0,1,2...) for dereferencing pointers */
   struct vector pos;      /* Position vector */
   int species;            /* 0 for species 1 or 1 for species 2 */
   int structure;          /* Defines the type of structure the particle belongs to (see README.md) */
   double diameter;        /* L/D for the cylindrical part */
   long neighbours[10];     /* Array of the neighbours for each particle */
   long cell;              /* Cell in cell-list scheme to which particle belongs */
   struct disc *next;      /* Pointer to next particle in the same cell */
};

struct mystat {               /* Define statistics counters */
   double sum;
   double sum2;
   long samples;
   double mean;
   double rms;
};
#ifndef SEED_DEFINITION
#define SEED_DEFINITION 1 
extern long seed;                  /* Seed for random number generator */
#endif

