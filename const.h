#ifndef  CONST_INC
#define  CONST_INC

#define Qe  1.602176487E-19
#define EPSILON 0.000000001

#define SQR(a) ((a)*(a))
#define CUB(a) ((a)*(a)*(a))
#define MAX( x , y ) ( (x) > (y) ? (x) : (y) )
#define TWOPI (2.*M_PI)

/* IMPORTANT PARAMETERS FOR THE RUNS!! */
#define NUM_REDEP 10   /* was 1?  10!! */
#define SUBC 1        /* was 10!! */

#ifdef DEBUG
#define ASSERT(a) assert(a)
#else
#define ASSERT(a)
#endif

#endif   /* ----- #ifndef CONST_INC  ----- */
