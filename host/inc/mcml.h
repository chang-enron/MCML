#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#define PI 3.1415926
#define WEIGHT 1E-4     /* Critical weight for roulette. */
#define CHANCE 0.1      /* Chance of roulette survival. */
#define STRLEN 256      /* String length. */

#define Boolean char

#define SIGN(x) ((x)>=0 ? 1:-1)
#define MaxSteps  4096
#define MaxLayers 10


#define USE_FLOAT 1
#if USE_FLOAT
    #define REAL float
#else
    #define REAL double
#endif

typedef struct {
    REAL z0, z1;    /* z coordinates of a layer. [cm] */
    REAL n;         /* refractive index of a layer. */
    REAL mua;       /* absorption coefficient. [1/cm] */
    REAL mus;       /* scattering coefficient. [1/cm] */
    REAL g;         /* anisotropy. */

    REAL cos_crit0, cos_crit1;  
} LayerStruct;

typedef struct {
    char   out_fname[STRLEN]; /* output file name. */
    char   out_fformat;       /* output file format. */

    long     num_photons;       /* to be traced. */
    REAL Wth;               /* play roulette if photon */
    REAL dz;                /* z grid separation.[cm] */ 
    REAL dr;                /* r grid separation.[cm] */
    REAL da;                /* alpha grid separation. */
    short nz;                   /* array range 0..nz-1. */
    short nr;                   /* array range 0..nr-1. */
    short  na;                   /* array range 0..na-1. */
    short   num_layers;         /* number of layers. */
    LayerStruct layerspecs[MaxLayers];  /* layer parameters. */ 
} cl_InputStruct;
typedef struct {
    REAL x, y ,z;   /* Cartesian coordinates.[cm] */
    REAL ux, uy, uz;/* directional cosines of a photon. */
    REAL w;         /* weight. */
    Boolean dead;       /* 1 if photon is terminated. */
    REAL s;         /* current step size. [cm]. */
    REAL sleft;     /* step size left. dimensionless [-]. */
    unsigned char layer;        /* index to layer where the photon */
} cl_PhotonStruct;

typedef struct {
    unsigned char x, y;
    REAL  w;
} trace;

typedef struct {
    short  total_steps;
    Boolean  Rd_valid;
    Boolean  Tt_valid;
    trace  data[MaxSteps]; // for Drop
    trace  Rdra; // for RecordR
    trace  Ttra;    // for RecordT
} tmpOutStruct;

typedef struct {
    REAL    Rsp;    /* specular reflectance. [-] */
    REAL ** Rd_ra;  /* 2D distribution of diffuse */
    /* reflectance. [1/(cm2 sr)] */
    REAL *  Rd_r;   /* 1D radial distribution of diffuse */
    /* reflectance. [1/cm2] */
    REAL *  Rd_a;   /* 1D angular distribution of diffuse */
    /* reflectance. [1/sr] */
    REAL    Rd;     /* total diffuse reflectance. [-] */

    REAL ** A_rz;   /* 2D probability density in turbid */
    /* media over r & z. [1/cm3] */
    REAL *  A_z;    /* 1D probability density over z. */
    /* [1/cm] */
    REAL *  A_l;    /* each layer's absorption */
    /* probability. [-] */
    REAL    A;      /* total absorption probability. [-] */

    REAL ** Tt_ra;  /* 2D distribution of total */
    /* transmittance. [1/(cm2 sr)] */
    REAL *  Tt_r;   /* 1D radial distribution of */
    /* transmittance. [1/cm2] */
    REAL *  Tt_a;   /* 1D angular distribution of */
    /* transmittance. [1/sr] */
    REAL    Tt;     /* total transmittance. [-] */
} OutStruct;


typedef struct {
    char   out_fname[STRLEN]; /* output file name. */
    char   out_fformat;       /* output file format. */
    long     num_photons;       /* to be traced. */
    REAL Wth;               /* play roulette if photon */
    REAL dz;                /* z grid separation.[cm] */ 
    REAL dr;                /* r grid separation.[cm] */
    REAL da;                /* alpha grid separation. */
    short nz;                   /* array range 0..nz-1. */
    short nr;                   /* array range 0..nr-1. */
    short na;                   /* array range 0..na-1. */
    short   num_layers;         /* number of layers. */
    LayerStruct *layerspecs;  /* layer parameters. */ 
} InputStruct;

typedef struct {
    REAL x, y ,z;   /* Cartesian coordinates.[cm] */
    REAL ux, uy, uz;/* directional cosines of a photon. */
    REAL w;         /* weight. */
    Boolean dead;       /* 1 if photon is terminated. */
    REAL s;         /* current step size. [cm]. */
    REAL sleft;     /* step size left. dimensionless [-]. */
    short layer;
} PhotonStruct;

REAL  *AllocVector(short, short);
REAL **AllocMatrix(short, short,short, short);
void     FreeVector(REAL *, short, short);
void     FreeMatrix(REAL **, short, short, short, short);
void     nrerror(char *);

