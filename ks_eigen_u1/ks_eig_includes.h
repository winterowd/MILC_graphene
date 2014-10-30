/****************** ks_eig_includes.h ******************************/
/*
*  Include files for Kogut-Susskind dynamical improved action application
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/io_scidac.h"
#include "../include/generic.h"
#include "../include/generic_ks_u1.h"
#include "../include/dirs.h"
#include "../include/generic_ks_u1_macros.h"

#ifdef FN
#define dslash_site dslash_fn_site
#define dslash_field dslash_fn_field
#endif
#ifdef EO
#define dslash_site dslash_eo_site
#define dslash_field dslash_eo_field
#endif

/* prototypes for functions in high level code */
int setup(void);
int readin(int prompt);
void measure_chirality_u1(complex *src, double *chirality, int parity);

EXTERN gauge_file *startlat_p;
