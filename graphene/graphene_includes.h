/****************** graphene_includes.h ******************************/
/*
*  Include files for Kogut-Susskind dynamical improved action application
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/io_scidac.h"
#include "../include/generic_ks_u1.h"
#include "../include/generic.h"
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

/* prototypes for functions in this directory */
int setup(void);
int readin(int prompt);
int update_u1();
void update_h_u1( Real eps );
void update_u_u1( Real eps );
double hmom_action_u1( );
double fermion_action_u1( );


/* d_action_u1.c */
double d_action_u1();
void gauge_field_copy_u1(field_offset src,field_offset dest);
void gauge_field_copy_nc(field_offset src,field_offset dest);
double fermion_action_u1();
double hmom_action_u1(void);

void f_measure( field_offset phi_off, field_offset xxx_off, Real mass );
void gauge_field_copy_u1(field_offset src,field_offset dest);
void clear_latvec_u1(field_offset v,int parity);
void copy_latvec_u1(field_offset src,field_offset dest,int parity);
void scalar_mult_add_latvec_u1(field_offset src1,field_offset src2,
			    Real scalar,field_offset dest,int parity);
void scalar2_mult_add_latvec_u1(field_offset src1,Real scalar1,
			     field_offset src2,Real scalar2,
			     field_offset dest,int parity);
void scalar_mult_latvec_u1(field_offset src,Real scalar,
			field_offset dest,int parity);
void rephase(int flag);


