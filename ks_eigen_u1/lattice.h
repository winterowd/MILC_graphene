#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD QCD program, version 4
   This file defines global scalars and the fields in the lattice.

   Directory for dynamical improved KS action.  Allow:
	naik term in fermion action
	fat links in fermion action
	arbitrary paths in gauge action (eg Symanzik imp.)

   If "FN" is defined,
     Includes storage for Naik improvement (longlink[4], templongvec[4],
     gen_pt[16], etc.
     Includes storage for "fat links"  (fatlink[4])
*/


#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/generic_ks_u1.h" /* For ferm_links_t and ks_action_paths */
#include "../include/random.h"
#include "../include/io_lat.h"    /* For gauge_file */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

/* The lattice is an array of sites.  */
#define MOM_SITE   /* If there is a mom member of the site struct */
typedef struct {
    /* The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	int index;
#ifdef SITERAND
	/* The state information for a random number generator */
	double_prn site_prn;
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;
#endif

/* ------------------------------------------------------------ */
/*   Now come the physical fields, program dependent            */
/* ------------------------------------------------------------ */
	/* gauge field */
	complex link[4];	/* the fundamental field */
        Real potential[4];
#ifdef HMC_ALGORITHM
 	complex old_link[4];
        Real old_potential[4];
	/* For accept/reject */
#endif

	/* antihermitian momentum matrices in each direction */
 	complex mom[4];

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

	/* 3 element complex vectors */
 	complex phi;		/* Gaussian random source vector */
 	complex resid;	/* conjugate gradient residual vector */
 	complex cg_p;	/* conjugate gradient change vector */
 	complex xxx;		/* solution vector = Kinverse * phi */
 	complex ttt;		/* temporary vector, for K*ppp */
 	complex g_rand;	/* Gaussian random vector*/
	/* Use trick of combining xxx=D^adj D)^(-1) on even sites with
	   Dslash time this on odd sites when computing fermion force */
	
#ifdef PHI_ALGORITHM
 	complex old_xxx;	/* For predicting next xxx */
#endif
#ifdef SPECTRUM
	complex propmat[3];	/* For three source colors */
	complex propmat2[3];	/* nl_spectrum() */
	complex tempmat2;
	/* for spectrum_imp() */
	/**complex quark_source, quark_prop, anti_prop;**/
#define quark_source propmat2[0]
#define quark_prop propmat2[1]
#define anti_prop propmat2[2]
#endif

	/* temporary vectors and matrices */
	complex tempvec[4];	/* One for each direction */
#ifdef FN
	complex templongvec[4];	/* One for each direction */
        complex templongv1;
#endif
	complex tempmat1,staple;
} site;


/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars 
   beta is overall gauge coupling factor (6/g^2)
   mass is quark mass
   u0 is tadpole improvement factor, perhaps (plaq/3)^(1/4)
*/
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
EXTERN	int iseed;		/* random number seed */
EXTERN	int niter,nrestart,nflavors;
EXTERN  Real mass,u0_s,u0_t, v_Fermi;
#ifdef EXT_FIELD
EXTERN  Real N_b;
EXTERN  Real B_ext;
#endif
EXTERN	Real rsqmin,rsqprop;
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, RELOAD_BINARY,
			   RELOAD_CHECKPOINT, FRESH */
EXTERN	char startfile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN	int total_iters;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */
        /* source time, increment for it, and number of source slices */

/* Eigenvalue related global variables */
EXTERN  int Nvecs ; /* number of eigenvectors */
EXTERN  Real eigenval_tol ; /* Tolerance for the eigenvalue computation */
EXTERN  Real error_decr ; /* error decrease per Rayleigh minimization */
EXTERN  int MaxIter ; /* max  Rayleigh iterations */
EXTERN  int Restart ; /* Restart  Rayleigh every so many iterations */
EXTERN  int Kiters ; /* Kalkreuter iterations */
/*******/

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 16
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for definition of the quark action */
EXTERN ferm_links_u1_t    fn_links, fn_links_dmdu0;
EXTERN ks_action_paths ks_act_paths,ks_act_paths_dmdu0;
EXTERN int test_flag;
#endif /* _LATTICE_H */
