#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice.

   Directory for dynamical improved KS action.  Allow:
	arbitrary paths in quark action 
	arbitrary paths in gauge action (eg Symanzik imp.)

   If "FN" is defined,
     Includes storage for Naik improvement (longlink[4], templongvec[4],
     gen_pt[16], etc.

     MODIFIED FOR U(1) by C.W. 10/10
*/

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/generic_ks_u1.h" /* For fn_links and ks_act_paths */

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
        Real temp1;  /* work space for non-compact gauge force 03/11 */
        Real stple;

#ifdef EXT_FIELD
        complex ext_link_fat[4];
        complex ext_link_lng[4];
        Real ext_potential_fat[4];
        Real ext_potential_lng[4];
#endif

#if  defined(HYBRIDS)
	complex tmplink[4];	/* three link straight paths */
#endif

#ifdef HMC_ALGORITHM
 	complex old_link[4];
        Real old_potential[4];
	/* For accept/reject */
#endif

	/* purely complex momentum matrices in each direction */
 	complex mom[4];

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

	/* 3 element complex vectors */
#ifdef ONEMASS
 	complex phi;		/* Gaussian random source vector */
        complex theta;          /* Auxiliary source vector for chiral susc. */
 	complex xxx;		/* solution vector = Kinverse * phi */
        complex yyy;            /* solution vector = Kinverse * xxx */
#ifdef PHI_ALGORITHM
 	complex old_xxx;	/* For predicting next xxx */
#endif
#else
 	complex phi1;	/* Gaussian random source vector, mass1 */
 	complex phi2;	/* Gaussian random source vector, mass2 */
 	complex xxx1;	/* solution vector = Kinverse * phi, mass1 */
 	complex xxx2;	/* solution vector = Kinverse * phi, mass2 */
#endif
 	complex resid;	/* conjugate gradient residual vector */
 	complex cg_p;	/* conjugate gradient change vector */
 	complex ttt;		/* temporary vector, for K*ppp */
 	complex g_rand;	/* Gaussian random vector*/
	/* Use trick of combining xxx=D^adj D)^(-1) on even sites with
	   Dslash times this on odd sites when computing fermion force */
	
#ifdef HYBRIDS
        complex field_strength[6];
#endif
#ifdef SPECTRUM
        complex tempmat1,staple;
	complex propmat;	/* For three source colors ???? */
	complex propmat2[3];	/* nl_spectrum() */
	complex tempmat2;
	/* for spectrum_imp() */
	/**complex quark_source, quark_prop, anti_prop;**/
  //#define quark_source propmat2[0]
  //#define quark_prop propmat2[1]
  //#define anti_prop propmat2[2]
#endif

	/* temporary vectors and matrices */
	complex tempvec[4];	/* One for each direction */
#ifdef FN
	complex templongvec[4];	/* One for each direction */
        complex templongv1;
#endif
#ifdef DM_DU0
	complex dMdu_x;	/* temp vector for <psi-bar(dM/du0)psi>
				   calculation in 'f_meas.c' */
#endif
#ifdef NPBP_REPS
 	complex M_inv;	/* temp vector for M^{-1} g_rand */
#endif
#if defined(CHEM_POT) || defined(D_CHEM_POT)
 	complex dM_M_inv;	/* temp vector for dM/dmu M^{-1} g_rand */
        complex deriv[6];
#endif
} site;

/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars 
   beta is overall gauge coupling factor
   mass, mass1 and mass2 are quark masses
   u0 is tadpole improvement factor, perhaps (plaq/3)^(1/4)
*/
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
#ifdef FIX_NODE_GEOM
EXTERN  int node_geometry[4];  /* Specifies fixed "nsquares" (i.e. 4D
			    hypercubes) for the compute nodes in each
			    coordinate direction.  Must be divisors of
			    the lattice dimensions */
#ifdef FIX_IONODE_GEOM
EXTERN int ionode_geometry[4]; /* Specifies fixed "nsquares" for I/O
			     partitions in each coordinate direction,
			     one I/O node for each square.  The I/O
			     node is at the origin of the square.
			     Must be divisors of the node_geometry. */
#endif
#endif
EXTERN	int iseed;		/* random number seed */
EXTERN	int warms,trajecs,steps,niter,nrestart,propinterval;
EXTERN  int npbp_reps_in;
EXTERN  int prec_pbp;  /* Precisiong of pbp measurements */
EXTERN  int dyn_flavors[MAX_DYN_MASSES]; 
#ifdef ONEMASS
EXTERN  int nflavors;
#else
EXTERN	int nflavors1,nflavors2;  /* number of flavors of types 1 and 2 */
#endif
EXTERN  int nlight_flavors;
EXTERN	Real epsilon;
EXTERN  Real beta,u0_s,u0_t,v_Fermi;
#ifdef EXT_FIELD
EXTERN  Real N_b;
EXTERN  Real B_ext;
#endif
EXTERN  int test_flag;
EXTERN  int n_dyn_masses; // number of dynamical masses
#ifdef ONEMASS
EXTERN  Real mass;
#else
EXTERN  Real mass1,mass2;
#endif
EXTERN	Real rsqmin,rsqprop;
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, RELOAD_BINARY,
			   RELOAD_CHECKPOINT, FRESH */
EXTERN	int saveflag;	/* do with lattice: FORGET, SAVE, SAVE_BINARY,
			   SAVE_CHECKPOINT */
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
#ifdef NON_COMPACT
EXTERN	char startfile2[MAXFILENAME],savefile2[MAXFILENAME];
#endif
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int total_iters;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */
EXTERN  int source_start, source_inc, n_sources;
        /* source time, increment for it, and number of source slices */
EXTERN  char spectrum_request[MAX_SPECTRUM_REQUEST]; /* request list for spectral measurements */
/* parameters for spectrum_multimom */
EXTERN  int spectrum_multimom_nmasses;
EXTERN  Real spectrum_multimom_low_mass;
EXTERN  Real spectrum_multimom_mass_step;
/* parameters for fpi */
EXTERN  int fpi_nmasses;
EXTERN  Real fpi_mass[MAX_FPI_NMASSES];

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN gauge_file *startlat_p;

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
EXTERN ferm_links_u1_t        fn_links;
EXTERN ks_action_paths ks_act_paths;
EXTERN ferm_links_u1_t        fn_links_dmdu0;
EXTERN ks_action_paths ks_act_paths_dmdu0;

/* debugging variables */
EXTERN int num_accepts;
EXTERN int num_rejects;

#endif /* _LATTICE_H */
