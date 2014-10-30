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

     MODIFIED FOR U(1) by C.W. 06/10
*/

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/generic_ks_u1.h" /* For ferm_links_t and ks_action_paths */
#include "../include/random.h"    /* For double_prn */
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */ /*Why did DeTar
							 delete this? */
#include "params.h"

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

/* The lattice is an array of sites.  */
typedef struct {
    /* The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	int index;
	/* The state information for a random number generator */
	double_prn site_prn;
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;

/* ------------------------------------------------------------ */
/*   Now come the physical fields, program dependent            */
/* ------------------------------------------------------------ */
	/* gauge field */
	complex link[4];	/* the fundamental field */
        complex old_link[4];   	/* For accept/reject */
  
        Real potential[4];         /* Scalar potential for fundamental field */
        Real old_potential[4];     /* For accept/reject */

	/* purely complex momentum matrices in each direction */
 	complex mom[4];

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

	/*complex numbers */
 	complex phi[MAX_N_PSEUDO]; /* Gaussian random source, each pseudoferm */
 	complex g_rand;	/* Gaussian random vector*/
	/* Use trick of combining xxx=D^adj D)^(-1) on even sites with
	   Dslash times this on odd sites when computing fermion force */
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
   dyn_flavors are the number of flavors renormalizing the gauge action 
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
EXTERN  params par_buf;
EXTERN	int iseed;		/* random number seed */
EXTERN  Real beta,u0;
EXTERN  Real v_Fermi;
EXTERN	int warms,trajecs,steps,niter,nrestart,propinterval;
EXTERN  int niter_md[MAX_N_PSEUDO], niter_fa[MAX_N_PSEUDO], niter_gr[MAX_N_PSEUDO];
EXTERN  int prec_md[MAX_N_PSEUDO], prec_fa[MAX_N_PSEUDO], prec_gr[MAX_N_PSEUDO];
EXTERN  int prec_ff;
EXTERN	Real epsilon;
EXTERN	Real rsqmin_md[MAX_N_PSEUDO], rsqmin_fa[MAX_N_PSEUDO], rsqmin_gr[MAX_N_PSEUDO];
EXTERN  Real rsqprop;
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, RELOAD_BINARY,
			   RELOAD_CHECKPOINT, FRESH */
EXTERN	int saveflag;	/* do with lattice: FORGET, SAVE, SAVE_BINARY,
			   SAVE_CHECKPOINT */
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN	int total_iters;
        /* source time, increment for it, and number of source slices */

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
EXTERN ferm_links_t        fn_links;
EXTERN ks_action_paths     ks_act_paths;

ferm_links_u1_t fn_links_u1; //added for debugging matrix inversion for U(1)

EXTERN int n_pseudo;
EXTERN int phases_in; /* 1 if KS and BC phases absorbed into matrices in site structure */

#endif /* _LATTICE_H */
