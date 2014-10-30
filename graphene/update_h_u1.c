/****** update_h_u1.c  -- ******************/
/* updates momentum matrices for improved U(1) action */
/* D.T. & J.H., naik term    8/96
*  D.T., fat link fermion term 5/97
*  D.T. general quark action 1/98
*  D.T. two types of quarks 3/99
*  T.D. and A.H. improved gauge updating spliced in 5/97
*  C.W. updated to U(1) pending fermion force 7/10
*  C.W. fermion force completed working on imp_gauge_force 10/10
* MIMD version 7 */

#include "graphene_includes.h"
//#include "ks_imp_includes.h"	/* definitions files and prototypes */

void update_h_u1( Real eps ){
  int ff_prec = PRECISION;  /* Just use prevailing precision for now */
  int i, dir;
  site *s;
#ifdef FN
    free_fn_links_u1(&fn_links);
    free_fn_links_u1(&fn_links_dmdu0);
#endif
    /* gauge field force */
    rephase(OFF);
#ifdef NON_COMPACT
    gauge_force(eps); 
#else
    imp_gauge_force(eps,F_OFFSET(mom)); //compact gauge force
#endif
    rephase(ON);
    /* fermionic force */
    /* First compute M*xxx in temporary vector xxx_odd */
    /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
    /*for(dir=XUP;dir<=TUP;dir++)
      FORALLMYSITES(i,s) {
	if(dir!=TUP && s->link[dir].real != 1.)
	  printf("spatial link[%d], numbr. %d is non-zero\n");
	
	  }*/
    eo_fermion_force_oneterm_u1( eps, ((Real)nflavors)/4., F_OFFSET(xxx),
      ff_prec, &fn_links, &ks_act_paths );
#ifdef DRUT_DEBUG
    for(dir=XUP;dir<=TUP;dir++)
      FORALLSITES(i,s) {
	if(dir != TUP) {
	/*printf("spatial link[%d], numbr. %d is non-zero (%e,%e)\n", dir, i, 
	  s->link[dir].real, s->link[dir].imag);*/
	  s->mom[dir].imag = 0.;
	}
      }
#endif
    /* keep N_f/4 as in graphene can just keep N_f=4 for unquenched and 
       N_f=0 for quenched C.W. 7/12/2011 */
    
} /* update_h_u1 */


