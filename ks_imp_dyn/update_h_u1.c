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

//#include "graphene_includes.h"
#include "ks_imp_includes.h"	/* definitions files and prototypes */

void update_h_u1( Real eps ){
  int ff_prec = PRECISION;  /* Just use prevailing precision for now */
#ifdef FN
    free_fn_links_u1(&fn_links_u1);
    free_fn_links_u1(&fn_links_dmdu0);
#endif
    /* gauge field force */
    rephase(OFF);
    imp_gauge_force(eps,F_OFFSET(mom));
    rephase(ON);
    /* fermionic force */
    /* First compute M*xxx in temporary vector xxx_odd */
    /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
    eo_fermion_force_oneterm_u1( eps, ((Real)nflavors)/4., F_OFFSET(xxx),
			      ff_prec, &fn_links_u1, &ks_act_paths );

} /* update_h_u1 */


