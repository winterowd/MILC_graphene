/*************************** ranmom_u1.c *******************************/
/* MIMD version 7 */
/* Added 07/10 C.W. */
/* Produce Gaussian random momenta for the gauge fields. */

#include "generic_includes.h"
#include "../include/generic_ks_u1_macros.h"
#include <defines.h>                 /* For SITERAND */

void ranmom_u1(){
register int i,dir;
register site *s;
    FORALLSITES(i,s){
      //node0_printf("site (%d, %d, %d, %d) \n",s->x,s->y,s->z,s->t);
      for(dir=XUP;dir<=TUP;dir++){
#ifdef SCHROED_FUN
	    if(dir==TUP || s->t>0){
#endif
#ifdef SITERAND
	      /*	random_anti_hermitian( (anti_hermitmat *)&(s->mom[dir]),
			&(s->site_prn) ); */
	      s->mom[dir].real = 0.0;
	      s->mom[dir].imag = gaussian_rand_no_u1(&s->site_prn);
	      //node0_printf("SITERAND!!!\n");
	      //node0_printf("site_prn = %e\n",s->site_prn);
#else
	      /*        random_anti_hermitian( (anti_hermitmat *)&(s->mom[dir]),
		&node_prn );  */
	      s->mom[dir].real = 0.0;
	      s->mom[dir].imag = gaussian_rand_no_u1(&node_prn);
#endif
#ifdef SCHROED_FUN
	    }
	    else{
		s->mom[dir].real = 0.0;
		s->mom[dir].imag = 0.0;
	    }
#endif
	
#ifdef DRUT_DEBUG
	    if(dir != TUP) {
	      s->mom[dir].real = 0.0;
	      s->mom[dir].imag = 0.0;
	    }
#endif
	    //node0_printf("mom[%d] = %e, ", dir, s->mom[dir].imag);
      }
      //node0_printf("\n");
    }
}

