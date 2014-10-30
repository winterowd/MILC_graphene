/********************* rephase_u1.c *******************************/

#include "generic_ks_includes_u1.h"

/********* phaseset() - set up KS phase vectors **********/
/* ANTIPERIODIC bc's in t and PERIODIC in x,y,z */

void phaseset() {
register site *sit; /* pointer to current site */
register int i;
    /*	phase of link(i,mu) = sum(nu<mu) { -1^i[nu] }		*/
    /*	all t phases for t=nt-1 time slice get extra minus sign	*/
    /*	   to give antiperiodic boundary conditions		*/

#ifdef PERIODICBC
    node0_printf("with periodic boundary conditions in time\n");
#endif

    FORALLMYSITES(i,sit){
	sit->phase[TUP] = 1.0;
	if( (sit->t)%2 == 1) sit->phase[XUP] = -1.0;
	  else sit->phase[XUP] = 1.0;
	if( (sit->x)%2 == 1) sit->phase[YUP] = -sit->phase[XUP];
	  else sit->phase[YUP] = sit->phase[XUP];
#ifdef FOUR_DIM
	if( (sit->y)%2 == 1) sit->phase[ZUP] = -sit->phase[YUP];
	  else sit->phase[ZUP] = sit->phase[YUP];
#else
	sit->phase[ZUP] = 1.0;
#endif


#ifndef PERIODICBC
	if( sit->t == nt-1) {
	    /* antiperiodic boundary conditions in Euclidean time */
	    sit->phase[TUP] = -sit->phase[TUP];
	}
#endif
    }
}

/************************** rephase() ******************************/
/* put Kogut-Sussind and boundary condition phase factors into or
   out of lattice */
void rephase( int flag ){
register int i,dir;
register site *s;
    /* Check to make sure we are going in expected direction */
    if( !( (flag==ON && phases_in==OFF) || (flag==OFF && phases_in==ON) ) ){
        node0_printf("DUMMY: you fouled up the phases\n");
        terminate(1);
    }
    FORALLMYSITES(i,s){
      FORALLMYUPDIR(dir) { //for(dir=XUP;dir<=TUP;dir++){
	s->link[dir].real *= s->phase[dir];
	s->link[dir].imag *= s->phase[dir];
      }
    }
    
    phases_in = flag;
} /* rephase.c */

