/*************** d_action_u1.c ****************************************/
/* MIMD version 7 */
/* Added 8/1/2010, C.W.
/* Measure total action, as needed by the hybrid Monte Carlo
   algorithm.  When this routine is called the conjugate gradient
   should already have been run on the even sites, so that in the case
   of two masses, the vectors xxx1 and xxx2 contain (M_adjoint*M)^(-1)
   * phi1 and phi2 and in the case of one mass the vector xxx contains
   (M_adjoint*M)^(-1) * phi.
*/
//#include "graphene_includes.h"
#include "ks_imp_includes.h"	/* definitions files and prototypes */
Real ahmat_mag_sq(anti_hermitmat *pt);

/*DEBUG*/
double old_g, old_h, old_f, old_a;
/*ENDDEBUG*/

double d_action_u1(){
double hmom_action(),fermion_action();
    double ssplaq,stplaq,g_action,h_action,f_action;

    d_plaquette_u1(&ssplaq,&stplaq); 
    ssplaq *= -1.0; stplaq *= -1.0;
    g_action = -beta*volume*(ssplaq+stplaq);
    node0_printf("PLAQUETTE ACTION: %e\n",g_action);

    rephase(OFF);
    g_action = (beta/3.0)*imp_gauge_action_u1(); /* leave as compact for now 10/10 CW */
					     
    rephase(ON);
    h_action = hmom_action_u1();
    f_action = fermion_action_u1();

    node0_printf("ACTION: g,h,f = %e  %e  %e  %e\n",
    g_action, h_action, f_action, g_action+h_action+f_action );

/*DEBUG*/
node0_printf("DG = %e, DH = %e, DF = %e, D = %e\n",
g_action-old_g, h_action-old_h, f_action-old_f,
g_action+h_action+f_action-old_a);
old_g=g_action; old_h=h_action; old_f=f_action;
old_a=g_action+h_action+f_action;
/*ENDDEBUG*/

    return(g_action+h_action+f_action);
}

/* fermion contribution to the action */
double fermion_action_u1() {
register int i;
register site *s;
register complex cc;
double sum;
    sum=0.0;
    FORMYEVENSITES(i,s){
	/* phi is defined on even sites only */
#ifdef ONEMASS
      CMULJ_( s->phi, s->xxx, cc );
      sum += (double)cc.real;
#else
      CMULJ_( s->phi1, s->xxx1, cc );
      sum += (double)cc.real;
      CMULJ_( s->phi2, s->xxx2, cc );
      sum += (double)cc.real;
#endif
    }
    g_doublesum( &sum );
    return(sum);
}

/* gauge momentum contribution to the action */
double hmom_action_u1() {
register int i,dir;
register site *s;
double sum;

    sum=0.0;
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){
	  sum += (double)( (s->mom[dir].imag)*(s->mom[dir].imag) );
	  //should we subtract anything here?? 08/01/2010 CW
	}
    }
    g_doublesum( &sum );
    return(sum);
}

/* copy a gauge field - an array of four complex numbers */
void gauge_field_copy_u1(field_offset src,field_offset dest){
register int i,dir,src2,dest2;
register site *s;
    FORALLSITES(i,s){
	src2=src; dest2=dest;
        for(dir=XUP;dir<=TUP; dir++){
	  CADD(*(complex *)F_PT(s,src2), 0.0, *(complex *)F_PT(s,dest2));
	  src2 += sizeof(complex);
	  dest2 += sizeof(complex);
	}
    }
#ifdef FN
    free_fn_links_u1(&fn_links_u1);
    //free_fn_links(&fn_links_dmdu0);
#endif
}
