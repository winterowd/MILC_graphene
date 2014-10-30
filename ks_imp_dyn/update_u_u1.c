/*************** update_u_u1.c ************************************/
/* MIMD version 7 (Modified 7/28/2010)*/

/* update the link matrices					*
*  								*
*  Go to sixth order in the exponential of the momentum		*
*  matrices, since unitarity seems important.			*
*  Evaluation is done as:					*
*	exp(H) * U = ( 1 + H + H^2/2 + H^3/3 ...)*U		*
*	= U + H*( U + (H/2)*( U + (H/3)*( ... )))		*
*								*
*/

#include "ks_imp_includes.h"	/* definitions files and prototypes */

void update_u_u1( Real eps ){

register int i,dir;
register site *s;
 complex *link,temp1,temp2,htemp,temp_mul;
register Real t2,t3,t4,t5,t6;
/**TEMP**
Real gf_x,gf_av,gf_max;
int gf_i,gf_j;
**END TEMP **/

/**double dtime,dtime2,dclock();**/
/**dtime = -dclock();**/

/* Temporary by-hand optimization until pgcc compiler bug is fixed */
t2 = eps/2.0;
t3 = eps/3.0;
t4 = eps/4.0;
t5 = eps/5.0;
t6 = eps/6.0;

/** TEMP **
gf_av=gf_max=0.0;
**END TEMP**/
#ifdef FN
    free_fn_links_u1(&fn_links);
    free_fn_links_u1(&fn_links_dmdu0);
#endif
    FORALLSITES(i,s){
	for(dir=XUP; dir <=TUP; dir++){
	  //uncompress_anti_hermitian( &(s->mom[dir]) , &htemp );
	  htemp.real = 0.0;
	  htemp.imag = s->mom[dir].imag;
	  link = &(s->link[dir]);
	  CMUL( htemp, *link, temp1);
	  
	  CMULREAL(temp1, t6, temp_mul);
	  CADD(*link, temp_mul, temp2);
	  CMUL(htemp, temp2, temp1);
	  
	  CMULREAL(temp1, t5, temp_mul);
	  CADD(temp_mul, *link, temp2);
	  CMUL(htemp, temp2, temp1);
	  
	  CMULREAL(temp1, t4, temp_mul);
	  CADD(temp_mul, *link, temp2);
	  CMUL(htemp, temp2, temp1);
	  	  
	  CMULREAL(temp1, t3, temp_mul);
	  CADD(temp_mul, *link, temp2);
	  CMUL(htemp, temp2, temp1);
	  	  
	  CMULREAL(temp1, t2, temp_mul);
	  CADD(temp_mul, *link, temp2);
	  CMUL(htemp, temp2, temp1);
	  
	  CMULREAL(temp1, eps, temp_mul);
	  CADD(temp_mul, *link, temp2);
	  link->real = temp2.real;
	  link->imag = temp2.imag;
	}
    }

/**dtime += dclock();
node0_printf("LINK_UPDATE: time = %e  mflops = %e\n",
dtime, (double)(5616.0*volume/(1.0e6*dtime*numnodes())) );**/
} /* update_u */
