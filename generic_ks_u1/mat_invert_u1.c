/******** mat_invert_u1.c *************/
/* MIMD version 7*/
/* DT 6/97
* separated from spectrum_mom.c 12/97
* Modify check_invert() to compute magnitude of error 11/98 DT
*/
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash" has
   been defined to be the appropriate "dslash_fn" or "dslash_eo"
*/
#include "generic_ks_includes_u1.h"
#include "../include/dslash_ks_redefine.h"

/* dst = M src. No parity selection here */

//Added "_u1" to function names (MUST CHANGE IN HEADER FILE "generic_ks_u1.h") 4/11/2010

void ks_dirac_op_u1( complex *src, complex *dst, Real mass, 
		  ferm_links_t_u1 *fn){
    register int i;
    register site *s;

    dslash_fn_field_u1( src, dst, EVENANDODD, fn);
    FORALLMYSITES(i,s){
      CMULREAL( src[i], +2.0*mass, src[i]);
      CADD( dst[i], src[i], dst[i]);
    }
}

/* dst = Madj src. No parity selection here */

void ks_dirac_adj_op_u1( complex *src, complex *dst, Real mass,
		      ferm_links_t_u1 *fn){
    register int i;
    register site *s;

    dslash_fn_field_u1( src, dst, EVENANDODD, fn);
    FORALLMYSITES(i,s){
      CMULREAL( dst[i], -1.0, dst[i]);
      CMULREAL( src[i], +2.0*mass, src[i]);
      CADD( dst[i], src[i], dst[i]);
    }
}

/* dst = Madj M src with parity selection */

void ks_dirac_opsq_u1( complex *src, complex *dst, Real mass, int parity,
		    ferm_links_t_u1 *fn){
    register int i;
    register site *s;
    int otherparity = 0;
    Real msq_x4 = 4.0*mass*mass;
    complex *tmp;
    complex temp_mul;

    tmp = (complex *)malloc(sites_on_node * sizeof(complex));
    if(tmp==NULL){
      printf("ks_dirac_op(%d): no room for tmp\n",this_node);
      terminate(1);
    }

    switch(parity){
    case(EVEN): otherparity=ODD; break;
    case(ODD):  otherparity=EVEN; break;
    case(EVENANDODD): otherparity=EVENANDODD;
    }

    dslash_field_fn_u1( src, tmp, otherparity, fn);
    dslash_field_fn_u1( tmp, dst, parity, fn);
    FORMYSITESANDPARITY(i,s,parity){
      CMULREAL( dst[i], -1.0, dst[i]);
      CMULREAL( src[i], msq_x4, temp_mul);
      CADD( dst[i], temp_mul, dst[i]);
    }

    free(tmp);
}

/* This algorithm solves even and odd sites separately */

int mat_invert_cg_field(complex *src, complex *dst, 
			 quark_invert_control *qic,
			 Real mass, ferm_links_t_u1 *fn ){
    int cgn;
    complex *tmp;

    tmp = (complex *)malloc(sites_on_node * sizeof(complex));
    if(tmp==NULL){
      printf("mat_invert_cg_field(%d): no room for tmp\n",this_node);
      terminate(1);
    }

    /* "Precondition" both even and odd sites */
    /* temp <- M_adj * src */
    ks_dirac_adj_op_u1( src, tmp, mass, fn);

    /* We don't call with EVENANDODD anymore because we are
       transitioning to the QOP/QDP standard */

    /* dst_e <- (M_adj M)^-1 temp_e  (even sites only) */
    qic->parity = EVEN;
    cgn = ks_congrad_field_u1( tmp, dst, qic, mass, fn );

    /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
    qic->parity = ODD;
    cgn = ks_congrad_field_u1( tmp, dst, qic, mass, fn );

    free(tmp);

    // check_invert_field2( dst, tmp, mass, 1e-6, fn);
    // check_invert_field( dst, src, mass, 1e-6, fn);
    return cgn;
}


/* Invert using Leo's UML trick */

/* Our M is     (  2m		D_eo   )
		( -D_eo^adj	2m )
    Note D_oe = -D_eo^adj 

   define U =	( 2m		-D_eo )
		(  0		2m    )

	  L =   ( 2m		0  )
		( D_eo^adj	2m )

   observe UML =    2m ( 4m^2+D_eo D_eo^adj	0    )
		       (  0			4m^2 )
   and the upper left corner is what our conjugate gradient(even)
   inverts, except for a funny minus sign in our cong. grad. and a
   factor of 2m

   Use M^-1 phi = L L^-1 M^-1 U^-1 U phi
		= L  (UML)^-1 U phi
		= L  -(1/2m)(congrad) U phi
*/
         

/* This algorithm solves even sites, reconstructs odd and then polishes
   to compensate for loss of significance in the reconstruction
*/

int mat_invert_uml_field_u1(complex *src, complex *dst, 
			 quark_invert_control *qic,
			 Real mass, ferm_links_t_u1 *fn ){
    int cgn;
    register int i;
    register site *s;
    complex *tmp, *ttt, *temp_abs;

    tmp = (complex *)malloc(sites_on_node * sizeof(complex));
    if(tmp==NULL){
      printf("mat_invert_uml_field(%d): no room for tmp\n",this_node);
      terminate(1);
    }
    ttt = (complex *)malloc(sites_on_node * sizeof(complex));
    if(ttt==NULL){
      printf("mat_invert_uml_field(%d): no room for ttt\n",this_node);
      terminate(1);
    }
    /* "Precondition" both even and odd sites */
    /* temp <- M_adj * src */
    dslash_fn_field_u1( src, ttt, EVENANDODD, fn);
    FORALLMYSITES(i,s){
      CMULREAL( src[i], -2.0*mass, src[i]);
      CADD( ttt[i], src[i], tmp[i]);
      CMULREAL( tmp[i], -1.0, tmp[i]);
    }

    /* dst_e <- (M_adj M)^-1 temp_e  (even sites only) */
    qic->parity     = EVEN;
    cgn = ks_congrad_field_u1( tmp, dst, qic, mass, fn ); //update instances in all files!!

    /* reconstruct odd site solution */
    /* dst_o <-  1/2m (Dslash_oe*dst_e + src_o) */
    dslash_fn_field_u1( dst, ttt, ODD, fn );
    FORMYODDSITES(i,s){
      CSUB( src[i], ttt[i], dst[i]);
      CMULREAL( dst[i], 1.0/(2.0*mass), dst[i]);
    }

    /* Polish off odd sites to correct for possible roundoff error */
    /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
    qic->parity = ODD;
    cgn = ks_congrad_field_u1( tmp, dst, qic, mass, fn );

    // check_invert_field( dst, src, mass, 1e-6, fn);
    free(tmp);
    free(ttt);

    return cgn;
}


/* FOR TESTING: multiply src by matrix and check against dest */
void check_invert_field( complex *src, complex *dest, Real mass,
			 Real tol, ferm_links_t_u1 *fn){
    register int i,k,flag;
    register site *s;
    Real r_diff, i_diff;
    double sum,sum2,dflag,dmaxerr,derr;
    complex *tmp;

    tmp = (complex *)malloc(sites_on_node * sizeof(complex));
    if(tmp==NULL){
      printf("check_invert_field(%d): no room for tmp\n",this_node);
      terminate(1);
    }

    /* Compute tmp = M src */
    ks_dirac_op_u1( src, tmp, mass, fn);

    sum2=sum=0.0;
    dmaxerr=0;
    flag = 0;
    FORALLMYSITES(i,s){
	for(k=0;k<3;k++){
	    r_diff = dest[i].c[k].real - tmp[i].c[k].real;
	    i_diff = dest[i].c[k].imag - tmp[i].c[k].imag;
	    if( fabs(r_diff) > tol || fabs(i_diff) > tol ){
	      printf("site %d color %d  expected ( %.4e , %.4e ) got ( %.4e , %.4e )\n",
		     i,k,
		     dest[i].c[k].real, dest[i].c[k].imag,
		     tmp[i].c[k].real, tmp[i].c[k].imag);
	      flag++;
	    }
	    derr = r_diff*r_diff + i_diff*i_diff;
	    if(derr>dmaxerr)dmaxerr=derr;
 	    sum += derr;
	}
	//not sure about this format
	sum2 += (double)cabs_sq(dest+i);
    }
    g_doublesum( &sum );
    g_doublesum( &sum2 );
    dflag=flag;
    g_doublesum( &dflag );
    g_doublemax( &dmaxerr );
    if(this_node==0){
      printf("Inversion checked, frac. error = %e\n",sqrt(sum/sum2));
      printf("Flagged comparisons = %d\n",(int)dflag);
      printf("Max err. = %e frac. = %e\n",sqrt(dmaxerr),
	     sqrt(dmaxerr*volume/sum2));
      fflush(stdout);
    }
    free(tmp);
}

/* FOR TESTING: multiply src by Madj M and check against dest */
void check_invert_field2_u1( complex *src, complex *dest, Real mass,
			  Real tol, ferm_links_t_u1 *fn){
    register int i,k,flag;
    register site *s;
    Real r_diff, i_diff;
    double sum,sum2,dflag,dmaxerr,derr;
    complex *tmp;

    tmp = (complex *)malloc(sites_on_node * sizeof(complex));
    if(tmp==NULL){
      printf("check_invert_field(%d): no room for tmp\n",this_node);
      terminate(1);
    }

    /* Compute tmp = (Madj M) src */
    ks_dirac_opsq_u1( src, tmp, mass, EVENANDODD, fn);

    sum2=sum=0.0;
    dmaxerr=0;
    flag = 0;
    FORALLMYSITES(i,s){
	for(k=0;k<3;k++){
	    r_diff = dest[i].c[k].real - tmp[i].c[k].real;
	    i_diff = dest[i].c[k].imag - tmp[i].c[k].imag;
	    if( fabs(r_diff) > tol || fabs(i_diff) > tol ){
	      printf("site %d color %d  expected ( %.4e , %.4e ) got ( %.4e , %.4e )\n",
		     i,k,
		     dest[i].c[k].real, dest[i].c[k].imag,
		     tmp[i].c[k].real, tmp[i].c[k].imag);
	      flag++;
	    }
	    derr = r_diff*r_diff + i_diff*i_diff;
	    if(derr>dmaxerr)dmaxerr=derr;
 	    sum += derr;
	}
	//not sure about this format
	sum2 += (double)cabs_sq(dest+i);
    }
    g_doublesum( &sum );
    g_doublesum( &sum2 );
    dflag=flag;
    g_doublesum( &dflag );
    g_doublemax( &dmaxerr );
    if(this_node==0){
      printf("Inversion checked, frac. error = %e\n",sqrt(sum/sum2));
      printf("Flagged comparisons = %d\n",(int)dflag);
      printf("Max err. = %e frac. = %e\n",sqrt(dmaxerr),
	     sqrt(dmaxerr*volume/sum2));
      fflush(stdout);
    }
    free(tmp);
}

