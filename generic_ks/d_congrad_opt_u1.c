/************************* d_congrad_opt_u1.c ***********************************/
/* MIMD version 7 */
/* Added 7/29/2010, C.W.
/* Some data parallel operations for d_congrad */
/* C. DeTar 12/05  split redundant code from various d_congrads */

#include "generic_ks_includes_u1.h"
#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
//CHANGE MACROS
/* clear an complex in the lattice */
void clear_latvec_u1(field_offset v, int parity){
register int i,j;
register site *s;
register complex *vv;
    switch(parity){
	case EVEN: FOREVENSITESDOMAIN(i,s){
		vv = (complex *)F_PT(s,v);
		vv->real = vv->imag = 0.0;
            } break;
	case ODD: FORODDSITESDOMAIN(i,s){
		vv = (complex *)F_PT(s,v);
		vv->real = vv->imag = 0.0; 
            } break;
	case EVENANDODD: FORALLSITESDOMAIN(i,s){
		vv = (complex *)F_PT(s,v);
		vv->real = vv->imag = 0.0; 
            } break;
    } 
}

/* copy an complex in the lattice */
void copy_latvec_u1(field_offset src, field_offset dest, int parity){
register int i;
register site *s;
register complex *spt,*dpt;
    switch(parity){
	case EVEN: FORMYEVENSITESDOMAIN(i,s){
		s = &(lattice[i]);
		spt = (complex *)F_PT(s,src);
		dpt = (complex *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
	case ODD: FORMYODDSITESDOMAIN(i,s){
		s = &(lattice[i]);
		spt = (complex *)F_PT(s,src);
		dpt = (complex *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
	case EVENANDODD: FORALLMYSITESDOMAIN(i,s){
		s = &(lattice[i]);
		spt = (complex *)F_PT(s,src);
		dpt = (complex *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
    } 
}
/*leave the rest of routines for later as currently they are of no need CW 
  7/30/2010 */
/* scalar multiply an SU3 vector in the lattice */

void scalar_mult_latvec_u1( field_offset src, Real scalar,
			 field_offset dest, int parity)
{
  register int i;
  register site *s;
  register complex *spt,*dpt, temp_mul;
   
    switch(parity){
       case EVEN: FORMYEVENSITESDOMAIN(i,s){
      spt = (complex *)F_PT(s,src);
      dpt = (complex *)F_PT(s,dest);
      //scalar_mult_su3_vector( spt , scalar , dpt );
		CMULREAL(*spt, scalar, *dpt);
	    } break;
	case ODD: FORMYODDSITESDOMAIN(i,s){
		spt = (complex *)F_PT(s,src);
		dpt = (complex *)F_PT(s,dest);
		//scalar_mult_su3_vector( spt , scalar , dpt );
		CMULREAL(*spt, scalar, *dpt);
	    } break;
	case EVENANDODD: FORALLMYSITESDOMAIN(i,s){
		spt = (complex *)F_PT(s,src);
		dpt = (complex *)F_PT(s,dest);
		//scalar_mult_su3_vector( spt , scalar , dpt );
		CMULREAL(*spt, scalar, *dpt);
	    } break;
    } 
}

/* scalar multiply and add a complex vector in the lattice */
void scalar_mult_add_latvec_u1( field_offset src1, field_offset src2,
			     Real scalar, field_offset dest, int parity ){
  register int i;
  register site *s;
  register complex *spt1,*spt2,*dpt, temp_mul;
  

	FORMYPARITYDOMAIN(i,s,parity){
               spt1 = (complex *)F_PT(s,src1);
                spt2 = (complex *)F_PT(s,src2);
                dpt = (complex *)F_PT(s,dest);
		if(i < loopend-FETCH_UP){
		  prefetch_VVV( (complex *)F_PT((s+FETCH_UP),src1), 
				(complex *)F_PT((s+FETCH_UP),src2),
				(complex *)F_PT((s+FETCH_UP),dest));
		}
                //scalar_mult_add_su3_vector( spt1 , spt2 , scalar , dpt);
		CMULREAL( *spt2, scalar, temp_mul);
		CADD( temp_mul, *spt1, *dpt)
	} END_LOOP
}

/* d_congrad_opt_u1.c */
