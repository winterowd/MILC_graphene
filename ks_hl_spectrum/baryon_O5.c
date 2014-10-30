/***************** baryon_twpt.c ********************************************/

/* MIMD version 7 */
/* Heechang Na 2006/7 */

/* Calculate Two point function of baryon with two staggered light quarks, and
   one wilson heavy quark. The 2-pt for the operator O5 is;

    (2)                 -i*p*x   
   C    (p,t) = SUM    e       (-4)     epsilon     epsilon         
    mu,nu        x                             a b c       a' b' c'

                  a a'        b b'      c c'
               * G    (x,0)  G   (x,0) G    (x,0)
                  1 chi       2 chi     H mu,nu

  where a,b, and c are color indices, and mu, nu are spin indices.
  Since G_1, and G_2 are single component staggered propagators, 
  they don't have spin indices.
  G_H is wilson heavy quark propagator, which has spin and color index. */


#include "ks_hl_spectrum_includes.h"


void ks_baryon_2point_O5(field_offset ks_1, field_offset ks_2, 
			 field_offset heavy_quark, 
			 double_complex *propagator[4][4])
{
  int i, j, k, t, my_x, my_y, my_z;
  site *s;
  double_complex epx, epx1;
  double pi, mom[3];
  int eps;
  int a, b, c, ap, bp, cp; /*color indices. bp -> b prime */
  int si,sf; /*dirac indices for heavy quark propagator */
  
  su3_matrix *ks_1st_quark;
  su3_matrix *ks_2nd_quark;
  wilson_propagator *quark;
  double_complex prop;
  double_complex prop_matrix[4][4], prop_tot[4][4];
  
  int p[3] ={0,0,0};
  
  pi = 4.0 * atan( 1.);
  mom[0] = -2.*pi/(double)nx;
  mom[1] = -2.*pi/(double)ny;
  mom[2] = -2.*pi/(double)nz;
  
  FORALLSITES(i,s){
    t = s->t;
    my_x = s->x;
    my_y = s->y;
    my_z = s->z;
    
    /*initialize propagator buffer*/
    /*accumulate at propagator */
    for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
	prop_tot[si][sf].real = 0.0;
	prop_tot[si][sf].imag = 0.0;
      }
    /* asign each field */
    ks_1st_quark = (su3_matrix *)F_PT(s, ks_1);
    ks_2nd_quark = (su3_matrix *)F_PT(s, ks_2);
    quark = (wilson_propagator *)F_PT(s, heavy_quark);
    
    /*product of three propagators with sumation of colors */
    for(a=0;a<3;a++){ ap = a;
      for(j=-1;j<=1;j=j+2){
	b = bp = (a+j+3)%3;
	c = cp = (a-j+3)%3;
	eps = 1;
	
	CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
	prop.real *= eps;
	prop.imag *= eps;
	for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
	    
	    CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);
	    
	    prop_tot[si][sf].real += prop_matrix[si][sf].real;
	    prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
	  }
	
      }
    }
    
    for(a=0;a<3;a++){ ap = a;
      for(j=-1;j<=1;j=j+2){
	b = cp = (a+j+3)%3;
	c = bp = (a-j+3)%3;
	eps = -1;
	
	CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
	prop.real *= eps;
	prop.imag *= eps;
	for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
	    
	    CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);
	    
	    prop_tot[si][sf].real += prop_matrix[si][sf].real;
	    prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
	  }
	
      }
    }
    
    for(a=0;a<3;a++){ bp = a;
      for(j=-1;j<=1;j=j+2){
	b = ap = (a+j+3)%3;
	c = cp = (a-j+3)%3;
	eps = -1;
	
	CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
	prop.real *= eps;
	prop.imag *= eps;
	for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
	    
	    CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);
	    
	    prop_tot[si][sf].real += prop_matrix[si][sf].real;
	    prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
	  }
	
      }
    }
    
    for(a=0;a<3;a++){ bp = a;
      for(j=-1;j<=1;j=j+2){
	b = cp = (a+j+3)%3;
	c = ap = (a-j+3)%3;
	eps = 1;
	
	CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
	prop.real *= eps;
	prop.imag *= eps;
	for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
	    
	    CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);
	    
	    prop_tot[si][sf].real += prop_matrix[si][sf].real;
	    prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
	  }
	
      }
    }
    
    
    for(a=0;a<3;a++){ cp = a;
      for(j=-1;j<=1;j=j+2){
	b = bp = (a+j+3)%3;
	c = ap = (a-j+3)%3;
	eps = -1;
	
	CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
	prop.real *= eps;
	prop.imag *= eps;
	for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
	    
	    CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);
	    
	    prop_tot[si][sf].real += prop_matrix[si][sf].real;
	    prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
	  }
	
      }
    }
    
    for(a=0;a<3;a++){ cp = a;
      for(j=-1;j<=1;j=j+2){
	b = ap = (a+j+3)%3;
	c = bp = (a-j+3)%3;
	eps = 1;
	
	CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
	prop.real *= eps;
	prop.imag *= eps;
	for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
	    
	    CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);
	    
	    prop_tot[si][sf].real += prop_matrix[si][sf].real;
	    prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
	  }
	
      }
    }
    
    
    /* exp(i*vector_p.vector_x) */
    epx1.real = 0.0;
    epx1.imag = (mom[0]*(double)p[0]*(double)my_x + mom[1]*(double)p[1]*(double)my_y
		 + mom[2]*(double)p[2]*(double)my_z);
    epx = dcexp(&epx1);
    
    /* 4*(-1)^(...) factor */
    /* epx has information of exponential and 4*(-1)^(..) */
    
    if((1)%2 == 1) {
      epx.real *= 4.0;
      epx.imag *= 4.0;
    }
    else{
      epx.real *= -4.0;
      epx.imag *= -4.0;
    }   /*In case of O_5 */
    
    
    /*if((my_x)%2 == 1) {
      epx.real *= -4.0;
      epx.imag *= -4.0;
      }
      else{
      epx.real *= 4.0;
      epx.imag *= 4.0;
      }*/ /*In case of O_mu */
    
    
    for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
	CMUL(epx, prop_tot[si][sf], prop_tot[si][sf]);
	propagator[si][sf][t].real += prop_tot[si][sf].real;
	propagator[si][sf][t].imag += prop_tot[si][sf].imag;
      }
  }/*foralllattice loop*/
  
  return 0;
}



