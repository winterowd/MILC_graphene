/******** spectrum2.c *************/
/* MIMD version 7 */
/* Spectrum for Kogut-Susskind pointlike hadrons, wall source 
* MIMD version 7 with SG additions
* updated 2/21/95, by SG.
* modified 7/96 DT, do both Naik and conventional spectrum, wall
*	sources only. 
* 7/29/96 DT fixed names of rho and rho2 propagators.
*   rho is VT (rho-b_1), rho2 is PV (rho-a_1)
*   also pi and pi2 -> pi_ps_prop, pi_sc_prop
* modified 8/96 DT changed source and source_inc
* modified 5/97 DT for improved actions
* modified 2/99 DT for new output format
*
* This version does arbitrary number of wall sources 
* This version DOES NOT fix the gauge -- you should do that before 
*	calling it
* 5/9/06 CD Restored point source capability (zero momentum only)
*     source_start >= 0 with
*      source_start + n_sources*source_inc < nt invokes corner wall src
*     source_start >= nt with
*      source_start + n_sources*source_inc < 2*nt invokes point src
*     Caution: any other choice mixes point and corner wall sources.
*/
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash_site" has
   been defined to be the appropriate "dslash_fn_site" or "dslash_eo_site"
*/

#include "generic_ks_includes_u1.h"

int spectrum2_u1( Real vmass, field_offset temp1, field_offset temp2,
	       ferm_links_u1_t *fn){ 
  /* return the C.G. iteration number */
  double *pi_ps_prop,*pi_ps_prop_temp,*pi_sc_prop,*rho_pv_prop,*rho_vt_prop;
  complex *fermion_prop_x, *scal_prop_x, *scal_prop_t;
  complex *fermion_prop_ta, *fermion_prop_tb, *fermion_prop_tc, *fermion_prop_td;
  Real vmass_x2, pi;
  site *s;
  register complex cc;
  Real finalrsq;
  int isrc;
  register int i,x,y,z,t,cgn;
  register int t_source,t_off;
  int source_type = 0;  /* 1 corner wall 2 point 3 mixed */
  char *source_string[4] = {"GOOFED","CORNER","POINT","MIXED"};

  vmass_x2 = 2.*vmass;
  pi = M_PI;
  z=0; //fermions on z=0 plane
  cgn=0;
  scal_prop_x = (complex *)malloc(nx*sizeof(complex));
  scal_prop_t = (complex *)malloc(nt*sizeof(complex));
  fermion_prop_x = (complex *)malloc(nx*sizeof(complex));
  fermion_prop_ta = (complex *)malloc(nt*sizeof(complex));
  fermion_prop_tb = (complex *)malloc(nt*sizeof(complex));
  fermion_prop_tc = (complex *)malloc(nt*sizeof(complex));
  fermion_prop_td = (complex *)malloc(nt*sizeof(complex));
  pi_ps_prop = (double *)malloc(nt*sizeof(double) );   /* "pi" */
  pi_ps_prop_temp = (double *)malloc(nt*sizeof(double) );
  pi_sc_prop = (double *)malloc(nt*sizeof(double) );   /* "pi2" */
  rho_vt_prop = (double *)malloc(nt*sizeof(double) );  /* "rho" */
  rho_pv_prop = (double *)malloc(nt*sizeof(double) );  /* "rho2" */
  for( t=0; t<nt; t++){
    pi_ps_prop[t] = pi_ps_prop_temp[t] = pi_sc_prop[t] = rho_pv_prop[t] = rho_vt_prop[t] = scal_prop_t[t].real = scal_prop_t[t].imag = 0.0;
    fermion_prop_ta[t].real = fermion_prop_ta[t].imag = fermion_prop_tb[t].real = fermion_prop_tb[t].imag = 0.;
    fermion_prop_tc[t].real = fermion_prop_tc[t].imag = fermion_prop_td[t].real = fermion_prop_td[t].imag = 0.;
  }
  for( x=0; x<nx; x++) {
    fermion_prop_x[x].real = fermion_prop_x[x].imag = scal_prop_x[x].real = scal_prop_x[x].imag = 0.0;
  }
  

  for(t_source=source_start, isrc=0; t_source<2*nt && isrc < n_sources; 
	++isrc, t_source += source_inc ) {
    if(this_node==0)printf("spectrum(): source time = %d\n",t_source);
    
        /* initialize temp1 and temp2 */
    clear_latvec_u1( temp1, EVENANDODD);
    clear_latvec_u1( temp2, EVENANDODD); 
    clear_latvec_u1( F_OFFSET(ttt), EVENANDODD); //debugging
    clear_latvec_u1( F_OFFSET(propmat), EVENANDODD);

    if(t_source < nt){  /* wall source */
      for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
		if( node_number(x,y,z,t_source) != mynode() )continue;
		i=node_index(x,y,0,t_source);
		((complex *)(F_PT(&lattice[i],temp1)))->real
		  = -1.0;
	      }
	  source_type |= 1;
    }
    else{  /* point source at origin */
      if( node_number(0,0,0,t_source%nt) == mynode() )
	{ 
	  i=node_index(0,0,0,t_source%nt);
	  ((complex *)(F_PT(&lattice[i],temp1)))->real
	    = (-nx*ny*nz/8.);
	}
      source_type |= 2;
    }
    
    //debugging test for dslash_field and dslash_site
    /*complex *t_src, *t_dest;
    //msg_tag *tags1[16], *tags2[16];
    t_src  = (complex *)malloc(sizeof(complex)*sites_on_node);
    t_dest = (complex *)malloc(sizeof(complex)*sites_on_node);
    int dir;
    FORALLMYUPDIR(dir) {
    FORALLMYSITES(i,s){
      t_src[i]  = *((complex *)F_PT(s,temp1) );
      t_dest[i].real = t_dest[i].imag = 0.0;
      node0_printf("site %d dir; %d lng: %e %e\n", i, dir, fn->lng[4*i+dir].real,
		   fn->lng[4*i+dir].imag);
    }
    }
    dslash_fn_field_u1( t_src, t_dest, ODD, fn);
    dslash_fn_field_u1( t_dest, t_dest, EVEN, fn);
    //dslash_fn_site_u1( temp1, F_OFFSET(ttt), ODD, fn);
    //scalar_mult_latvec_u1( temp1, vmass_x2, F_OFFSET(ttt), EVEN);
    node0_printf("DSLASH_SITE TEST!\n");
    FORALLMYSITES(i,s) {
      node0_printf("site %d %d %d %d phi: %e %e ttt: %e %e\n", s->x, s->y, 
		   s->z, s->t, s->phi.real, s->phi.imag, s->ttt.real, s->ttt.imag);
		   }
    node0_printf("DSLASH_FIELD TEST!\n");
    FORALLMYSITES(i,s) {
      node0_printf("site %d %d %d %d dest: %e %e\n", s->x, s->y, 
		   s->z, s->t, t_dest[i].real, t_dest[i].imag);
    }
    free(t_src); free(t_dest);*/

    /* do a C.G. (source in temp1, result in temp2) */
    if(t_source%2 == 0) {
      cgn += ks_congrad_u1( temp1, temp2, vmass,
				   niter, nrestart, rsqprop, PRECISION, 
			 EVEN, &finalrsq, fn);
      /* Multiply by -Madjoint */
      dslash_fn_site_u1( temp2, F_OFFSET(ttt), ODD, fn);
      scalar_mult_latvec_u1( temp2, -vmass_x2, F_OFFSET(ttt),
			EVEN);
      /**copy_latvec( temp1, F_OFFSET(g_rand), EVENANDODD );
	 checkmul();**/
      /**check_invert( F_OFFSET(ttt), temp1 );**/
    }
    else {
      cgn += ks_congrad_u1( temp1, temp2, vmass,
			 niter, nrestart, rsqprop, PRECISION, 
			 ODD, &finalrsq, fn);
      /* Multiply by -Madjoint */
      dslash_fn_site_u1( temp2, F_OFFSET(ttt), EVEN, fn);
      scalar_mult_latvec_u1( temp2, -vmass_x2, F_OFFSET(ttt),
			  ODD);
    }

    /* fill the hadron matrix */
    copy_latvec_u1( F_OFFSET(ttt), F_OFFSET(propmat), EVENANDODD);
    // } /* end loop on icol */
    
    /*node0_printf("KS_CONGRAD RESULT:\n");
    FORALLMYSITES(i,s) {
      node0_printf("site %d %d %d %d xxx: %e %e ttt: %e %e\n", s->x, s->y, 
		   s->z, s->t, s->xxx.real, s->xxx.imag, s->ttt.real, s->ttt.imag);
		   }*/

    //multiply by M to get back source vector???? (debugging test) 1/15/14
    /*complex temp_mul;
    t_src  = (complex *)malloc(sizeof(complex)*sites_on_node);
    t_dest = (complex *)malloc(sizeof(complex)*sites_on_node);
    FORALLMYSITES(i,s){
      t_src[i]  = *((complex *)F_PT(s,F_OFFSET(ttt)) );
      t_dest[i].real = t_dest[i].imag = 0.0;
    }
    dslash_fn_field_u1( t_src, t_dest, EVENANDODD, fn);
    FORALLMYSITES(i,s){
      CMULREAL(t_src[i], vmass_x2, temp_mul);
      CADD(t_dest[i], temp_mul, t_dest[i]);
    }
    //dslash_fn_site_u1( F_OFFSET(ttt), temp2, EVENANDODD, fn);
    scalar_mult_add_latvec_u1( temp2, F_OFFSET(ttt), vmass_x2, temp2,
      EVENANDODD);
    node0_printf("CHECK INVERSE RESULT:\n");
    FORALLMYSITES(i,s) {
      node0_printf("site %d %d %d %d xxx: %e %e\n", s->x, s->y, 
		   s->z, s->t, t_dest[i].real, t_dest[i].imag);
      if( i != 0 && t_dest[i].real > .00000001)
	node0_printf("ALERT! NON_ZERO ELEMENT!\n");
      node0_printf("site %d %d %d %d xxx: %e %e\n", s->x, s->y, 
	s->z, s->t, s->xxx.real, s->xxx.imag);
    }    
    free(t_src); free(t_dest);*/

    /* measure the meson propagator */
    for(t=0; t<nt; t++){
      /* define the time value offset t from t_source */
      t_off = (t+t_source)%nt;
	    
      for(x=0;x<nx;x++)for(y=0;y<ny;y++) { //for(z=0;z<nz;z++){
	if( node_number(x,y,0,t_off) != mynode() )continue;
	i=node_index(x,y,0,t_off);
	CMULJ_( lattice[i].propmat, lattice[i].propmat, cc);
	/*cc = su3_dot( &lattice[i].propmat[icol],
	  &lattice[i].propmat[icol] );*/
	
	scal_prop_x[x].real += lattice[i].xxx.real;
	scal_prop_x[x].imag += lattice[i].xxx.imag;
	
	scal_prop_t[t].real += lattice[i].xxx.real;
	scal_prop_t[t].imag += lattice[i].xxx.imag;
	
	//collect sums of all four sites on a square
	if( x%2 == 0 && y%2 == 0) { 
	  fermion_prop_ta[t].real += lattice[i].propmat.real;
	  fermion_prop_ta[t].imag += lattice[i].propmat.imag;
	}
	if( x%2 == 0 && y%2 == 1) { 
	  fermion_prop_tb[t].real += lattice[i].propmat.real;
	  fermion_prop_tb[t].imag += lattice[i].propmat.imag;
	}
	if( x%2 == 1 && y%2 == 1) { 
	  fermion_prop_tc[t].real += lattice[i].propmat.real;
	  fermion_prop_tc[t].imag += lattice[i].propmat.imag;
	}
	if( x%2 == 1 && y%2 == 0) { 
	  fermion_prop_td[t].real += lattice[i].propmat.real;
	  fermion_prop_td[t].imag += lattice[i].propmat.imag;
	}
	
	fermion_prop_x[x].real += cos(pi*t/nt)*lattice[i].propmat.real - sin(pi*t/nt)*lattice[i].propmat.imag;
	fermion_prop_x[x].imag += cos(pi*t/nt)*lattice[i].propmat.imag + sin(pi*t/nt)*lattice[i].propmat.real;
	
	pi_ps_prop[t] += cc.real;
	pi_ps_prop_temp[t] += cc.real;
	
	if( (x+y)%2==0)rho_pv_prop[t] += cc.real;
	else	   rho_pv_prop[t] -= cc.real;
	if( (y+z)%2==0)rho_pv_prop[t] += cc.real;
	else	   rho_pv_prop[t] -= cc.real;
	if( (z+x)%2==0)rho_pv_prop[t] += cc.real;
	else	   rho_pv_prop[t] -= cc.real;
	
	if( x%2==0)rho_vt_prop[t] += cc.real;
	else       rho_vt_prop[t] -= cc.real;
	if( y%2==0)rho_vt_prop[t] += cc.real;
	else       rho_vt_prop[t] -= cc.real;
	if( z%2==0)rho_vt_prop[t] += cc.real;
	else       rho_vt_prop[t] -= cc.real;
	
	if( (x+y+z)%2==0)pi_sc_prop[t] += cc.real;
	else	     pi_sc_prop[t] -= cc.real;
	
      } /* x,y,z */
      
    } /* nt-loop */
    for(t=0; t<nt; t++) {
      printf("DEBUG: %d %d %e\n", isrc, t, pi_ps_prop_temp[t]);
      pi_ps_prop_temp[t]=0; 
    }
    
  } /* end loop on t_source */
  
  /* dump the propagators */
  g_veccomplexsum( fermion_prop_x, nx);
  g_veccomplexsum( fermion_prop_ta, nt);
  g_veccomplexsum( fermion_prop_tb, nt);
  g_veccomplexsum( fermion_prop_tc, nt);
  g_veccomplexsum( fermion_prop_td, nt);
  g_veccomplexsum( scal_prop_x, nx);
  g_veccomplexsum( scal_prop_t, nt);
  g_vecdoublesum( pi_ps_prop, nt );
  g_vecdoublesum( rho_pv_prop, nt );
  g_vecdoublesum( rho_vt_prop, nt );
  g_vecdoublesum( pi_sc_prop, nt );
  if( this_node==0 ){
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",vmass,vmass);
    printf("SOURCE: %s\n",source_string[source_type]);
    printf("SINKS: PION_PS PION_SC RHO_VT RHO_PV\n");
    for(t=0;t<nt;t++) printf("%d %e 0.0 %e 0.0 %e 0.0 %e 0.0\n",t,
	   pi_ps_prop[t]/n_sources, pi_sc_prop[t]/n_sources,
	   rho_vt_prop[t]/n_sources, rho_pv_prop[t]/n_sources);
    for(x=0; x<nx; x++) { 
      printf("SPAT_FERM_PROP: %d %e %e\n", x, 
	     fermion_prop_x[x].real/n_sources, fermion_prop_x[x].imag/n_sources);
    }
    for(t=0;t<nt;t++) { 
      printf("TEMP_FERM_PROP (A): %d %e %e\n", t, fermion_prop_ta[t].real/n_sources, fermion_prop_ta[t].imag/n_sources);
      printf("TEMP_FERM_PROP (B): %d %e %e\n", t, fermion_prop_tb[t].real/n_sources, fermion_prop_tb[t].imag/n_sources);
      printf("TEMP_FERM_PROP (C): %d %e %e\n", t, fermion_prop_tc[t].real/n_sources, fermion_prop_tc[t].imag/n_sources);
      printf("TEMP_FERM_PROP (D): %d %e %e\n", t, fermion_prop_td[t].real/n_sources, fermion_prop_td[t].imag/n_sources);
    }

    for(x=0; x<nx; x++) { 
      printf("SCALS: %d %e %e\n", x, 
	     scal_prop_x[x].real/n_sources, scal_prop_x[x].imag/n_sources);
    }
    for(t=0;t<nt;t++) printf("SCALT: %d %e %e\n", t, 
	      scal_prop_t[t].real/n_sources, scal_prop_t[t].imag/n_sources);
    

    //printf("SPAT_ENDPROP\n");

    fflush(stdout);
  }
  free(fermion_prop_x); free(fermion_prop_ta); free(fermion_prop_tb);
  free(fermion_prop_tc); free(fermion_prop_td); free( pi_ps_prop ); 
  free( pi_sc_prop ); free( rho_pv_prop ); free( rho_vt_prop );
  free( pi_ps_prop_temp );
  return(cgn);
} /* spectrum */

