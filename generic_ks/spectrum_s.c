/******** spectrum_s.c *************/
/* MIMD version 6 */
/* For measuring propagation IN THE Z DIRECTION (screening spectrum) */

/* Spectrum for Kogut-Susskind pointlike hadrons, wall source */
/* Wall source is modulated by a cosine at the lowest Matsubara frequency */
/* 03/15 CW Modified for graphene code */

#include "generic_ks_includes_u1.h"

int spectrum_s(Real vmass, int src_flag, ferm_links_u1_t *fn) /* return the C.G. iteration number */
{
  Real piprop,pi2prop,rhoprop0,rhoprop1,rho2prop0,rho2prop1;
  complex ferm_prop;
  Real vmass_x2;
  register complex cc;
  Real finalrsq, th;
  register int i,x,y,t,cgn;
  int source_type;
  char *source_string[2] = {"WALL", "POINT"};

  /* Fix ZUP Coulomb gauge - gauge links only*/
  //rephase( OFF );
  //gaugefix(ZUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL);
  //rephase( ON );
  piprop=pi2prop=rhoprop0=rhoprop1=rho2prop0=rho2prop1=ferm_prop.real=ferm_prop.imag=0.;
  
  vmass_x2 = 2.*vmass;
  cgn=0;
  /* Phase increment for minimum Matsubara frequency */
  th = PI/nt;

  
  /* initialize ttt and xxx */
  clear_latvec_u1( F_OFFSET(ttt), EVENANDODD);
  clear_latvec_u1( F_OFFSET(xxx), EVENANDODD);
  if(src_flag == 1) { //construct point source at origin
    if( node_number(0,0,0,0) == mynode() )
      { 
	  i=node_index(0,0,0,0);
	  lattice[i].ttt.real = (-nx*ny*nz/8.);
	}
    source_type=1;
  }
  else { //wall source  
    for(y=0;y<ny;y+=2)for(t=0;t<nt;t+=2)
      /**for(x=0;x<1;x+=2)for(y=0;y<1;y+=2)for(t=0;t<1;t+=2)**/
      {
	if( node_number(0,y,0,t) != mynode() )continue;
	i=node_index(0,y,0,t);
	/* Modulate source with Matsubara phase */
	lattice[i].ttt.real = -cos((double)t*th);
	lattice[i].ttt.imag = sin((double)t*th);
      }
    source_type=0;
  }
  /* Multiply by -Madjoint */
  //load_ferm_links(&fn_links); //this needs to be modified (how?)
  dslash_fn_site_u1( F_OFFSET(ttt), F_OFFSET(phi), ODD, fn);
  scalar_mult_latvec_u1( F_OFFSET(ttt), -vmass_x2, F_OFFSET(phi), EVEN);
  /* do a C.G. */
  //load_ferm_links(&fn_links); //this needs to be modified (how?)
  cgn += ks_congrad_u1(F_OFFSET(phi),F_OFFSET(xxx),vmass,
		       niter, nrestart, rsqprop, PRECISION, EVENANDODD, &finalrsq,
		       fn);
      
  /* fill the hadron matrix */
  copy_latvec_u1( F_OFFSET(xxx), F_OFFSET(propmat), EVENANDODD);
    
  if(this_node==0)printf("SPAT_MES_SOURCE: %s\n",source_string[source_type]);
  /* measure the meson propagator */
  for(x=0; x<nx; x++)
    {
      /* clear meson propgators */
      piprop=rhoprop0=rhoprop1=pi2prop=rho2prop0=rho2prop1=0.;
      
      for(y=0;y<ny;y++)for(t=0;t<nt;t++)
	{
	  if( node_number(x,y,0,t) != mynode() )continue;
	  i=node_index(x,y,0,t);
	  
	  CMULJ_( lattice[i].propmat, lattice[i].propmat, cc);

	  //if(source_type==1) //project to lowest Matsubara mode for point source
	    //cc.real = cc.real*cos((double)t*th);

	  if(source_type==1) { //project to lowest Matsubara mode for point source
	    ferm_prop.real += cos((double)th*t)*lattice[i].propmat.real - sin((double)th*t)*lattice[i].propmat.imag;
	    ferm_prop.imag += cos((double)th*t)*lattice[i].propmat.imag + sin((double)th*t)*lattice[i].propmat.real;
	  }
	  else { //phase already in wall source
	    ferm_prop.real += lattice[i].propmat.real;
	    ferm_prop.imag += lattice[i].propmat.imag;
	  }

	  piprop += cc.real;
	  
	  if( (x+y)%2==0)rhoprop0 += cc.real;
	  else	   rhoprop0 -= cc.real;
	  if( (y+t)%2==0)rhoprop1 += cc.real;
	  else	   rhoprop1 -= cc.real;
	  if( (t+x)%2==0)rhoprop1 += cc.real;
	  else	   rhoprop1 -= cc.real;
	  
	  if( x%2==0)rho2prop1 += cc.real;
	  else       rho2prop1 -= cc.real;
	  if( y%2==0)rho2prop1 += cc.real;
	  else       rho2prop1 -= cc.real;
	  if( t%2==0)rho2prop0 += cc.real;
	  else       rho2prop0 -= cc.real;
	  
	  if( (x+y+t)%2==0)pi2prop += cc.real;
	  else	     pi2prop -= cc.real;
	  
	}
      
      g_sync();
      
      /* dump meson propagators */
      g_floatsum( &piprop );
      g_floatsum( &rhoprop0 );
      g_floatsum( &rhoprop1 );
      g_floatsum( &rho2prop0 );
      g_floatsum( &rho2prop1 );
      g_floatsum( &pi2prop );
      if(mynode()==0)printf("MES_SCREEN  %d  %e  %e  %e  %e  %e  %e\n",x,
			    (double)piprop,(double)rhoprop0,
			    (double)rhoprop1,
			    (double)pi2prop,(double)rho2prop0,
                            (double)rho2prop1);
      if(this_node==0)printf("FERMION_PROP_X %d %e %e\n",x, (double)ferm_prop.real, (double)ferm_prop.imag);
    } /* nx-loop */
  
  return(cgn);
} /* spectrum */
