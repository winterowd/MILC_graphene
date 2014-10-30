/***************** control_O5c **************************************/

/* Main procedure for heavy-light-light baryons */
/* Just one operator, namely O5 */

/* MIMD version 7 */

/* Baryon two point function using Staggered light quark propagators
 * from Fermilab and Wilson heavy quark propagators, based on
 * Ludmila's ks_hl_spectrum code. */
/* 2.25.2006 Heechang */

/* This code is modified to print every dirac components of two-point
   functions.  But, it uses only one clover quark for heavy.
   11.14.2006 Heechang. */
/* This code take care of the correct spin indices for source and
   sink.  But, it still have the -gamma_mu convention as a convention.
   12.7.2007 Heechang */

#define CONTROL
#include <string.h>
#include "ks_hl_spectrum_includes.h"

/*--------------------------------------------------------------------*/

int main(int argc,char *argv[])
{
  int prompt , k, ns, i;
  int c1,c2,s1,s2;
  double starttime,endtime,dtime;
  site *s;
  double space_vol;
  
  int status, num_prop,t,color,spin, color1, spin1;
  
  int key[4];
  
#define restrict rstrict /* C-90 T3D cludge */
  int restrict[4];
  
  double_complex *(prop[4][4]);
  //complex *prop_smear[35];
  ks_prop_file *spf1;
  ks_prop_file *spf2;
  w_prop_file *fp_in_w;
  wilson_propagator *qdest;
  wilson_propagator qtemp1;
  ks_quark_source ksqs, ksqs2;
  wilson_quark_source wqs;
  
  key[XUP] = 1;
  key[YUP] = 1;
  key[ZUP] = 1;
  key[TUP] = 0;
  
  initialize_machine(argc,argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();
  prompt = setup(); 
  setup_restrict_fourier(key, restrict);
  space_vol = (double)(nx*ny*nz);

  /* Initialize the source structures */
  init_wqs(&wqs);
  init_ksqs(&ksqs);
  init_ksqs(&ksqs2);

  while( readin(prompt) == 0){
    
    starttime=dclock();
    
    
    /**************************************************************/
    /*Allocate storage space for propagators*/
    for(spin=0;spin<4;spin++) 
      for(spin1=0;spin1<4;spin1++)
	{
	  prop[spin][spin1] = (double_complex *)malloc(nt*sizeof(double_complex));
	  for(t=0;t<nt;t++){
	    prop[spin][spin1][t].real = 0.0; 
	    prop[spin][spin1][t].imag = 0.0; 
	  }
	}
    
    /**************************************************************/
    /*load fermi-style staggered propagator*/
    /*Baryon needs two staggered propagators. -H */
    
    
    /* First staggered fermion (stag_propagator) */
    reload_ksprop_to_field(ks_prop_startflag, 
			   start_ks_prop_file, &ksqs, 
			   F_OFFSET(prop), 1);
    FORALLSITES(i,s){
      for(color = 0; color < 3; color++)
	for(k = 0; k < 3; k++)
	  s->stag_propagator.e[color][k] = s->prop[color].c[k];
    }
    
    /* Second staggered fermion (stag_propagator2) */
    reload_ksprop_to_field(ks_prop_startflag2, 
			   start_ks_prop_file2, &ksqs2, 
			   F_OFFSET(prop), 1);
    FORALLSITES(i,s){
      for(color = 0; color < 3; color++)
	for(k = 0; k < 3; k++)
	  s->stag_propagator2.e[color][k] = s->prop[color].c[k];
    }
    
    /* Load Wilson propagator for each kappa (quark_propagator) */

    for(k=0; k<num_kap; k++){
      kappa = kap[k];
      wpf = r_open_wprop(startflag_w[k], startfile_w[k]);
      for(spin=0;spin<4;spin++)
	for(color=0;color<3;color++){
	  if(reload_wprop_sc_to_field(startflag_w[k], wpf,
				      &wqs, spin, color, psi, 1) != 0)
	    terminate(1);
	  FORALLSITES(i,s){
	    copy_wvec(&psi[i],
		      &lattice[i].quark_propagator.c[color].d[spin]);
	  }
	}
      r_close_wprop(startflag_w[k],wpf);
      
      
      /*******************************************************************/
      /* Rotate the heavy quark */
      
      rotate_w_quark(F_OFFSET(quark_propagator), 
		     F_OFFSET(quark_propagator_copy), d1[k]);  
      
      // result in quark_propagator_copy
      
      /**************************************************************/
      /*Calculate and print out the spectrum */
      
      /************************************************/
      /*** two point function call ********************/
      /************************************************/
      ks_baryon_2point_O5(F_OFFSET(stag_propagator), 
			  F_OFFSET(stag_propagator2), 
			  F_OFFSET(quark_propagator), prop);
      /*ks_baryon_2point_O5() is new routine for baryon two point function 
	with two staggered quark propagators. 
	It is based on All_KS_hl_prop() in ks_hl_spectrum -Heechang */
      
      
      if(this_node==0) 
	printf("\n\nTr %d, p000_k_%f\n_________________________________\n",i, kap[k]);
      for(t=0;t<nt;t++)
	{
	  for(spin=0;spin<4;spin++) for(spin1=0;spin1<4;spin1++){
	      g_doublesum( &(prop[spin][spin1][t].real) );
	      g_doublesum( &(prop[spin][spin1][t].imag) );
	    }
	}
      g_sync();
      if(this_node==0){
	for(t=0;t<nt;t++){
	  printf("DatA %d ", t);
	  for(spin=0;spin<4;spin++) for(spin1=0;spin1<4;spin1++){
	      printf("%e %e ", prop[spin][spin1][t].real, prop[spin][spin1][t].imag);
	    }
	  printf("\n");
	}
      }
      
      for(t=0;t<nt;t++) for(spin=0;spin<4;spin++) for(spin1=0;spin1<4;spin1++){
	    prop[spin][spin1][t].real=0.0; prop[spin][spin1][t].imag=0.0; 
	  } 	
      
      
    } /*loop kappa*/
  } /*while(readin) */
}

