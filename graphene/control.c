/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for U(1) with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This version combines code for the PHI algorithm (approriate for 4
   flavors) and the R algorithm for "epsilon squared" updating of 
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */

#define CONTROL
#include "graphene_includes.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#endif

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int
main( int argc, char **argv )
{
  //test change for CVS
#ifdef DEBUG_DIAG
  complex *src, *dest;
  complex **t_fl = &fn_links.fat;
  complex **t_ll = &fn_links.lng;  
  register complex *fat1;
  register complex *long1;
  Real final_rsq;
  int iters;
#endif
  int meascount,traj_done;
  int prompt;
  int s_iters, avs_iters, avspect_iters, avbcorr_iters;
  complex *temp_field, *temp_link;
#ifdef NON_COMPACT
  Real *temp_pot;
  Real *temp;
#endif
  int i, dir;
  site *s;
  double dtime, dclock();
  
  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();

//  restore_random_state_scidac_to_site("randsave", F_OFFSET(site_prn));
//  restore_color_vector_scidac_to_site("xxx1save", F_OFFSET(xxx1),1);
//  restore_color_vector_scidac_to_site("xxx2save", F_OFFSET(xxx2),1);

#ifdef DEBUG_DIAG
  /*Debug site routines*/

  /* Allocate space for t_ll if NULL */
  if(*t_ll == NULL){
    *t_ll = (complex *)malloc(sites_on_node*4*sizeof(complex));
    if(*t_ll==NULL){
      printf("load_longlinks(%d): no room for t_ll\n", this_node);
      terminate(1);
    }
  }
 
  /* Allocate space for t_fl if NULL */
  if(*t_fl == NULL){
    *t_fl = (complex *)malloc(sites_on_node*4*sizeof(complex));
    if(*t_fl==NULL){
      printf("NODE %d: no room for t_fl\n",this_node);
      terminate(1);
    }
  }

  if(src == NULL){
    src = (complex *)malloc(sites_on_node*sizeof(complex));
    if(src == NULL){
      printf("NODE %d: no room for src\n",this_node);
      terminate(1);
    }
  }

  if(dest == NULL){
    dest = (complex *)malloc(sites_on_node*sizeof(complex));
    if(dest == NULL){
      printf("NODE %d: no room for dest\n",this_node);
      terminate(1);
    }
  }


  FORALLUPDIR(dir) { //set up fat and long links
    FORALLSITES(i,s){
      fat1 = (*t_fl) + 4*i + dir;
      long1 = (*t_ll) + 4*i + dir;
      fat1->real = dir; 
      fat1->imag = 2.0*dir;
      long1->real = 100*s->y; 
      long1->imag = 100*s->t;
    }
  }

  FORALLSITES(i, s) { //src is phi
    s->phi.real = (float) (i%10);
    s->phi.imag = (float) (i%10);
    //src[i].real = 10.0*s->x;
    //src[i].imag = 25.0*s->y;
  }

  /*dslash_fn_site_u1( F_OFFSET(phi), F_OFFSET(xxx), EVENANDODD, 
    &fn_links );*/
  /*if( readin(prompt) == 0 ) {
    printf("max_cg_iterations = %d\n", niter);
    fn_links.valid = 1;
    dslash_fn_field_u1(&src[0], &dest[0], ODD, &fn_links);
    //dslash_fn_field_u1(&dest[0], &dest[0], EVEN, &fn_links);
    dslash_fn_field_u1(&src[0], &dest[0], ODD, &fn_links);
    //dslash_fn_field_u1(&dest[0], &dest[0], EVEN, &fn_links);
    FORALLSITES(i,s) //print out result
      printf("dest[%d] = %f, %f\n", i, dest[i].real, 
	     dest[i].imag);  
    iters = ks_congrad_u1( F_OFFSET(phi), F_OFFSET(xxx), mass, 
			   niter, nrestart, rsqmin, PRECISION, 
			   EVEN, &final_rsq, &fn_links ); 
    
    FORALLSITES(i,s) //print out result
      printf("dest[%d] = %f, %f\n", i, lattice[i].xxx.real, 
	     lattice[i].xxx.imag);  
  }
  */

  //free(*t_fl);
  //free(*t_ll);
  //free(src);
  //free(dest);

#endif

  /* loop over input sets */
  while( readin(prompt) == 0) {
    
    /* perform warmup trajectories */
    dtime = -dclock();
    for( traj_done=0; traj_done < warms; traj_done++ ){
      update_u1();
    }
    node0_printf("WARMUPS COMPLETED\n"); fflush(stdout);
    
    /* perform measuring trajectories, reunitarizing and measuring 	*/
    meascount=num_accepts=num_rejects=0;  /* number of measurements 		*/
    avspect_iters = avs_iters = avbcorr_iters = 0;
    for( traj_done=0; traj_done < trajecs; traj_done++ ){ 
      
      /* do the trajectories */
      s_iters=update_u1();
            
      /* measure every "propinterval" trajectories */
      if( (traj_done%propinterval)==(propinterval-1) ){
	
	/* call gauge_variable fermion_variable measuring routines */
	/* results are printed in output file */
	rephase(OFF);
	g_measure_u1( );
	rephase(ON);

	/* Load fat and long links for fermion measurements */
	load_ferm_links_u1(&fn_links, &ks_act_paths);
#ifdef DM_DU0
	load_ferm_links_u1(&fn_links_dmdu0, &ks_act_paths_dmdu0);
#endif

	/* Measure pbp, etc */
#ifdef ONEMASS
	f_meas_imp_u1(F_OFFSET(phi),F_OFFSET(xxx),mass, &fn_links, 
		   &fn_links_dmdu0);
#else
	f_meas_imp_u1( F_OFFSET(phi1), F_OFFSET(xxx1), mass1, 
		    &fn_links, &fn_links_dmdu0);
	f_meas_imp_u1( F_OFFSET(phi2), F_OFFSET(xxx2), mass2,
		    &fn_links, &fn_links_dmdu0);
#endif

	/* Measure derivatives wrto chemical potential */
#ifdef D_CHEM_POT
#ifdef ONEMASS
	Deriv_O6( F_OFFSET(phi1), F_OFFSET(xxx1), F_OFFSET(xxx2), mass,
		  &fn_links, &fn_links_dmdu0);
#else
	Deriv_O6( F_OFFSET(phi1), F_OFFSET(xxx1), F_OFFSET(xxx2), mass1,
		  &fn_links, &fn_links_dmdu0);
	Deriv_O6( F_OFFSET(phi1), F_OFFSET(xxx1), F_OFFSET(xxx2), mass2,
		  &fn_links, &fn_links_dmdu0);
#endif
#endif

#ifdef SPECTRUM 
	//Fix gauge by requiring sum of scalar potential over entire 
	//lattice to be zero
	double avg_potential=0.;
	int t;
	//debugging free propagator so comment out for now
	for(t=0; t<nt; t++) {
	  FORALLSITES(i, s) {
	    if(s->t == t)
	      avg_potential += s->potential[TUP];
	  }
	  avg_potential = (double) (avg_potential/((double)nx*ny*nz));
	  FORALLSITES(i, s) {
	    if(s->t == t)
	      s->potential[TUP] = s->potential[TUP] - avg_potential;
	  }
	  avg_potential=0.;
	}
	/* Fix Landau gauge - gauge links only*/
	rephase( OFF );
	//NOT NECESSARY FOR TIME-LINKS ONLY MEASURING NEUTRAL MESONS
	//gaugefix_u1(8,(Real)1.8,2500,(Real)GAUGE_FIX_TOL);
	rephase( ON );
#ifdef FN
	invalidate_all_ferm_links_u1(&fn_links);
#ifdef DM_DU0
	invalidate_all_ferm_links_u1(&fn_links_dmdu0);
#endif
#endif
	/* Load fat and long links for fermion measurements */
	load_ferm_links_u1(&fn_links, &ks_act_paths);   
#ifdef DM_DU0
	load_ferm_links_u1(&fn_links_dmdu0, &ks_act_paths_dmdu0);
#endif	
	if(strstr(spectrum_request,",spectrum,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum2_u1(mass,F_OFFSET(phi),F_OFFSET(xxx),
					&fn_links);
#else
	  avspect_iters += spectrum2( mass1, F_OFFSET(phi1),
				      F_OFFSET(xxx1), &fn_links);
	  avspect_iters += spectrum2( mass2, F_OFFSET(phi1),
				      F_OFFSET(xxx1), &fn_links);
#endif
	}
      /*	
	if(strstr(spectrum_request,",spectrum_point,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_fzw(mass,F_OFFSET(phi),F_OFFSET(xxx),
					&fn_links);
#else
	  avspect_iters += spectrum_fzw( mass1, F_OFFSET(phi1),
					 F_OFFSET(xxx1), &fn_links);
	  avspect_iters += spectrum_fzw( mass2, F_OFFSET(phi1),
					 F_OFFSET(xxx1), &fn_links);
#endif
	}
	
	if(strstr(spectrum_request,",nl_spectrum,") != NULL){
#ifdef ONEMASS
	  avspect_iters += nl_spectrum(mass,F_OFFSET(phi),F_OFFSET(xxx),
				       F_OFFSET(tempmat1),F_OFFSET(staple),
				       &fn_links);
#else
	  avspect_iters += nl_spectrum( mass1, F_OFFSET(phi1), 
		F_OFFSET(xxx1), F_OFFSET(tempmat1),F_OFFSET(staple),
					&fn_links);
#endif
	}
	
	if(strstr(spectrum_request,",spectrum_mom,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_mom(mass,mass,F_OFFSET(phi),5e-3,
					&fn_links);
#else
	  avspect_iters += spectrum_mom( mass1, mass1, 
					 F_OFFSET(phi1), 1e-1,
					 &fn_links);
#endif
	}
	
	if(strstr(spectrum_request,",spectrum_multimom,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_multimom(mass,
					     spectrum_multimom_low_mass,
					     spectrum_multimom_mass_step,
					     spectrum_multimom_nmasses,
					     5e-3, &fn_links);
#else
	  avspect_iters += spectrum_multimom(mass1,
					     spectrum_multimom_low_mass,
					     spectrum_multimom_mass_step,
					     spectrum_multimom_nmasses,
					     5e-3, &fn_links);

#endif
	}
	
#ifndef ONEMASS
	if(strstr(spectrum_request,",spectrum_nd,") != NULL){
	  avspect_iters += spectrum_nd( mass1, mass2, 1e-1,
					&fn_links);
	}
#endif
	if(strstr(spectrum_request,",spectrum_nlpi2,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_nlpi2(mass,mass,F_OFFSET(phi),5e-3,
					  &fn_links );
#else
	  avspect_iters += spectrum_nlpi2( mass1, mass1, 
					   F_OFFSET(phi1),1e-1,
					   &fn_links );
	  avspect_iters += spectrum_nlpi2( mass2, mass2, 
					   F_OFFSET(phi1),1e-1,
					   &fn_links );
#endif
	}
	
	if(strstr(spectrum_request,",spectrum_singlets,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_singlets(mass, 5e-3, F_OFFSET(phi),
					     &fn_links);
#else
	  avspect_iters += spectrum_singlets(mass1, 5e-3, F_OFFSET(phi1),
					     &fn_links );
	  avspect_iters += spectrum_singlets(mass2, 5e-3, F_OFFSET(phi1),
					     &fn_links );
#endif
	}

	if(strstr(spectrum_request,",fpi,") != NULL)
	  {
	    avspect_iters += fpi_2( fpi_mass, fpi_nmasses, 2e-3,
				    &fn_links );
	  }
	
#ifdef HYBRIDS
	if(strstr(spectrum_request,",spectrum_hybrids,") != NULL){
#ifdef ONEMASS
	  avspect_iters += spectrum_hybrids( mass,F_OFFSET(phi),1e-1,
					     &fn_links);
#else
	  avspect_iters += spectrum_hybrids( mass1, F_OFFSET(phi1), 5e-3,
					     &fn_links);
	  avspect_iters += spectrum_hybrids( mass2, F_OFFSET(phi1), 2e-3,
					     &fn_links);
#endif
	}
#endif
	if(strstr(spectrum_request,",hvy_pot,") != NULL){
	  rephase( OFF );
	  hvy_pot( F_OFFSET(link[XUP]) );
	  rephase( ON );
	}
*/
#endif /* SPECTRUM */
	avs_iters += s_iters;
	++meascount;
	fflush(stdout);
      }
      node0_printf("trajectory finished!!!!\n\n");
    }	/* end loop over trajectories */
    
    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
    if(meascount>0)  {
      node0_printf("average cg iters for step= %e\n",
		   (double)avs_iters/meascount);
#ifdef SPECTRUM
      node0_printf("average cg iters for spectrum = %e\n",
		   (double)avspect_iters/meascount);
#endif
    }
    node0_printf("ACCEPTS: %d\nREJECTS: %d\nRATIO: %e\n", num_accepts, 
	   num_rejects, (double)num_accepts/trajecs);
    dtime += dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",dtime);
      printf("total_iters = %d\n",total_iters);
    }
    fflush(stdout);
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      rephase( OFF );
      temp_field = (complex *)malloc(sites_on_node*4*sizeof(complex));
      if(temp_field == NULL) {
	printf("Malloc failed to create temp_field in control\n");
	exit(1);
      }
#ifdef NON_COMPACT
      /* for non_compact gauge action also saving real gauge potential
	 in seperate file savefile2 */
      temp_pot = (Real *)malloc(sites_on_node*4*sizeof(Real));
      if(temp_pot == NULL) {
	printf("Malloc failed to create temp_pot in control\n");
	exit(1);
      }
      //put gauge potential into field structure
      FORALLSITES(i,s) {
	for(dir=XUP;dir<=TUP;dir++){
	  temp = temp_pot + 4*i + dir;
	  *temp = s->potential[dir];
	}
      }
      char default_file_xml_2[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>MILC ILDG archival gauge configuration</title>";
      char recinfo_2[] = "TEST";
      save_real_scidac_from_field( savefile2, default_file_xml_2, recinfo_2, 
			       QIO_SINGLEFILE, temp_pot, 4);
      free(temp_pot);
#endif
      //printf("inside save lattice\n");
      //Put gauge links into field structure
      FORALLSITES(i,s) {
	for(dir=XUP;dir<=TUP;dir++){
	  temp_link = temp_field + 4*i + dir;
	  temp_link->real = s->link[dir].real;
	  temp_link->imag = s->link[dir].imag;
	}
      }
      //printf("before calling scidac routine\n");
      //write field to output file
      char default_file_xml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>MILC ILDG archival gauge configuration</title>";
      char recinfo[] = "TEST";
      save_complex_scidac_from_field( savefile, default_file_xml, recinfo, 
			       QIO_SINGLEFILE, QIO_SERIAL, temp_field, 4);
      node0_printf("after calling scidac routine\n");
      /* NEED TO KNOW FOLLOWING ARGUMENTS:
	 char *fileinfo, char *recinfo */
      free(temp_field);
      /* modify reading and writing routines later 11/15/2010 C.W. */
      //save_lattice( saveflag, savefile, stringLFN );
      rephase( ON );
#ifdef HAVE_QIO
//       save_random_state_scidac_from_site("randsave", "Dummy file XML",
//        "Random number state", QIO_SINGLEFILE, F_OFFSET(site_prn));
//       save_color_vector_scidac_from_site("xxx1save", "Dummy file XML",
//        "xxx vector", QIO_SINGLEFILE, F_OFFSET(xxx1),1);
//       save_color_vector_scidac_from_site("xxx2save", "Dummy file XML",
//        "xxx vector", QIO_SINGLEFILE, F_OFFSET(xxx2),1);
#endif
    }
  }
#ifdef HAVE_QDP
  QDP_finalize();
#endif  
  normal_exit(0);
  return 0;
}
