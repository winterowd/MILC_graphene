/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This file is for lattice generation with the RHMC algorithm */

#define DEBUG_INVERTER
#define CONTROL
#include "graphene_includes.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"

#ifdef MILC_GLOBAL_DEBUG
#include "debug.h"
#endif /* MILC_GLOBAL_DEBUG */

/* For information */
#define NULL_FP -1

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int
main( int argc, char **argv )
{
  complex **t_fl = &fn_links_u1.fat;
  complex **t_ll = &fn_links_u1.lng;
  register complex *fat1;
  register complex *long1;
  register complex *src;
  register complex *dest;
  quark_invert_control qic;
  double norm = 0.0;
  int dir;
  register site *s;
  int i,meascount,traj_done;
  int prompt;
  int s_iters, avs_iters, avspect_iters, avbcorr_iters;
  double dtime, dclock();
  
  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#ifndef QDP_PROFILE
  QDP_profcontrol(0);
#endif
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();

#ifdef DEBUG_INVERTER
  printf("Need %zu for source and dest allocation\n", 
	2*sites_on_node * sizeof(complex));
  /*Allocate for test source and destination vector */
  src = (complex *)malloc(sites_on_node * sizeof(complex));
  dest = (complex *)malloc(sites_on_node * sizeof(complex));

  /*setup for cg routine*/
  qic.max = 100;
  qic.nrestart = 2;
  qic.resid = .001;
  qic.relresid = .001;

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

  //set values on each each site and initialize destination vector to zero
  FORALLMYSITES(i,s){
    
    src[i].real = 0.5;
    src[i].imag = 0.5;
    dest[i].real = 0.0;
    dest[i].imag = 0.0;
    //printf("At site (%d, %d, %d, %d)\n", s->x, s->y, s->z, s->t);
    //printf("Source vector: Real %f, Imag %f\n", src->real, src->imag);
  }

  FORMYODDSITES(i,s){
    norm += (double)cabs( (src+i) );
  } 
  printf("norm of source on odd sites: %f\n", norm);
  //for testing purposes setup unit fat and long links
  FORALLMYUPDIR(dir) {
    FORALLMYSITES(i,s){
      //printf("inside test loop \n");
      fat1 = (*t_fl) + 4*i + dir;
      long1 = (*t_ll) + 4*i + dir;
      fat1->real = 0.0; 
      fat1->imag = 0.0;
      long1->real = 0.0; 
      long1->imag = 0.0;
      //output for testing purposes
      //printf("At site (%d, %d, %d, %d) in direction %d\n", s->x, s->y, s->z, s->t, dir); 
      //printf("Long link: Real %f, Imag %f\n", long1->real, long1->imag);
      //printf("Fat link: Real %f, Imag %f\n", fat1->real, fat1->imag);
    }
  }
  //fermion links set to valid
  fn_links_u1.valid = 1;
 //call matrix inversion routine
  //printf("HI\n");
  //mat_invert_uml_field_u1(src, dest, &qic, 1.0, &fn_links_u1);
  init_path_table(&ks_act_paths);
  int construct = make_path_table(&ks_act_paths, NULL);
  printf("Path table constructed: %d, Num paths: %d\n", construct,
	 ks_act_paths.num_q_paths);
  funnylat_u1();
  phases_in = ON;
  //load_ferm_links_u1(&fn_links_u1, &ks_act_paths);
  load_longlinks_u1(&fn_links_u1, &ks_act_paths);
  load_fatlinks_u1(&fn_links_u1, &ks_act_paths);
  FORALLMYUPDIR(dir) {
    FORALLMYSITES(i,s){
      fat1 = (*t_fl) + 4*i + dir;
      long1 = (*t_ll) + 4*i + dir;
      printf("At site (%d, %d, %d, %d) in direction %d\n", 
	     s->x, s->y, s->z, s->t, dir);
      printf("Fat Link in %d direction, Real: %f, Imag: %f\n", dir,
	     fat1->real, fat1->imag);
      printf("Long Link in %d direction, Real: %f, Imag: %f\n", dir,
	     long1->real, long1->imag);
    }
  }
  /*
  FORALLMYSITES(i,s){
    printf("At site (%d, %d, %d, %d), index %d\n", s->x, s->y, s->z, s->t, i);
    printf("Destination vector: Real %f, Imag %f\n", dest[i].real, 
	   dest[i].imag);
  }
  */

  //end here for debugging purposes
  return 0;

#endif
    
  /* loop over input sets */
  while( readin(prompt) == 0) {
    
    /* perform warmup trajectories */
#ifdef MILC_GLOBAL_DEBUG
    global_current_time_step = 0;
#endif /* MILC_GLOBAL_DEBUG */

    dtime = -dclock();
    for( traj_done=0; traj_done < warms; traj_done++ ){
      //update();  
    }
    node0_printf("WARMUPS COMPLETED\n"); fflush(stdout);
    
    /* perform measuring trajectories, reunitarizing and measuring 	*/
    meascount=0;		/* number of measurements 		*/
    avspect_iters = avs_iters = avbcorr_iters = 0;

    for( traj_done=0; traj_done < trajecs; traj_done++ ){ 
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
  {
  int isite, idir;
  site *s;
  FORALLSITES(isite,s) {
    for( idir=XUP;idir<=TUP;idir++ ) {
      lattice[isite].on_step_Y[idir] = 0;
      lattice[isite].on_step_W[idir] = 0;
      lattice[isite].on_step_V[idir] = 0;
    }
  }
  }
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
      /* do the trajectories */
      //s_iters=update(); 

      /* measure every "propinterval" trajectories */
      if( (traj_done%propinterval)==(propinterval-1) ){
	
	/* call gauge_variable fermion_variable measuring routines */
	/* results are printed in output file */
	//rephase(OFF); COMMENT OUT FOR TIME BEING
	//g_measure( );
	//rephase(ON);
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ
        g_measure_plaq( );
#endif
#ifdef MEASURE_AND_TUNE_HISQ
        g_measure_tune( );
#endif /* MEASURE_AND_TUNE_HISQ */
#endif /* MILC_GLOBAL_DEBUG */


	/**************************************************************/
	/* Compute chiral condensate and related quantities           */
	
	/* Make fermion links if not already done */
	
	for(i = 0; i < par_buf.num_pbp_masses; i++){
#ifdef HISQ
	  fn_links.hl.current_X_set = par_buf.ksp_pbp[i].naik_term_epsilon_index;
#endif
	  //load_ferm_links_u1(&fn_links, &ks_act_paths); 
/*why did DeTar delete all previous ks_act_paths*? */ 
#ifdef DM_DU0
#ifdef HISQ
	  fn_links_dmdu0.hl.current_X_set = 
	    par_buf.ksp_pbp[i].naik_term_epsilon_index;
#endif
	  //load_ferm_links_u1(&fn_links_dmdu0, &ks_act_paths_dmdu0 );
#endif
	  /*f_meas_imp_field( par_buf.npbp_reps, &par_buf.qic_pbp[i], 
	    par_buf.ksp_pbp[i].mass, &fn_links, &fn_links_dmdu0); gave error when compiling 4/20/2010 */

#ifdef D_CHEM_POT
	  Deriv_O6_field( par_buf.npbp_reps, &par_buf.qic_pbp[i],
			  par_buf.ksp_pbp[i].mass, &fn_links, &fn_links_dmdu0);
#endif
	}

	avs_iters += s_iters;
	++meascount;
	fflush(stdout);
      }
    }	/* end loop over trajectories */
    
    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
    if(meascount>0)  {
      node0_printf("average cg iters for step= %e\n",
		   (double)avs_iters/meascount);
    }
    
    dtime += dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",dtime);
      printf("total_iters = %d\n",total_iters);
    }
    fflush(stdout);
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      rephase( OFF );
      //save_lattice( saveflag, savefile, stringLFN );
      rephase( ON );
    }
  }
#ifdef HAVE_QDP
  QDP_finalize();
#endif  
  normal_exit(0);
  return 0;
}

