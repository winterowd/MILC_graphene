/******************** control.c *****************************/
/* MIMD version 7 */
/* Main procedure for U1 eigenvalues with improved dynamical fermions */

#define CONTROL
#include "ks_eig_includes.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"
#define LOOPEND
#include "../include/loopend.h"

EXTERN  gauge_header start_lat_hdr;     /* Input gauge field header */

int main( int argc, char **argv ){
  register site *s;
  int i, si;
  int prompt;
  double dtime;
  complex **eigVec ;
  complex *tmp ;
  double *eigVal ;
  int total_R_iters ;
  double chirality ;
  double eigen_sum;
  double eigen_sqr_sum;
  double IPR, rEV, iEV;

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
  
  /* loop over input sets */
  while( readin(prompt) == 0){
     
    dtime = -dclock();
    invalidate_all_ferm_links_u1(&fn_links);
    make_path_table(&ks_act_paths, &ks_act_paths_dmdu0);
    /* Load fat and long links for fermion measurements if needed */
    load_ferm_links_u1(&fn_links, &ks_act_paths);
    /* call fermion_variable measuring routines */
    /* results are printed in output file */
    f_meas_imp_u1( F_OFFSET(phi), F_OFFSET(xxx), mass,
		&fn_links, &fn_links_dmdu0);
    eigVal = (double *)malloc(Nvecs*sizeof(double));
    eigVec = (complex **)malloc(Nvecs*sizeof(complex*));
    for(i=0;i<Nvecs;i++)
      eigVec[i]=
	(complex*)malloc(sites_on_node*sizeof(complex));
    
    total_R_iters=Kalkreuter_u1(eigVec, eigVal, eigenval_tol, 
			     error_decr, Nvecs, MaxIter, Restart, 
			     Kiters, EVEN, &fn_links) ;
    tmp = (complex*)malloc(sites_on_node*sizeof(complex));
    for(i=0;i<Nvecs;i++)
      { 
	/* Construct to odd part of the vector.                 *
	 * Note that the true odd part of the eigenvector is    *
	 *  i/sqrt(eigVal) Dslash Psi. But since I only compute *
	 * the chirality the i factor is irrelevant (-i)*i=1!!  */
	dslash_fn_field_u1(eigVec[i], tmp, ODD, &fn_links) ;
	FORMYSITESANDPARITY(si,s,ODD){ 
	  CMULREAL( tmp[si], 1.0/sqrt(eigVal[i]), eigVec[i][si]);
	  /*scalar_mult_su3_vector( &(tmp[si]),
				  1.0/sqrt(eigVal[i]), 
				  &(eigVec[i][si]) ) ;*/
	} END_LOOP
	
	measure_chirality_u1(eigVec[i], &chirality, EVENANDODD);
	/* Here I divide by 2 since the EVEN vector is normalized to
	 * 1. The EVENANDODD vector is normalized to 2. I could have
	 * normalized the EVENANDODD vector to 1 and then not devide
	 * by to.  The measure_chirality routine assumes vectors
	 * normalized to 1.  */

	//compute IPR
	IPR=0.0;
	FORALLMYSITES(si,s) {
	  rEV = eigVec[i][si].real; iEV = eigVec[i][si].imag;
	  IPR += (rEV*rEV + iEV*iEV)*(rEV*rEV + iEV*iEV);
	  /*node0_printf("eigenvec %d %d %d %d %d %e %e\n", i, si, s->x, s->y, s->t, 
	    eigVec[i][si].real, eigVec[i][si].imag); */
	} 
	node0_printf("IPR %d %e\n", i, IPR);
	node0_printf("Chirality(%i): %g\n",i,chirality/2) ;
      } 
    free(tmp);

    /* Check that sum of eigenvalues is zero and sum of the squares
       of the positive eigenvalues is equal to the volume */
    eigen_sum = eigen_sqr_sum = 0;
    for(i=0; i<Nvecs; i++) {
      //eigen_sum += eigVal[i];
      if(eigVal[i] > 0)
	eigen_sqr_sum += eigVal[i];
    }
    printf("eigen_sqr_sum = %e\n", eigen_sqr_sum );
    /**
       for(i=0;i<Nvecs;i++)
       {
       sprintf(label,"DENSITY(%i)",i) ;
       print_densities(eigVec[i], label, ny/2,nz/2,nt/2, EVEN) ;
       }
    **/
    for(i=0;i<Nvecs;i++)
      free(eigVec[i]) ;
    free(eigVec) ;
    free(eigVal) ;
#ifdef FN
    invalidate_all_ferm_links_u1(&fn_links);
#endif
    fflush(stdout);
    
    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
    dtime += dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",dtime);
      printf("total_iters = %d\n",total_iters);
      printf("total Rayleigh iters = %d\n",total_R_iters);
    }
    fflush(stdout);
  }
  return 0;
}

