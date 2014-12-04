/************************** d_plaq4_u1.c *******************************/
/* MIMD version 7 */
/* This version mallocs the temporary su3_matrix */

/* Double precision version of "plaquette4.c" including optional
   Schroedinger functional - UMH - 1/27/00 */

/* Measure the average plaquette of the space-space and
   space-time plaquettes */

/* Updated to U(1) 10/10 C.W. */

#include "generic_includes.h"

void d_plaquette_u1(double *ss_plaq,double *st_plaq) {
/* su3mat is scratch space of size su3_matrix */
complex *temp_complex;
 register int i,j,k,index,dir1,dir2;
register site *s;
register complex *m1,*m4;
 complex mtmp, temp_mul;
 double ss_sum,st_sum, *st_site, temp;
 double *ss_debug, *st_debug;
double ss1[50] = {0.}; double ss2[50] = {0.}; 
 double st1[50] = {0.}; double st2[50] = {0.}; double st3[50] = {0.};
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;

    temp_complex = (complex *)malloc(sizeof(complex)*sites_on_node);
    st_site = (double *)malloc(sizeof(double)*sites_on_node);
    ss_debug = (double *)malloc(nt*sizeof(double));
    st_debug = (double *)malloc(nt*sizeof(double));
    for(i=0; i<nt; i++) 
      ss_debug[i] = st_debug[i] = 0.;

    if(temp_complex == NULL)
      {
	printf("plaquette: can't malloc temp_complex\n");
	fflush(stdout); terminate(1);
      }

    for(dir1=YUP;dir1<=TUP;dir1++){
	for(dir2=XUP;dir2<dir1;dir2++){
	  
	  //node0_printf("BEFORE GATHER D_PLAQ4_U1\n");
	  mtag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(complex),
				     dir1, EVENANDODD, gen_pt[0] );
	  mtag1 = start_gather_site( F_OFFSET(link[dir1]), sizeof(complex),
				     dir2, EVENANDODD, gen_pt[1] );
	  //node0_printf("AFTER GATHER D_PLAQ4_U1\n");
	  FORALLSITES(i,s){
	    m1 = &(s->link[dir1]);
	    m4 = &(s->link[dir2]);
	    CMULJ_( *m4, *m1, temp_complex[i] );
	    //mult_su3_an(m4,m1,&su3mat[i]);
	  }
	  //node0_printf("BEFORE WAIT D_PLAQ4_U1\n");
	  wait_gather(mtag0);
	  wait_gather(mtag1);
	  //node0_printf("AFTER WAIT D_PLAQ4_U1\n");
	  FORALLSITES(i,s){
#ifdef SCHROED_FUN
		if(dir1==TUP ){
		  if(s->t==(nt-1)){
		    CMUL( temp_complex[i], s->boundary[dir2], mtmp );
		    /*mult_su3_nn( &su3mat[i],
		      &(s->boundary[dir2]), &mtmp); */
		  }
		  else{
		    CMUL( temp_complex[i], *(complex *)(gen_pt[0][i]),
			  mtmp );
		    /*mult_su3_nn( &su3mat[i],
			(su3_matrix *)(gen_pt[0][i]), &mtmp); */
		  }
		  CMULJ_( *(complex *)(gen_pt[1][i]), mtmp, temp_mul );
		  st_sum += (double)temp_mul.real;
		  /* st_sum +=
		     realtrace_su3((su3_matrix *)(gen_pt[1][i]), &mtmp); */
		}
		else if(s->t > 0){
		  CMUL( temp_complex[i], *(complex *)(gen_pt[0][i]), 
			mtmp );
		  CMULJ_( *(complex *)(gen_pt[1][i]), mtmp, temp_mul );
		  ss_sum += (double)temp_mul.real;
		  /*mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
		    &mtmp);
		    ss_sum +=
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]), &mtmp); */
		}
#else
		CMUL( temp_complex[i], *(complex *)(gen_pt[0][i]), mtmp );
		/* mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
		   &mtmp); */

		if(dir1==TUP ) {
		  CMULJ_( *(complex *)(gen_pt[1][i]), mtmp, temp_mul );
		  st_debug[s->t] += (double)temp_mul.real;
		  st_sum += (double)temp_mul.real;
		  if(dir2==XUP) {
		    st1[s->z] += (double)temp_mul.real;
		    if(s->z == 0)
		      st_site[i] = (double)temp_mul.real;
		  }
		  if(dir2==YUP) {
		    st2[s->z] += (double)temp_mul.real;
		    if(s->z == 0)
		      st_site[i] += (double)temp_mul.real;
		  }
		  if(dir2==ZUP) {
		    st3[s->z] += (double)temp_mul.real;
		    if(s->z == 0)
		      st_site[i] += (double)temp_mul.real;
		  }
		  st_site[i] = (1./3.)*st_site[i];
		  /*st_sum += (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp); */
		}
		else {          
		  CMULJ_( *(complex *)(gen_pt[1][i]), mtmp, temp_mul );
		  ss_debug[s->t] += (double)temp_mul.real;
		  ss_sum += (double)temp_mul.real;
		  if(dir1 == ZUP || dir2 == ZUP)
		    ss1[s->z] += (double)temp_mul.real;
		  else
		    ss2[s->z] += (double)temp_mul.real;
		  /* ss_sum += (double)
		     realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp); */
		  //if(s->z >= 5)
		  //node0_printf("Z > 5!!!!\n");
		}
#endif
	    }

	    cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	}
    }
    g_doublesum( &ss_sum );
    g_doublesum( &st_sum );


    //store average of st plaq at (x,y) for all t at t=0  
    /*
    for(i=0; i<nx; i++){
      for(j=0; j<ny; j++) {
	temp=0.;
	for(k=0; k<nt; k++){
	  index = node_index(i,j,0,k);
	  temp += st_site[index];
	}
	index = node_index(i,j,0,0);
	st_site[index] = 1./((double)nt)*(temp);
      }
    }
    
    //print out st at each plaquette on z=0 plane
    FORALLSITES(i,s) if(s->z==0 && s->t==0)
      printf("TEST_ST: x: %d y: %d st: %e\n", s->x, s->y, st_site[i]);
    */

#ifdef SCHROED_FUN
    *ss_plaq = ss_sum /((Real)(3*nx*ny*nz*(nt-1)));
#else
    *ss_plaq = ss_sum /((Real)(3*nx*ny*nz*nt));
#endif
    *st_plaq = st_sum /((double)(3*nx*ny*nz*nt));
    //print out the xz/yz oriented and xy oriented ss plaquettes
    for(i=0; i<nz; i++) { 
      g_doublesum(&ss1[i]);
      g_doublesum(&ss2[i]);
      g_doublesum(&st1[i]);
      g_doublesum(&st2[i]);
      g_doublesum(&st3[i]);
      ss1[i] *= 1./(2.*nx*ny*nt);
      ss2[i] *= 1./(nx*ny*nt);
      st1[i] *= 1./(nx*ny*nt);
      st2[i] *= 1./(nx*ny*nt);
      st3[i] *= 1./(nx*ny*nt);
      node0_printf("PLAQS: %e %e %e %e %e\n", ss1[i], ss2[i], st1[i], st2[i],
		   st3[i]);
    } 
    for(i=0; i<nt; i++) {
      g_doublesum(&ss_debug[i]);
      g_doublesum(&st_debug[i]);
      node0_printf("DEBUG_PLAQS: %d %e %e\n", i, 
		   (1./((double)(3*nx*ny*nz)))*ss_debug[i], 
		   (1./((double)(3*nx*ny*nz)))*st_debug[i]); 
    }
    free(ss_debug); free(st_debug);
    free(temp_complex);
    free(st_site);
} /* d_plaquette4 */

