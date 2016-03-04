/**************** f_meas_u1.c ***************************************/
/* MIMD version 7 */
/* UH 11/01 Include measurements for quark number susceptibilities */
/*          Note: this does not work for p4-action! */
/* UH 11/1/01 write complex stochastic estimators */
/* TB 10/01 Include measurements of dM/du0 for EOS */
/* CD 7/14/01 allow for multiple stochastic estimators NPBP_REPS */
/* DT 12/97 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.
*/

/* Measure fermionic observables:
    psi-bar-psi (separately on even and odd sites)
    fermion action

    This routine uses g_rand, phi, xxx, and other vectors used by
    the matrix inversion.
*/

#include "generic_ks_includes_u1.h"	/* definitions files and prototypes */

void shift_field_haldane(int dir, complex *src, complex *dest, complex *links, int forward) {

  register int i;
  register site *s;
  msg_tag *tag[2];
  complex *tvec;

  tvec = (complex *)malloc(sites_on_node*sizeof(complex));
  if(forward==1) {
    
    tag[0] = start_gather_field( src, sizeof(complex), dir, EVENANDODD, gen_pt[0] );
    wait_gather(tag[0]);
    
    FORALLMYSITES(i, s) {
      CMUL(links[4*i+dir], *(complex *)gen_pt[0][i], dest[i]);
    }
  }
  else{
  
    
    FORALLMYSITES(i, s) {    
      CMULJ_(links[4*i+dir], src[i], tvec[i]);
    }

    tag[1] = start_gather_field(tvec, sizeof(complex), OPP_DIR(dir), EVENANDODD, gen_pt[1]);
    
    wait_gather(tag[1]);

    FORALLMYSITES(i, s) {
      dest[i].real = *(complex *)gen_pt[1][i].real;
      dest[i].imag = *(complex *)gen_pt[1][i].imag;
    }

  }

  free(tvec);
  cleanup_gather(tag[0]);
  cleanup_gather(tag[1]);
  
}//sym_shift_field()

void shift_field_path_haldane(int n, int *d, complex *src, complex *dest, complex *links, int *forward) {

  int i, j;
  site *s; 
  complex *tvec;
  
  tvec = (complex *)malloc(sites_on_node*sizeof(complex));

  for(j=0; j<n; j++) { //loop over path length
    if(j==0)
      shift_field_haldane(d[j], src, tvec, links, forward[j]);
    else
      shift_field_haldane(d[j], dest, tvec, links, forward[j]);
    FORALLMYSITES(i, s) {
      dest[i].real = tvec[i].real;
      dest[i].imag = tvec[i].imag;
    }
  }

  free(tvec);

}//shift_field_path()

void three_link_shift_haldane(complex *src, complex *dest, complex *links, int *eta) {

  register int i, j;
  register site *s;
  int *forward[3];
  complex *tvec;

  tvec = (complex *)malloc(sites_on_node*sizeof(complex));

  /* All permutation with appropriate sign */
  struct {
    int d[3];
    Real sign ;
  } p[6]={{{XUP,YUP,TUP},+1.0/6.0},
	  {{YUP,TUP,XUP},+1.0/6.0},
	  {{TUP,XUP,YUP},+1.0/6.0},
	  {{XUP,TUP,YUP},+1.0/6.0},
	  {{YUP,XUP,TUP},+1.0/6.0},
	  {{TUP,YUP,XUP},+1.0/6.0}}; /* The factor of 6 accounts for the *
				* multiplicity of the permutations */
  for(i=0; i<3; i++) {
    if(eta[i]==0)
      forward[i]=1;
    else
      forward[i]=0;
  }
  
  FORALLSITES(i, s) {
    dest[i].real = dest[i].imag = 0.;
  }
  for(j=0; j<6; j++) { //loop over paths
    shift_field_path_haldane(3, p[j].d, src, tvec, links, forward);
    FORALLMYSITES(i, s) {
      CMULREAL(tvec[i], p[j].sign, tvec[i]);
      CADD(dest[i], tvec[i], dest[i]);
    }
  }
  
  free(tvec);

}//three_link_shift()

void f_meas_imp_u1( field_offset phi_off, field_offset xxx_off, Real mass,
		    ferm_links_u1_t *fn, ferm_links_u1_t *fn_dmdu0){
    Real r_psi_bar_psi_even, i_psi_bar_psi_even;
    Real r_psi_bar_psi_odd, i_psi_bar_psi_odd;
    Real r_haldane_even, i_haldane_even;
    Real r_haldane_odd, i_haldane_odd;
    Real r_ferm_action;
    /* local variables for accumulators */
    register int i, dir;
    register site *st;
    double rfaction;
    double_complex pbp_e, pbp_o;
    double_complex haldane_e, haldane_o;
    //double_complex temp_e, temp_o, xsc_e, xsc_o;
    complex *t_fatlink, *t_longlink, temp_mul;
    complex *temp_vec1, *temp_vec2, *temp_vec3;
    complex *links;
    //complex *g_rand_temp, *temp_invert, *zzz;
    quark_invert_control qic;
    int my_volume;
    int xdisp, ydisp, tdisp;
    int xcorner, ycorner, tcorner;
    int xoppcorner, yoppcorner, toppcorner;
    int eta[3];
#ifdef NPBP_REPS
    double pbp_pbp;
#endif
    complex cc;
#ifdef FOUR_DIM
    my_volume=volume;cleanup_gather(tag0);
      cleanup_gather(tag1);
#else
    my_volume=(double)(nx*nt*ny);
#endif

#ifdef DM_DU0
    double r_pb_dMdu_p_even, r_pb_dMdu_p_odd;
#endif

#ifdef CHEM_POT
#ifndef FN		/* FN is assumed for quark number susc. */
BOMB THE COMPILE
#endif
#ifndef NPBP_REPS	/* Need multiple repetitions for susceptibilities! */
BOMB THE COMPILE
#endif
    msg_tag *tag0, *tag1, *tag2, *tag3;
    double_complex pb_dMdmu_p_e, pb_dMdmu_p_o;
    double_complex pb_d2Mdmu2_p_e, pb_d2Mdmu2_p_o;
    Real r_pb_dMdmu_p_e, i_pb_dMdmu_p_e;
    Real r_pb_dMdmu_p_o, i_pb_dMdmu_p_o;
    Real r_pb_d2Mdmu2_p_e, i_pb_d2Mdmu2_p_e;
    Real r_pb_d2Mdmu2_p_o, i_pb_d2Mdmu2_p_o;
    double MidM_MidM;
#endif

#ifdef NPBP_REPS
    int npbp_reps = npbp_reps_in;  /* Number of repetitions of stochastic
                                   estimate */
    int prec = prec_pbp;  /* Precision of the inversion */
#else
    int npbp_reps = 1;   /* Default values */
    int prec = PRECISION;
#endif
    int jpbp_reps;

    t_fatlink = fn->fat;
    t_longlink = fn->lng;

    //g_rand_temp = (complex *)malloc(sites_on_node*sizeof(complex));
    //temp_invert = (complex *)malloc(sites_on_node*sizeof(complex));
    //zzz = (complex *)malloc(sites_on_node*sizeof(complex));

    /* set up inverter control */
    qic.prec       = prec;
    qic.max        = niter;
    qic.nrestart   = nrestart;
    qic.resid      = rsqprop;
    qic.relresid   = 0;

    //allocate vectors for haldane mass calculation
    temp_vec1 = (complex *)malloc(sites_on_node*sizeof(complex));
    temp_vec2 = (complex *)malloc(sites_on_node*sizeof(complex));
    temp_vec3 = (complex *)malloc(sites_on_node*sizeof(complex));
    links = (complex *)malloc(sites_on_node*4*sizeof(complex));

    //copy links to field vector without phases 
    //(assumes phases are in upon entry to function)
    FORALLUPDIR(dir) {
      FORALLSITES(i, st) {	
	links[4*i+dir].real = st->link[dir].real*st->phase[dir];
	links[4*i+dir].imag = st->link[dir].imag*st->phase[dir];
      }
    }

    for(jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++){

      for(xdisp=0; xdisp<nx; xdisp+=2) for(ydisp=0; ydisp<ny; ydisp+=2) for(tdisp=0; tdisp<nt; tdisp+=2) {
            for(xcorner=0, xoppcorner=1; xcorner<2; xcorner++, xoppcorner--) for(ycorner=0, yoppcorner=1; ycorner<2; ycorner++, yoppcorner--) for(tcorner=0, toppcorner=1; tcorner<2; tcorner++, toppcorner--) {
      
      rfaction = (double)0.0;
      pbp_e = pbp_o = dcmplx((double)0.0,(double)0.0);
      haldane_e = haldane_o = dcmplx((double)0.0,(double)0.0);
      
      /* Make random source, and do inversion */
      /* generate g_rand random; phi_off = M g_rand */
      /* and theta_off = M g_rand_temp */
#ifndef Z2RSOURCE
      //grsource_imp_u1( theta_off, mass, EVENANDODD, fn );
      /* save first gaussian random vector  put theta into field and 
	 clear other vector */
      /*FORALLMYSITES(i, st) {
	g_rand_temp[i].real = st->g_rand.real;
	g_rand_temp[i].imag = st->g_rand.imag;
	temp_invert[i].real = st->theta.real;
	temp_invert[i].imag = st->theta.imag;
	zzz[i].real = zzz[i].imag = 0.;
	}*/
      clear_latvec_u1( phi_off, EVENANDODD);
      grsource_imp_u1( phi_off, mass, EVENANDODD, fn );
#else
      //z2rsource_imp( phi_off, mass, EVENANDODD, fn ); LEAVE OUT FOR NOW 10/10
#endif
      FORALLSITES(i,st) { //copy g_rand to temp_vec1, clear temp_vec2 and temp_vec3
	eta[0] = st->x - 2*((int)(st->x)/2);
	eta[1] = st->y - 2*((int)(st->y)/2);
	eta[2] = st->t - 2*((int)(st->t)/2);
	//if( (eta[0]==xcorner) && (eta[1]==ycorner) && (eta[2]==tcorner) && (st->z==0)) { //only do source at one corner of cube for now (02/04/16)
	if( (st->x==(xdisp+xcorner)) && (st->y==(ydisp+ycorner)) && (st->t==(tdisp+tcorner)) && (st->z==0)) { 
	  //temp_vec1[i].real = st->g_rand.real;
	  //temp_vec1[i].imag = st->g_rand.imag;
	  temp_vec1[i].real = 1.0;
	  temp_vec1[i].imag = 0.0;
	  printf("point_source at %d %d %d\n", st->x, st->y, st->t);
          printf("opp_corner at %d %d %d\n", xdisp+xoppcorner, ydisp+yoppcorner, tdisp+toppcorner);
	  printf("eta %d %d %d\n", eta[0], eta[1], eta[2]);
	}
	else {
	  temp_vec1[i].real = 0.0;
	  temp_vec1[i].imag = 0.0;
	}
	temp_vec3[i].real = temp_vec2[i].real = temp_vec3[i].imag = temp_vec2[i].imag = 0.;
      }
      
      /* phi_off = M g_rand (still) */
      /* theta_off = M g_rand_temp
      /* xxx_off = M^{-1} g_rand */
      /* NEED TO DOUBLE CHECK 7/27/2011 */
      clear_latvec_u1(xxx_off,EVENANDODD);
      mat_invert_uml_u1( F_OFFSET(g_rand), xxx_off, phi_off, mass, 
			       prec, fn );     

      //call routine to shift temp_vec1 and put result in temp_vec2
      eta[0]=xcorner; eta[1]=ycorner; eta[2]=tcorner
      three_link_shift_haldane(temp_vec1, temp_vec2, links, eta);

      FORALLMYSITES(i,st) {
        if((st->x==(xdisp+xoppcorner)) && (st->y==(ydisp+yoppcorner)) && (st->t==(tdisp+toppcorner))) {
          printf("temp_vec2 %d %d %d %d %e %e\n", st->x, st->y, st->z, st->t, temp_vec2[i].real, temp_vec2[i].imag);
        }
      }
      
      //invert on shifted source
      mat_invert_uml_field_u1( temp_vec2, temp_vec3, &qic, mass, fn);
      

      /*complex *t_src, *t_dest;
      site *s;
      Real vmass_x2 = 2.*mass;
      t_src  = (complex *)malloc(sizeof(complex)*sites_on_node);
      t_dest = (complex *)malloc(sizeof(complex)*sites_on_node);
      FORALLMYSITES(i,s){
	t_src[i]  = *((complex *)F_PT(s,xxx_off) );
	t_dest[i].real = t_dest[i].imag = 0.0;
      }
      dslash_fn_field_u1( t_src, t_dest, EVENANDODD, fn);
      FORALLMYSITES(i,s){
	CMULREAL(t_src[i], vmass_x2, temp_mul);
	CADD(t_dest[i], temp_mul, t_dest[i]);
      }
      node0_printf("CHECK INVERSE RESULT:\n");
      FORALLMYSITES(i,s) {
	node0_printf("site %d %d %d %d\n t_dest: %e %e g_rand %e %e\n", s->x, 
		     s->y, s->z, s->t, t_dest[i].real, t_dest[i].imag,
		     s->g_rand.real, s->g_rand.imag);
	
      }
      free(t_src); free(t_dest);*/
      //mat_invert_uml_field_u1( g_rand_temp, zzz, &qic, mass, fn );
      //clear_latvec_u1(yyy_off,EVENANDODD);
      //mat_invert_uml_u1( xxx_off, yyy_off, theta_off, mass, prec, fn); 
#ifdef DM_DU0
      //r_pb_dMdu_p_even = r_pb_dMdu_p_odd = (double)0.0;
      /* dMdu_x = dM/du0 M^{-1} g_rand */
      /*ddslash_fn_du0_site( xxx_off, F_OFFSET(dMdu_x), EVENANDODD, 
	fn, fn_dmdu0 );  LEAVE OUT FOR NOW 10/10 */
#endif

#ifdef CHEM_POT
      pb_dMdmu_p_e = pb_dMdmu_p_o = dcmplx((double)0.0,(double)0.0);
      pb_d2Mdmu2_p_e = pb_d2Mdmu2_p_o = dcmplx((double)0.0,(double)0.0);

      /* Start gathers from positive t-direction */
      tag0 = start_gather_site( xxx_off, sizeof(complex), TUP,
	EVENANDODD, gen_pt[0] );
      tag1 = start_gather_site( xxx_off, sizeof(complex), T3UP,
	EVENANDODD, gen_pt[1] );

      FORALLMYSITES(i,st){
	CMULJ_( t_fatlink[4*i+TUP], *(complex *)F_PT(st,xxx_off), 
		st->templongvec[TUP] );
	/* mult_adj_su3_mat_vec( &(t_fatlink[4*i+TUP]),
	  (su3_vector *)F_PT(st,xxx_off), &(st->tempvec[TUP]) ); */
	CMULJ_( t_longlink[4*i+TUP], *(complex *)F_PT(st,xxx_off),
		st->templongvec[TUP] );
	/* mult_adj_su3_mat_vec( &(t_longlink[4*i+TUP]),
	   (su3_vector *)F_PT(st,xxx_off), &(st->templongvec[TUP]) ); */
      }

      /* Start gathers from negative t-direction */
      tag2 = start_gather_site( F_OFFSET(tempvec[TUP]), sizeof(complex),
	OPP_DIR(TUP), EVENANDODD, gen_pt[2] );
      tag3 = start_gather_site( F_OFFSET(templongvec[TUP]), sizeof(complex),
	OPP_3_DIR(T3UP), EVENANDODD, gen_pt[3] );

      /* Wait gathers from positive t-direction and multiply by matrix */
      wait_gather(tag0);
      wait_gather(tag1);

      FORALLMYSITES(i,st){
	CMUL( t_fatlink[4*i+TUP], *(complex *)gen_pt[0][i], 
	      st->tempvec[0] );
	/* mult_su3_mat_vec( &(t_fatlink[4*i+TUP]),
	   (su3_vector *)gen_pt[0][i], &(st->tempvec[0]) ); */
	CMUL( t_longlink[4*i+TUP], *(complex *)gen_pt[1][i],
	      (st->templongvec[0] );
	/* mult_su3_mat_vec( &(t_longlink[4*i+TUP]),
	   (su3_vector *)gen_pt[1][i], &(st->templongvec[0]) ); */
      }

      /* Wait gathers from negative t-direction */
      wait_gather(tag2);
      wait_gather(tag3);
#endif

      /* fermion action = phi.xxx */
      /* psi-bar-psi on even sites = g_rand.xxx (fermion only live
       on hyperplane z=0) C.W. 10/10 */
      /* (Tr{D^-1})^2 and Tr{D^-2} needed for chiral susc. */
      FORMYEVENSITES(i,st){
	CMULJ_( *(complex *)F_PT(st,phi_off), *(complex *)F_PT(st,xxx_off),
		cc );
	/* cc = su3_dot( (su3_vector *)F_PT(st,phi_off),
	   (su3_vector *)F_PT(st,xxx_off) ); */
	rfaction += cc.real;
	CMULJ_( st->g_rand, *(complex *)F_PT(st,xxx_off), cc);
	/* cc = su3_dot( &(st->g_rand), (su3_vector *)F_PT(st,xxx_off) ); */
	CSUM(pbp_e, cc);
	CMULJ_( temp_vec1[i], temp_vec3[i], cc);
	CSUM(haldane_e, cc)

	//CMULJ_( g_rand_temp[i], zzz[i], cc);
	//CSUM(temp_e, cc); //temp_e contains other term for Tr^2
#ifdef DM_DU0
	/* r_pb_dMdu_p_even = g_rand * dM/du0 M^{-1} g_rand |even*/
	CMULJ_( st->g_rand, st->dMdu_x, cc );
	/* cc = su3_dot( &(st->g_rand), &(st->dMdu_x) ); */
	r_pb_dMdu_p_even += cc.real;
#endif

#ifdef CHEM_POT
	CMULJ_( st->g_rand, st->tempvec[0], cc);
	//cc = su3_dot( &(st->g_rand), &(st->tempvec[0]) );
	CSUM(pb_dMdmu_p_e, cc);
	CSUM(pb_d2Mdmu2_p_e, cc);
	CMULJ_( st->g_rand, *(complex *)gen_pt[2][i], cc );
	//cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[2][i] );
	CSUM(pb_dMdmu_p_e, cc);
	CSUB(pb_d2Mdmu2_p_e, cc, pb_d2Mdmu2_p_e);
	CMULJ_( st->g_rand, st->templongvec[0], cc );
	//cc = su3_dot( &(st->g_rand), &(st->templongvec[0]) );
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_dMdmu_p_e, cc);
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_d2Mdmu2_p_e, cc);
	CMULJ_( st->g_rand, *(complex *)gen_pt[3][i], cc );
	//cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[3][i] );
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_dMdmu_p_e, cc);
	CMULREAL(cc, 3.0, cc);
	CSUB(pb_d2Mdmu2_p_e, cc, pb_d2Mdmu2_p_e);
	CADD( st->tempvec[0], *(complex *)gen_pt[2][i], st->tempvec[0] );
	/* add_su3_vector( &(st->tempvec[0]), (su3_vector *)gen_pt[2][i],
	   &(st->tempvec[0]) ); */
	CADD( st->templongvec[0], *(complex *)gen_pt[3][i], 
	      st->templongvec[0] );
	/* add_su3_vector( &(st->templongvec[0]), (su3_vector *)gen_pt[3][i],
	   &(st->templongvec[0]) ); */
	CMULREAL( st->templongvec[0], 3.0, temp_mul );
	CADD( temp_mul, st->tempvec[0], st->dM_M_inv );
	/* scalar_mult_add_su3_vector( &(st->tempvec[0]), &(st->templongvec[0]),
	   3.0, &(st->dM_M_inv) ); */
#endif
      }

      /* psi-bar-psi on odd sites */
      FORMYODDSITES(i,st){
	CMULJ_( st->g_rand, *(complex *)F_PT(st,xxx_off), cc );
	/* cc = su3_dot( &(st->g_rand), (su3_vector *)F_PT(st,xxx_off) ); */
	CSUM(pbp_o, cc);
	CMULJ_( temp_vec1[i], temp_vec3[i], cc );
	CSUM(haldane_o, cc);

	//CMULJ_( g_rand_temp[i], zzz[i], cc);
	//CSUM(temp_o, cc); //temp_o contains other term for Tr^2
#ifdef DM_DU0
	/* r_pb_dMdu_p_odd = g_rand * dM/du0 M^{-1} g_rand |odd*/
	CMULJ_( st->g_rand, st->dMdu_x, cc );
	/* cc = su3_dot( &(st->g_rand), &(st->dMdu_x) ); */
	r_pb_dMdu_p_odd += cc.real;
#endif

#ifdef CHEM_POT
	CMULJ_( st->g_rand, st->tempvec[0], cc);
	/*cc = su3_dot( &(st->g_rand), &(st->tempvec[0]) );*/
	CSUM(pb_dMdmu_p_o, cc);
	CSUM(pb_d2Mdmu2_p_o, cc);
	CMULJ_( st->g_rand, *(complex *)gen_pt[2][i], cc );
	/*cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[2][i] );*/
	CSUM(pb_dMdmu_p_o, cc);
	CSUB(pb_d2Mdmu2_p_o, cc, pb_d2Mdmu2_p_o);
	CMULJ_( st->g_rand, st->templongvec[0], cc );
	/*cc = su3_dot( &(st->g_rand), &(st->templongvec[0]) );*/
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_dMdmu_p_o, cc);
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_d2Mdmu2_p_o, cc);
	CMULJ_( st->g_rand, *(complex *)gen_pt[3][i], cc );
	/*cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[3][i] );*/
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_dMdmu_p_o, cc);
	CMULREAL(cc, 3.0, cc);
	CSUB(pb_d2Mdmu2_p_o, cc, pb_d2Mdmu2_p_o);
	CADD( st->tempvec[0], *(complex *)gen_pt[2][i], st->tempvec[0] );
	/*add_su3_vector( &(st->tempvec[0]), (su3_vector *)gen_pt[2][i],
	  &(st->tempvec[0]) );*/
	CADD( st->templongvec[0], *(complex *)gen_pt[3][i], 
	      st->templongvec[0] );
	/*add_su3_vector( &(st->templongvec[0]), (su3_vector *)gen_pt[3][i],
	  &(st->templongvec[0]) );*/
	CMULREAL( st->templongvec[0], 3.0, temp_mul );
	CADD( temp_mul, st->tempvec[0], st->dM_M_inv );
	/*scalar_mult_add_su3_vector( &(st->tempvec[0]), &(st->templongvec[0]),
	  3.0, &(st->dM_M_inv) );*/
#endif
      }
      //g_dcomplexsum( &temp_o );
      //g_dcomplexsum( &temp_e );      
      g_dcomplexsum( &pbp_o );
      g_dcomplexsum( &pbp_e );
      g_dcomplexsum( &haldane_o);
      g_dcomplexsum( &haldane_e);
      g_doublesum( &rfaction );
      
#ifdef DM_DU0
      g_doublesum( &r_pb_dMdu_p_even );
      g_doublesum( &r_pb_dMdu_p_odd );
      r_pb_dMdu_p_even *= (2.0/(double)my_volume);
      r_pb_dMdu_p_odd *= (2.0/(double)my_volume);
      node0_printf("PB_DMDU_P: mass %e  %e  %e ( %d of %d )\n", mass,
		   r_pb_dMdu_p_even, r_pb_dMdu_p_odd, jpbp_reps+1, npbp_reps);
#endif

      r_psi_bar_psi_odd =  pbp_o.real*(2.0/(double)my_volume) ;
      i_psi_bar_psi_odd =  pbp_o.imag*(2.0/(double)my_volume) ;
      r_psi_bar_psi_even =  pbp_e.real*(2.0/(double)my_volume) ;
      i_psi_bar_psi_even =  pbp_e.imag*(2.0/(double)my_volume) ;
      r_haldane_odd =  haldane_o.real*(1.0/(double)my_volume);
      i_haldane_odd =  haldane_o.imag*(1.0/(double)my_volume);
      r_haldane_even =  haldane_e.real*(1.0/(double)my_volume);
      i_haldane_even =  haldane_e.imag*(1.0/(double)my_volume);
      r_ferm_action =  rfaction*(1.0/(double)my_volume) ;
      node0_printf("PBP: mass %e     %e  %e  %e  %e ( %d of %d )\n", mass,
		   r_psi_bar_psi_even, r_psi_bar_psi_odd,
		   i_psi_bar_psi_even, i_psi_bar_psi_odd,
		   jpbp_reps+1, npbp_reps);
      node0_printf("HALDANE: mass %e     %e  %e  %e  %e ( %d of %d )\n", mass,
		   r_haldane_even, r_haldane_odd,
		   i_haldane_even, i_haldane_odd,
		   jpbp_reps+1, npbp_reps);
      node0_printf("FACTION: mass = %e,  %e ( %d of %d )\n", mass,
		   r_ferm_action, jpbp_reps+1, npbp_reps);

#ifdef CHEM_POT
      /* free up the buffers */
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
      cleanup_gather(tag3);

      g_dcomplexsum( &pb_dMdmu_p_e );
      g_dcomplexsum( &pb_dMdmu_p_o );
      g_dcomplexsum( &pb_d2Mdmu2_p_e );
      g_dcomplexsum( &pb_d2Mdmu2_p_o );

      r_pb_dMdmu_p_e =  pb_dMdmu_p_e.real*(2.0/(double)my_volume) ;
      i_pb_dMdmu_p_e =  pb_dMdmu_p_e.imag*(2.0/(double)my_volume) ;
      r_pb_dMdmu_p_o =  pb_dMdmu_p_o.real*(2.0/(double)my_volume) ;
      i_pb_dMdmu_p_o =  pb_dMdmu_p_o.imag*(2.0/(double)my_volume) ;
      r_pb_d2Mdmu2_p_e =  pb_d2Mdmu2_p_e.real*(2.0/(double)my_volume) ;
      i_pb_d2Mdmu2_p_e =  pb_d2Mdmu2_p_e.imag*(2.0/(double)my_volume) ;
      r_pb_d2Mdmu2_p_o =  pb_d2Mdmu2_p_o.real*(2.0/(double)my_volume) ;
      i_pb_d2Mdmu2_p_o =  pb_d2Mdmu2_p_o.imag*(2.0/(double)my_volume) ;
      node0_printf("PB_DMDMU_P: mass %e  %e  %e  %e  %e ( %d of %d )\n", mass,
		   r_pb_dMdmu_p_e, r_pb_dMdmu_p_o,
		   i_pb_dMdmu_p_e, i_pb_dMdmu_p_o,
		   jpbp_reps+1, npbp_reps);
      node0_printf("PB_D2MDMU2_P: mass %e  %e  %e  %e  %e ( %d of %d )\n", mass,
		   r_pb_d2Mdmu2_p_e, r_pb_d2Mdmu2_p_o,
		   i_pb_d2Mdmu2_p_e, i_pb_d2Mdmu2_p_o,
		   jpbp_reps+1, npbp_reps);
#endif

#ifdef NPBP_REPS
      pbp_pbp = (double)0.0;
      clear_latvec_u1(F_OFFSET(M_inv), EVENANDODD);
      FORALLMYSITES(i,st){
	CMULREAL( *(complex *)F_PT(st,xxx_off), 1.0, st->M_inv );
	/* su3vec_copy( (su3_vector *)F_PT(st,xxx_off), &(st->M_inv) ); */
      }
      clear_latvec_u1(xxx_off,EVENANDODD);
      mat_invert_uml_u1( F_OFFSET(M_inv), xxx_off, phi_off, mass, 
			       prec, fn );
      FORALLMYSITES(i,st){
	CMULJ_( st->g_rand, *(complex *)F_PT(st,xxx_off), cc );
	/*cc = su3_dot( &(st->g_rand), (su3_vector *)F_PT(st,xxx_off) );*/
	pbp_pbp += cc.real;
      }
      g_doublesum( &pbp_pbp );
      pbp_pbp =  pbp_pbp*(1.0/(double)my_volume) ;
      node0_printf("TR_MM_INV: mass %e,  %e ( %d of %d )\n", mass,
		   pbp_pbp, jpbp_reps+1, npbp_reps);
#endif

#ifdef CHEM_POT
      clear_latvec_u1(xxx_off,EVENANDODD);
      mat_invert_uml_u1( F_OFFSET(dM_M_inv), xxx_off, phi_off, mass, 
			       prec, fn );

      /* Start gathers from positive t-direction */
      tag0 = start_gather_site( xxx_off, sizeof(complex), TUP,
	EVENANDODD, gen_pt[0] );
      tag1 = start_gather_site( xxx_off, sizeof(complex), T3UP,
	EVENANDODD, gen_pt[1] );

      FORALLMYSITES(i,st){
	CMULJ_( t_fatlink[4*i+TUP], *(complex *)F_PT(st,xxx_off), 
		st->tempvec[TUP] );
	/* mult_adj_su3_mat_vec( &(t_fatlink[4*i+TUP]),
	   (su3_vector *)F_PT(st,xxx_off), &(st->tempvec[TUP]) ); */
	CMULJ_( t_longlink[4*i+TUP], *(complex *)F_PT(st,xxx_off),
		st->templongvec[TUP] );
	/*mult_adj_su3_mat_vec( &(t_longlink[4*i+TUP]),
	  (su3_vector *)F_PT(st,xxx_off), &(st->templongvec[TUP]) );*/
      }

      /* Start gathers from negative t-direction */
      tag2 = start_gather_site( F_OFFSET(tempvec[TUP]), sizeof(complex),
	OPP_DIR(TUP), EVENANDODD, gen_pt[2] );
      tag3 = start_gather_site( F_OFFSET(templongvec[TUP]), sizeof(complex),
	OPP_3_DIR(T3UP), EVENANDODD, gen_pt[3] );

      /* Wait gathers from positive t-direction and multiply by matrix */
      wait_gather(tag0);
      wait_gather(tag1);

      FORALLMYSITES(i,st){
	CMUL( t_fatlink[4*i+TUP], *(complex *)gen_pt[0][i], st->tempvec[0] );
	/*mult_su3_mat_vec( &(t_fatlink[4*i+TUP]),
	  (su3_vector *)gen_pt[0][i], &(st->tempvec[0]) );*/
	CMUL( t_longlink[4*i+TUP], *(complex *)gen_pt[1][i], 
	      st->templongvec[0] );
	/*mult_su3_mat_vec( &(t_longlink[4*i+TUP]),
	  (su3_vector *)gen_pt[1][i], &(st->templongvec[0]) );*/
      }

      /* Wait gathers from negative t-direction */
      wait_gather(tag2);
      wait_gather(tag3);

      MidM_MidM = (double)0.0;

      FORALLMYSITES(i,st){
	CMULJ_( st->g_rand, st->tempvec[0], cc );
	/*cc = su3_dot( &(st->g_rand), &(st->tempvec[0]) );*/
	MidM_MidM += cc.real;
	CMULJ_( st->g_rand, *(complex *)gen_pt[2][i], cc );
	/*cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[2][i] ); */
	MidM_MidM += cc.real;
	CMULJ_( st->g_rand, st->templongvec[0], cc );
	/*cc = su3_dot( &(st->g_rand), &(st->templongvec[0]) );*/
	MidM_MidM += 3.0 * cc.real;
	CMULJ_( st->g_rand, *(complex *)gen_pt[3][i], cc );
	/*cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[3][i] );*/
	MidM_MidM += 3.0 * cc.real;
      }

      g_doublesum( &MidM_MidM );
      MidM_MidM =  MidM_MidM*(1.0/(double)my_volume) ;
      node0_printf("TR_MidM_MidM: mass %e,  %e ( %d of %d )\n", mass,
		   MidM_MidM, jpbp_reps+1, npbp_reps);

      /* free up the buffers */
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
      cleanup_gather(tag3);
#endif
      }
      }
    }

      free(temp_vec1);
      free(temp_vec2);
      free(temp_vec3);
      free(links);
}

