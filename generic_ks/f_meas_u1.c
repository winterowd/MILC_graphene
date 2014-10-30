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

void f_meas_imp_u1( field_offset phi_off, field_offset xxx_off, Real mass,
		    ferm_links_u1_t *fn, ferm_links_u1_t *fn_dmdu0){
    Real r_psi_bar_psi_even, i_psi_bar_psi_even;
    Real  r_psi_bar_psi_odd, i_psi_bar_psi_odd;
    Real r_ferm_action;
    /* local variables for accumulators */
    register int i;
    register site *st;
    double rfaction;
    double_complex pbp_e, pbp_o;
    //double_complex temp_e, temp_o, xsc_e, xsc_o;
    complex *t_fatlink, *t_longlink, temp_mul;
    //complex *g_rand_temp, *temp_invert, *zzz;
    quark_invert_control qic;
    int my_volume;
#ifdef NPBP_REPS
    double pbp_pbp;
#endif
    complex cc;
#ifdef FOUR_DIM
    my_volume=volume;
#else
    my_volume=nx*nt*ny;
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

    for(jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++){
      rfaction = (double)0.0;
      pbp_e = pbp_o = dcmplx((double)0.0,(double)0.0);
      
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
      /* phi_off = M g_rand (still) */
      /* theta_off = M g_rand_temp
      /* xxx_off = M^{-1} g_rand */
      /* NEED TO DOUBLE CHECK 7/27/2011 */
      clear_latvec_u1(xxx_off,EVENANDODD);
      mat_invert_uml_u1( F_OFFSET(g_rand), xxx_off, phi_off, mass, 
			       prec, fn );     

      
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
      r_ferm_action =  rfaction*(1.0/(double)my_volume) ;
      node0_printf("PBP: mass %e     %e  %e  %e  %e ( %d of %d )\n", mass,
		   r_psi_bar_psi_even, r_psi_bar_psi_odd,
		   i_psi_bar_psi_even, i_psi_bar_psi_odd,
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

