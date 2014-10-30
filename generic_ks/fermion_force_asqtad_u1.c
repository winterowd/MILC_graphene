/****** fermion_force_asqtad_u1.c  -- ******************/
/* Added 08/01/2010
/* MIMD version 7 */
/* fermion force optimized for the Asqtad action 
 * formerly fermion_force_asqtad3.c
 * Uses restart-gathers and a bit more memory for better performance
 * D.T. 1/28/98, starting from gauge_stuff.c
 * K.O. 3/99 Added optimized fattening for Asq actions
 * D.T. 4/99 Combine force calculations for both mass quarks
 * K.O. 4/99 Optimized force for Asq action
 * S.G. 7/01, modified to use t_longlink and t_fatlink
 * C.D. 10/02, consolidated quark_stuff.c and quark_stuff_tmp.c
 *
 * J.O. 3/04 Rearranged loops for optimization
 * J.O. C.D. 3/04 Copied forward links for optimization and 
 *                kept mtags open for restart_gather
 *                Worked with pointers where possible to avoid copying.
 * C.D. 6/04 Corrected memory leak
 * C.D. 3/05 Separated from quark_stuff4.c
 * D.T. 12/05 First try at RHMC version
 * C.D. 5/07  Renamed from fermion_force_asqtad3.c
 *            Collected all asqtad: single, double, multi in this file.
 *
 * In this directory, assume all paths connect even to odd sites, etc.
 * Tabulate "backwards" paths (e.g. "XDOWN" is backward path to "XUP")
 * as separate parity transforms of the fundamental paths.  They will
 * generally need a negative sign in Dslash.  See bottom for a long
 * comment on sign conventions.
 */

/* External entry points in this file
   
    eo_fermion_force_oneterm
    eo_fermion_force_twoterms
    fermion_force_asqtad_block
    fermion_force_asqtad_multi
   
 */

/* Compile with fermion_force_multi.c */

/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED - C. DeTar
 * Fermion force: 253935 for eo_fermion_force_oneterm()
 * Fermion force: 433968 for eo_fermion_force_twoterms()
 */

#include "generic_ks_includes_u1.h"	/* definitions files and prototypes */
#include <string.h>

/* This routine is valid only for Asqtad, so requires the FN flag */
#ifndef FN
BOMB THE COMPILE
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

/* Forward declarations */

static void 
u_shift_fermion_u1(complex *src, complex *dest, int dir ) ;

static void 
add_force_to_mom_u1(complex *back, complex *forw, int dir,  Real coef);

static void 
side_link_force_u1(int mu, int nu, Real coeff, complex *Path,
		complex *Path_nu, complex *Path_mu, 
		complex *Path_numu) ;

static void
insert_vFermi(complex *src, complex *dest);

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

/**********************************************************************/
/*   Version for a single set of degenerate flavors                   */
/**********************************************************************/

/* Optimized force code for the Asq and Asqtad actions                 *
 * I assume that path 0 is the one link path 2 the 3-staple            *
 * path 3 the 5-staple path 4 the 7-staple and path 5 the Lapage term. *
 * Path 1 is the Naik term.                                            */
#define Pmu          tempvec[0] 
#define Pnumu        tempvec[1]
#define Prhonumu     tempvec[2]
#define P7           tempvec[3]
#define P7rho        tempvec[4]              
#define P7rhonu      tempvec[5]
#define P5           tempvec[6]
#define P3           tempvec[7]
#define P5nu         tempvec[3]
#define P3mu         tempvec[3]
#define Popmu        tempvec[4]
#define Pmumumu      tempvec[4]

static void 
eo_fermion_force_oneterm_field_u1( Real eps, Real weight, complex *temp_x,
				int prec, ferm_links_u1_t *fn, 
				ks_action_paths *ap )
{
  /* prec is ignored for now */
  /* note CG_solution and Dslash * solution are combined in "x_off" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need x_off transported from both ends of path. */
  /* For example weight = nflavors/4 */
  register int i ;
  register site *s;
  int mu,nu,rho,sig,j;
  int DirectLinks[8] ;
  Real ferm_epsilon, coeff;
  Real OneLink, Lepage, Naik, FiveSt, ThreeSt, SevenSt ;
  complex *tempvec[8] ;
  complex *temp_mul;
  complex temp;
  Real *act_path_coeff;
  
#ifdef FFTIME
  int nflop = 253935;
  double dtime;

  dtime=-dclock();
#endif
  ferm_epsilon = 2.0*weight*eps;
  
  /*  node0_printf("STARTING fn_fermion_force_oneterm() nterms = 1\n");*/

  /* Load path coefficients from table */
  act_path_coeff = ap->act_path_coeff;

  /* Initialize the DirectLink flags */
  for(mu=0;mu<8;mu++)
    DirectLinks[mu] = 0 ;

  /* Allocate temporary vectors */
  for(mu=0;mu<8;mu++)
    tempvec[mu] = (complex *)malloc( sites_on_node*sizeof(complex) );
  temp_mul = (complex *)malloc( sites_on_node*sizeof(complex) );

  for(sig=0;sig<8;sig++)
    {
      /* Path coefficients times fermion epsilon */
      /*OneLink = act_path_coeff[0]*ferm_epsilon ; 
      Naik    = act_path_coeff[1]*ferm_epsilon ;
      ThreeSt = act_path_coeff[2]*ferm_epsilon ;
      FiveSt  = act_path_coeff[3]*ferm_epsilon ;
      SevenSt = act_path_coeff[4]*ferm_epsilon ;
      Lepage  = act_path_coeff[5]*ferm_epsilon ;*/
      /* *************************************** */
	
      //multiply in appropriate tadpole factor for direction sig
      /*if(sig==TUP || sig==TDOWN) {
	ThreeSt /= u0_t;
	FiveSt  /= u0_t;
	SevenSt /= u0_t;
	Lepage  /= u0_t;
      }
      else {
	ThreeSt /= u0_s;
	FiveSt  /= u0_s;
	SevenSt /= u0_s;
	Lepage  /= u0_s;
	}*/

      for(mu=0;mu<8;mu++)if((mu!=sig)&&(mu!=OPP_DIR(sig)))
	{
	  //multiply in appropriate tadpole factor for direction mu
	  if(mu==TUP || mu==TDOWN) {
	    Naik    = act_path_coeff[1]*ferm_epsilon/(u0_t*u0_t*u0_t);
	    OneLink = act_path_coeff[0]*ferm_epsilon/u0_t;
	    /* ThreeSt /= (u0_t*u0_t);
	       FiveSt  /= (u0_t*u0_t);
	       SevenSt /= (u0_t*u0_t); */
	    Lepage  = act_path_coeff[5]*ferm_epsilon/(u0_t*u0_t*u0_t*u0_t*u0_s);
	    ThreeSt = act_path_coeff[2]*ferm_epsilon/(u0_t*u0_t*u0_s);
	  }
	  else {
	    Naik    = act_path_coeff[1]*ferm_epsilon/(u0_s*u0_s*u0_s);
	    OneLink = act_path_coeff[0]*ferm_epsilon/u0_s;
	    /*ThreeSt /= (u0_s*u0_s);
	      FiveSt  /= (u0_s*u0_s);
	      SevenSt /= (u0_s*u0_s);*/
	    if(sig==TUP || sig==TDOWN) {
	      Lepage  = act_path_coeff[5]*ferm_epsilon/(u0_s*u0_s*u0_s*u0_s*u0_t);
	      ThreeSt = act_path_coeff[2]*ferm_epsilon/(u0_s*u0_s*u0_t);
	    }
	    else {
	      Lepage  = act_path_coeff[5]*ferm_epsilon/(u0_s*u0_s*u0_s*u0_s*u0_s);
	      ThreeSt = act_path_coeff[2]*ferm_epsilon/(u0_s*u0_s*u0_s);
	    }
	  }
	  u_shift_fermion_u1(temp_x, Pmu, OPP_DIR(mu));
	  u_shift_fermion_u1(Pmu, P3, sig);
	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu]:         x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 1 link in the path: - (numbering starts form 0) */
	      if( sig==XUP || sig==YUP ) {
		//CMULREAL(*P3, v_Fermi, temp_mul);
		insert_vFermi(P3, temp_mul);
		add_force_to_mom_u1(temp_mul, Pmu, sig, -ThreeSt) ; }
	      else {
#ifndef FOUR_DIM
		if(sig != ZUP && sig != ZDOWN)
#endif
		  add_force_to_mom_u1(P3, Pmu, sig, -ThreeSt) ; }
	    }
	  for(nu=0;nu<8;nu++)if((nu!=mu )&&(nu!=OPP_DIR(mu ))&&
				(nu!=sig)&&(nu!=OPP_DIR(sig)))
	    {
	      //multiply in appropriate tadpole factor for direction nu
	      if(nu==TUP || nu==TDOWN || mu==TUP || mu==TDOWN) {
		FiveSt  = act_path_coeff[3]*ferm_epsilon/(u0_t*u0_t*u0_s*u0_s*u0_s);
	      }
	      else {
		if(sig==TUP || sig==TDOWN)
		  FiveSt  = act_path_coeff[3]*ferm_epsilon/(u0_s*u0_s*u0_s*u0_s*u0_t);
		else
		  FiveSt  = act_path_coeff[3]*ferm_epsilon/(u0_s*u0_s*u0_s*u0_s*u0_s);
	      }
	     
	      u_shift_fermion_u1(Pmu, Pnumu, OPP_DIR(nu));
	      u_shift_fermion_u1(Pnumu, P5, sig);
	      if(GOES_FORWARDS(sig))
		{
		  /* Add the force F_sig[x+mu+nu]:      x--+             *
		   *                                   |   |             *
		   *                                   o   o             *
		   * the 2 link in the path: + (numbering starts form 0) */
		  if( sig==XUP || sig==YUP ) {
		    //CMULREAL(*P5, v_Fermi, temp_mul);
		    insert_vFermi(P5, temp_mul);
		    add_force_to_mom_u1(temp_mul, Pnumu, sig, FiveSt); }
		  else {
#ifndef FOUR_DIM
		    if(sig != ZUP && sig != ZDOWN)
#endif
		      add_force_to_mom_u1(P5, Pnumu, sig, FiveSt); }
		}
	      for(rho=0;rho<8;rho++)if((rho!=mu )&&(rho!=OPP_DIR(mu ))&&
				       (rho!=nu )&&(rho!=OPP_DIR(nu ))&&
				       (rho!=sig)&&(rho!=OPP_DIR(sig)))
		{
		  //multiply in appropriate tadpole factor for direction rho
		  if(rho==TUP || rho==TDOWN || mu==TUP || mu==TDOWN || nu==TUP || nu==TDOWN) 
		    SevenSt = act_path_coeff[4]*ferm_epsilon/(u0_t*u0_t*u0_s*u0_s*u0_s*u0_s*u0_s);
		  else {
		    SevenSt = act_path_coeff[4]*ferm_epsilon/(u0_s*u0_s*u0_s*u0_s*u0_s*u0_s*u0_t);
		  }	
		  u_shift_fermion_u1(Pnumu, Prhonumu, OPP_DIR(rho));
		  /* Length 7 paths */
		  u_shift_fermion_u1(Prhonumu, P7,sig);
		  if(GOES_FORWARDS(sig))
		    {
		      /* Add the force F_sig[x+mu+nu+rho]:  x--+             *
		       *                                   |   |             *
		       *                                   o   o             *
		       * the 3 link in the path: - (numbering starts form 0) */
		      if( sig==XUP || sig==YUP ) {
			//CMULREAL(*P7, v_Fermi, temp_mul);
			insert_vFermi(P7, temp_mul);
			add_force_to_mom_u1(temp_mul, Prhonumu, sig, -SevenSt ) ; }
		      else {
#ifndef FOUR_DIM
			if(sig != ZUP && sig != ZDOWN)
#endif
			  add_force_to_mom_u1(P7, Prhonumu, sig, -SevenSt ) ; }
		      
		    }
		  /*Add the force F_rho the 2(4) link in the path: +     */
		  u_shift_fermion_u1(P7, P7rho, rho);
		  side_link_force_u1(rho,sig,SevenSt, Pnumu, P7, Prhonumu, P7rho);
		  /* Add the P7rho vector to P5 */
		  if(FiveSt != 0) {
		    if(rho==TUP || rho==TDOWN)
		      coeff = (u0_t*u0_t*SevenSt)/(FiveSt*u0_s*u0_s); 
		    else
		       coeff = SevenSt/FiveSt;
		  }
		  else coeff = 0;
		  //if(FiveSt != 0)coeff = act_path_coeff[4]/act_path_coeff[3]; else coeff = 0;
		  FORALLSITES(i,s) {
		    CMULREAL(P7rho[i], coeff, temp);
		    CADD(temp, P5[i], P5[i]);
		  }
		}/* rho */
	      /* Length 5 paths */
	      /*Add the force F_nu the 1(3) link in the path: -     */
	      u_shift_fermion_u1(P5,P5nu, nu);
	      side_link_force_u1(nu,sig,-FiveSt,Pmu,P5, 
			      Pnumu,P5nu) ;
	      /* Add the P5nu vector to P3 */
	      if(ThreeSt != 0) {
		if(nu==TUP || nu==TDOWN) {
		  coeff = (u0_t*u0_t*FiveSt)/(ThreeSt*u0_s*u0_s);
		}
		else
		  coeff = FiveSt/ThreeSt ; 
	      }
	      else coeff = 0;
	      //if(ThreeSt != 0)coeff = act_path_coeff[3]/act_path_coeff[2]; else coeff = 0;
	      FORALLSITES(i,s) {
		CMULREAL(P5nu[i], coeff, temp);
		CADD(temp, P3[i], P3[i]);
	      }
	    }/* nu */

	  /* Now the Lepage term... It is the same with 5-link paths with
             nu=mu and FiveSt=Lepage. So Pnumu is really Pmumu */
	  u_shift_fermion_u1(Pmu, Pnumu, OPP_DIR(mu));
	  u_shift_fermion_u1(Pnumu, P5, sig);
	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu+nu]:      x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 2 link in the path: + (numbering starts form 0) */
	      if( sig==XUP || sig==YUP ) {
		//CMULREAL(*P5, v_Fermi, temp_mul);
		insert_vFermi(P5, temp_mul);
		add_force_to_mom_u1(temp_mul, Pnumu, sig, Lepage) ;  }
	      else {
#ifndef FOUR_DIM
		if(sig != ZUP && sig != ZDOWN)
#endif
		  add_force_to_mom_u1(P5, Pnumu, sig, Lepage) ; }
	    }
	  /*Add the force F_nu the 1(3) link in the path: -     */
	  u_shift_fermion_u1(P5,P5nu, mu);
	  side_link_force_u1(mu, sig, -Lepage, Pmu, P5, Pnumu, P5nu) ;
	  /* Add the P5nu vector to P3 */
	  if(ThreeSt != 0) coeff = Lepage/ThreeSt ; else coeff = 0;
	  //if(ThreeSt != 0) coeff = act_path_coeff[5]/act_path_coeff[2]; else coeff = 0;
	  FORALLSITES(i,s) {
	    CMULREAL(P5nu[i], coeff, temp);
	    CADD(temp, P3[i], P3[i]);
	  }
	  /* Length 3 paths (Not the Naik term) */
	  /*Add the force F_mu the 0(2) link in the path: +     */
	  if(GOES_FORWARDS(mu)) 
	    u_shift_fermion_u1(P3,P3mu, mu );
	  /* The above shift is not needed if mu is backwards */
	  side_link_force_u1(mu, sig, ThreeSt, temp_x, P3, Pmu, P3mu);

	  /* Finally the OneLink and the Naik term */
	  /* Check if this direction is not already done */
	  if( (!DirectLinks[mu]) ){
	    DirectLinks[mu]=1 ;
	    if(GOES_BACKWARDS(mu))/* Do only the forward terms in the Dslash */
	      {
		/* Because I have shifted with OPP_DIR(mu) Pmu is a forward *
		 * shift.                                                   */
		/* The one link */
		if( mu==XDOWN || mu==YDOWN ) { //graphene action (v_F=1 for QED)
		  //CMULREAL(*Pmu, v_Fermi, temp_mul);
		  insert_vFermi(Pmu, temp_mul);
		  add_force_to_mom_u1(temp_mul, temp_x, OPP_DIR(mu), OneLink) ; }
		else {
#ifndef FOUR_DIM
		  if(OPP_DIR(mu) != ZUP && OPP_DIR(mu) != ZDOWN)
#endif
		    add_force_to_mom_u1(Pmu, temp_x, OPP_DIR(mu), OneLink); }
		/* For the same reason Pnumu is the forward double link */

		/* Popmu is a backward shift */
		u_shift_fermion_u1(temp_x, Popmu, mu);
		/* The Naik */
		/* link no 1: - */
		if( mu==XDOWN || mu==YDOWN ) {
		  //CMULREAL(*Popmu, v_Fermi, temp_mul);
		  insert_vFermi(Popmu, temp_mul);
		  add_force_to_mom_u1(Pnumu, temp_mul, OPP_DIR(mu), -Naik) ; }
		else {
#ifndef FOUR_DIM
		  if(OPP_DIR(mu) != ZUP && OPP_DIR(mu) != ZDOWN)
#endif
		    add_force_to_mom_u1(Pnumu, Popmu, OPP_DIR(mu), -Naik); }
		/*Pmumumu can overwrite Popmu which is no longer needed */
		u_shift_fermion_u1(Pnumu, Pmumumu, OPP_DIR(mu));
		/* link no 0: + */
		if( mu==XDOWN || mu==YDOWN ) {
		  //CMULREAL(*Pmumumu, v_Fermi, temp_mul);
		  insert_vFermi(Pmumumu, temp_mul);
		  add_force_to_mom_u1(temp_mul, temp_x, OPP_DIR(mu), Naik); }
		else {
#ifndef FOUR_DIM
		  if(OPP_DIR(mu) != ZUP && OPP_DIR(mu) != ZDOWN)
#endif
		    add_force_to_mom_u1(Pmumumu, temp_x, OPP_DIR(mu), Naik); }
	      }
	    else /* The rest of the Naik terms */
	      {
		u_shift_fermion_u1(temp_x, Popmu, mu);
		/* link no 2: + */
		/* Pnumu is double backward shift */
		if( mu==XUP || mu==YUP ) {
		  //CMULREAL(*Popmu, v_Fermi, temp_mul);
		  insert_vFermi(Popmu, temp_mul);
		  add_force_to_mom_u1(temp_mul, Pnumu, mu, Naik) ; }
		else {
#ifndef FOUR_DIM
		  if(mu != ZUP && mu != ZDOWN)
#endif
		    add_force_to_mom_u1(Popmu, Pnumu, mu, Naik) ; }
	      }
	  }
	}/* mu */
      /* Here we have to do together the Naik term and the one link term */
    }/*sig */

  /* Free temporary vectors */
  for(mu=0;mu<8;mu++)
    free(tempvec[mu]) ;
  free(temp_mul);

#ifdef FFTIME
  dtime += dclock();
node0_printf("FFTIME:  time = %e (asqtad3) terms = 1 mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
/**printf("TLENGTH: %d\n",tlength);**/
#endif
} /* eo_fermion_force_oneterm (version 7) */
#undef Pmu          
#undef Pnumu        
#undef Prhonumu     
#undef P7           
#undef P7rho        
#undef P7rhonu      
#undef P5           
#undef P3           
#undef P5nu         
#undef P3mu         
#undef Popmu        
#undef Pmumumu      

void 
eo_fermion_force_oneterm_u1( Real eps, Real weight, field_offset x_off,
			  int prec, ferm_links_u1_t *fn, ks_action_paths *ap )
{
  int i ;
  site *s;
  complex *temp_x ;

  /*copy x_off to a temporary vector */
  temp_x = (complex *)malloc( sites_on_node*sizeof(complex) );
  if(temp_x == NULL){
    printf("eo_fermion_force_oneterm: No room for temporary\n");
    terminate(1);
  }

  FORALLSITES(i,s) temp_x[i] = *(complex *)F_PT(s,x_off) ;

  eo_fermion_force_oneterm_field_u1(eps, weight, temp_x, prec, fn, ap );

  free(temp_x);
}




/*   Covariant shift of the src fermion field in the direction dir  *
 *  by one unit. The result is stored in dest.                       */ 
static void 
u_shift_fermion_u1(complex *src, complex *dest, int dir ) {
  complex *tmpvec ; 
  msg_tag *mtag ;
  register site *s ;
  register int i ;
  
  if(GOES_FORWARDS(dir)) /* forward shift */
    {
      mtag = start_gather_field(src, sizeof(complex), 
				    dir, EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i,s)
	CMUL( s->link[dir], *(complex *)(gen_pt[0][i]), dest[i] );
      cleanup_gather(mtag);
    }
  else /* backward shift */
    {
      tmpvec = (complex *)malloc( sites_on_node*sizeof(complex) );
      FORALLSITES(i,s)
	CMULJ_( s->link[OPP_DIR(dir)], src[i], tmpvec[i] );
      mtag = start_gather_field(tmpvec, sizeof(complex), dir, 
				EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      /* copy the gen_pt to the dest */
      FORALLSITES(i,s)
	dest[i] = *(complex *)gen_pt[0][i];
      cleanup_gather(mtag);
      free(tmpvec) ;
    }
}

static void
insert_vFermi(complex *src, complex *dest) {

  register site *s ;
  register int i ;

  FORALLSITES(i,s)
    CMULREAL(src[i], v_Fermi, dest[i]);
  
}// insert_vFermi()

/* Add in contribution to the force */
/* Put antihermitian traceless part into momentum */

static void 
add_force_to_mom_u1(complex *back,complex *forw, int dir,Real coeff) {

  register site *s ;
  register int i ;  
  register Real tmp_coeff, temp;

  complex tmat, tmat2, temp_mul;

  //THIS MAY BE WRONG: FERMION FORCE IN +-z FOR SIDE LINKS ONLY C.W. 7/11/11
  //SPATIAL LINKS REMAIN UNITY!!!!
  /*#ifndef FOUR_DIM
  if( dir != TUP && dir != TDOWN ) 
    return;
    #endif*/
  
  if(GOES_BACKWARDS(dir))
    {
      //printf("dir = TDOWN\n");
      dir = OPP_DIR(dir) ; 
      coeff = -coeff ;
    }
  FORALLSITES(i,s){
    if(s->parity==ODD) 
      tmp_coeff = -coeff ;
    else
      tmp_coeff = coeff ;
    /*uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
      unnecessary for U(1) where momentum is already a complex number */
    //su3_projector(&(back[i]), &(forw[i]), &tmat);
    CMUL_J( back[i], forw[i], tmat );
    //scalar_mult_add_su3_matrix(&tmat2, &tmat,  tmp_coeff, &tmat2 );
    CMULREAL( tmat, tmp_coeff, temp_mul );
    /*
 #ifdef EXT_FIELD
    if( dir == YUP ) {
      if( coeff > 0 ) {
	if( coeff/(2.*epsilon) == act_path_coeff[1] ) //long link
	  CMUL( temp_mul, s->ext_link_lng[dir], temp_mul );
	else //fat link
	  CMUL( temp_mul, s->ext_link_fat[dir], temp_mul);
      }
      else {
	if( coeff/(2.*epsilon) == act_path_coeff[1] ) //long link
	  CMULJ_( temp_mul, s->ext_link_lng[dir], temp_mul );
	else //fat link
	  CMULJ_( temp_mul, s->ext_link_fat[dir], temp_mul);
      }
    }
    else if( dir == XUP ) {
      if( coeff > 0 ) {
      }
    }
    #endif */
    CADD( s->mom[dir], temp_mul, tmat2 );
    //make_anti_hermitian( &tmat2, &(s->mom[dir]) ); MAKE PURELY IMAGINARY!
    s->mom[dir].real = 0.0;
    s->mom[dir].imag = tmat2.imag;
  }
}//updated 8/4/2010




/*  This routine is needed in order to add the force on the side link *
 * of the paths in the Asq and Asqtad actions. It gets as inputs the  *
 * direction mu of the side link and the direction nu of the Dslash   *
 * term we are dealing with. Then it takes also 4 fermion fields:     *
 * Path: the piece of the path with no hop in the nu or mu direction  *
 * Path_nu: the piece of the path with a hop in nu  but not in mu     *
 * Path_mu: is Path times the link mu                                 *
 * Path_numu: is Path_nu times the link mu                            */

static void 
side_link_force_u1(int mu, int nu, Real coeff, complex *Path, 
		complex *Path_nu, complex *Path_mu, 
		complex *Path_numu) {

  complex *temp_mul;
  //graphene cannot have fat or long links in the +-z direction
#ifndef FOUR_DIM
  if(nu == ZUP || nu == ZDOWN)
    return;
#endif
  temp_mul = (complex *)malloc( sites_on_node*sizeof(complex) );
  if(temp_mul==NULL){
    printf("side_link_force_u1(%d): no room for temp_mul\n", this_node);
    terminate(1);
  }
  if(GOES_FORWARDS(mu))
    {
      /*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
      if(GOES_FORWARDS(nu)) {
	if( nu == XUP || nu == YUP )  {//graphene action (v_Fermi = 1 for QED)
	  //CMULREAL(*Path_numu, v_Fermi, temp_mul);
	  insert_vFermi(Path_numu, temp_mul);
	  add_force_to_mom_u1(temp_mul, Path, mu, coeff ) ; }
	else
	  add_force_to_mom_u1(Path_numu, Path, mu, coeff );
      }
      else {
	if( nu == XDOWN || nu == YDOWN ) {
	  //CMULREAL(*Path_numu, v_Fermi, temp_mul);
	  insert_vFermi(Path_numu, temp_mul);
	  add_force_to_mom_u1(Path, temp_mul, OPP_DIR(mu), -coeff );  }
	else
	  add_force_to_mom_u1(Path, Path_numu, OPP_DIR(mu), -coeff ); /* ? extra - */
	
      }
    }
  else /*GOES_BACKWARDS(mu)*/
    {
      /* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */ 
      //printf("mu is BACKWARDS\n");
      if(GOES_FORWARDS(nu)) {
	//printf("mu is BACKWARDS and nu is FORWARDS\n");
	if( nu == XUP || nu == YUP ) {
	  //CMULREAL(*Path_nu, v_Fermi, temp_mul);
	  insert_vFermi(Path_nu, temp_mul);
	  add_force_to_mom_u1(temp_mul, Path_mu, mu, -coeff) ; }
	else
	  add_force_to_mom_u1(Path_nu, Path_mu, mu, -coeff); /* ? extra - */
      }
      else {
	if( nu == XDOWN || nu == YDOWN ) {
	  //CMULREAL(*Path_nu, v_Fermi, temp_mul);
	  insert_vFermi(Path_nu, temp_mul);
	  add_force_to_mom_u1(Path_mu, temp_mul, OPP_DIR(mu), coeff) ; }
	else
	  add_force_to_mom_u1(Path_mu, Path_nu, OPP_DIR(mu), coeff);
      }
    }
  free(temp_mul);
}


/* LONG COMMENTS
   Here we have combined "xxx", (offset "x_off")  which is
(M_adjoint M)^{-1} phi, with Dslash times this vector, which goes in the
odd sites of xxx.  Recall that phi is defined only on even sites.  In
computing the fermion force, we are looking at

< X |  d/dt ( Dslash_eo Dslash_oe ) | X >
=
< X | d/dt Dslash_eo | T > + < T | d/dt Dslash_oe | X >
where T = Dslash X.

The subsequent manipulations to get the coefficent of H, the momentum
matrix, in the simulation time derivative above look the same for
the two terms, except for a minus sign at the end, if we simply stick
T, which lives on odd sites, into the odd sites of X

 Each path in the action contributes terms when any link of the path
is the link for which we are computing the force.  We get a minus sign
for odd numbered links in the path, since they connect sites of the
opposite parity from what it would be for an even numbered link.
Minus signs from "going around" plaquette - ie KS phases, are supposed
to be already encoded in the path coefficients.
Minus signs from paths that go backwards are supposed to be already
encoded in the path coefficients.

Here, for example, are comments reproduced from the force routine for
the one-link plus Naik plus single-staple-fat-link action:

 The three link force has three contributions, where the link that
was differentiated is the first, second, or third link in the 3-link
path, respectively.  Diagramatically, where "O" represents the momentum,
the solid line the link corresponding to the momentum, and the dashed
lines the other links:
 

	O______________ x ............ x ...............
+
	x..............O______________x.................
+
	x..............x..............O________________
Think of this as
	< xxx | O | UUUxxx >		(  xxx, UUUX_p3 )
+
	< xxx U | O | UUxxx >		( X_m1U , UUX_p2 )
+
	< xxx U U | O | Uxxx >		( X_m2UU , UX_p1 )
where "U" indicates parallel transport, "X_p3" is xxx displaced
by +3, etc.
Note the second contribution has a relative minus sign
because it effectively contributes to the <odd|even>, or M_adjoint,
part of the force when we work on an even site. i.e., for M on
an even site, this three link path begins on an odd site.

The staple force has six contributions from each plane containing the
link direction:
Call these diagrams A-F:


	x...........x		O____________x
		    .			     .
		    .			     .
		    .			     .
		    .			     .
		    .			     .
	O___________x		x............x
	   (A)			    (B)



	x	    x		O____________x
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	O___________x		x	     x
	   (C)			    (D)



	x...........x		O____________x
	.			.
	.			.
	.			.
	.			.
	.			.
	O___________x		x............x
	   (E)			    (F)

As with the Naik term, diagrams C and D have a relative minus
sign because they connect sites of the other parity.

Also note an overall minus sign in the staple terms relative to the
one link term because, with the KS phase factors included, the fat
link is  "U - w3 * UUU", or the straight link MINUS w3 times the staples.

Finally, diagrams B and E get one more minus sign because the link
we are differentiating is in the opposite direction from the staple
as a whole.  You can think of this as this "U" being a correction to
a "U_adjoint", but the derivative of U is iHU and the derivative
of U_adjoint is -iHU_adjoint.

*/
/* LONG COMMENT on sign conventions
In most of the program, the KS phases and antiperiodic boundary
conditions are absorbed into the link matrices.  This greatly simplfies
multiplying by the fermion matrix.  However, it requires care in
specifying the path coefficients.  Remember that each time you
encircle a plaquette, you pick up a net minus sign from the KS phases.
Thus, when you have more than one path to the same point, you generally
have a relative minus sign for each plaquette in a surface bounded by
this path and the basic path for that displacement.

Examples:
  Fat Link:
    Positive:	X-------X

    Negative     --------
	 	|	|
		|	|
		X	X

  Naik connection, smeared
    Positive:	X-------x-------x-------X

    Negative:	---------
		|	|
		|	|
		X	x-------x-------X

    Positive:	--------x--------
		|		|
		|		|
		X		x-------X

    Negative:	--------x-------x-------x
		|			|
		|			|
		X			X
*/



/* Comment on acceptable actions.
   We construct the backwards part of dslash by reversing all the
   paths in the forwards part.  So, for example, in the p4 action
   the forwards part includes +X+Y+Y

		X
		|
		|
		X
		|
		|
	X---->--X

  so we put -X-Y-Y in the backwards part.  But this isn't the adjoint
  of U_x(0)U_y(+x)U_y(+x+y).  Since much of the code assumes that the
  backwards hop is the adjoint of the forwards (for example, in
  preventing going to 8 flavors), the code only works for actions
  where this is true.  Roughly, this means that the fat link must
  be symmetric about reflection around its midpoint.  Equivalently,
  the paths in the backwards part of Dslash are translations of the
  paths in the forwards part.  In the case of the "P4" or knight's move
  action, this means that we have to have both paths
   +X+Y+Y and +Y+Y+X to the same point, with the same coefficients.
  Alternatively, we could just use the symmetric path +Y+X+Y.
*/
