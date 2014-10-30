/**************** fermion_links_helpers_u1.c *****************************/
/* MILC Version 7 */
/* Routines used by all of the fermion_links*.c files */
//Modified 5/18/2010 -CW

#include "generic_ks_includes_u1.h"	/* definitions files and prototypes */
#define IMP_QUARK_ACTION_INFO_ONLY
#include <quark_action.h>

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

#ifdef ASQ_OPTIMIZED_FATTENING
static void 
compute_gen_staple_field_u1(complex *staple, int mu, int nu, 
			 complex *link, complex *fatlink, Real coef);

static void 
compute_gen_staple_site_u1(complex *staple, int mu, int nu, 
			field_offset link, complex* fatlink, Real coef);
#endif

#ifdef EXT_FIELD
/* set up external magnetic field in the z direction according to the 
   conventions of Bali et. al. arXiv:1111.495 */
void create_ext_field() {

  register int i, dir;
  register site *s;
  double temp;
  complex link,link2,temp1,temp2,temp3,temp4,ptemp,
    ptemp2,temp_mul,temp_mul2;
  register Real t2,t3,t4,t5,t6;

  t2 = 1/2.0;
  t3 = 1/3.0;
  t4 = 1/4.0;
  t5 = 1/5.0;
  t6 = 1/6.0;

  FORALLMYSITES(i,s){
    for (dir=XUP; dir<=TUP; dir++){
      if(dir==XUP) {
	//boundary terms in the fat and 
	if(s->x == (nx - 1)) {
	  temp = (double)s->y*nx;
	  s->ext_potential_fat[dir] = B_ext*temp; //check this sign
	  s->ext_potential_lng[dir] = B_ext*temp;
	}
	else if(s->x == (nx - 2) || s->x == (nx - 3)) {
	  temp = (double)s->y*nx;
	  s->ext_potential_fat[dir] = 0.;
	  s->ext_potential_lng[dir] = B_ext*temp; //check this sign
	}
	else {
	  s->ext_potential_fat[dir] = 0.;
	  s->ext_potential_lng[dir] = 0.;
	}
      }//xup
      else if(dir==YUP) {
	temp = (double)s->x;
	s->ext_potential_fat[dir] = B_ext*temp;
	s->ext_potential_lng[dir] = 3.*B_ext*temp;
      }//yup
      else { // ext A_mu = 0 for mu={z,t}
	s->ext_potential_fat[dir] = 0.;
	s->ext_potential_lng[dir] = 0.;
      }
      ptemp.real = ptemp2.real = 0.0;
      ptemp.imag = s->ext_potential_fat[dir];
      ptemp2.imag = s->ext_potential_lng[dir];
      link.real = link2.real = 1.0;
      link.imag = link2.imag = 0.0;
            
      CMULREAL(ptemp, t6, temp_mul);
      CMULREAL(ptemp2, t6, temp_mul2);
      CADD(link, temp_mul, temp2);
      CADD(link2, temp_mul2, temp4);
      CMUL(ptemp, temp2, temp1);
      CMUL(ptemp2, temp4, temp3);

      CMULREAL(temp1, t5, temp_mul);
      CMULREAL(temp3, t5, temp_mul2);
      CADD(temp_mul, link, temp2);
      CADD(temp_mul2, link2, temp4);
      CMUL(ptemp, temp2, temp1);
      CMUL(ptemp2, temp4, temp3);

      CMULREAL(temp1, t4, temp_mul);
      CMULREAL(temp3, t4, temp_mul2);
      CADD(temp_mul, link, temp2);
      CADD(temp_mul2, link2, temp4);
      CMUL(ptemp, temp2, temp1);
      CMUL(ptemp2, temp4, temp3);

      CMULREAL(temp1, t3, temp_mul);
      CMULREAL(temp3, t3, temp_mul2);
      CADD(temp_mul, link, temp2);
      CADD(temp_mul2, link2, temp4);
      CMUL(ptemp, temp2, temp1);
      CMUL(ptemp2, temp4, temp3);

      CMULREAL(temp1, t2, temp_mul);
      CMULREAL(temp3, t2, temp_mul2);
      CADD(temp_mul, link, temp2);
      CADD(temp_mul2, link2, temp4);
      CMUL(ptemp, temp2, temp1);
      CMUL(ptemp2, temp4, temp3);
      
      CMULREAL(temp1, 1.0, temp_mul);
      CMULREAL(temp3, 1.0, temp_mul2);
      CADD(temp_mul, link, temp2); 
      CADD(temp_mul2, link2, temp4);
      s->ext_link_fat[dir].real = temp2.real;
      s->ext_link_fat[dir].imag = temp2.imag;
      s->ext_link_lng[dir].real = temp4.real;
      s->ext_link_lng[dir].imag = temp4.imag;
      temp = (double)cabs(&(s->ext_link_fat[dir]));
      CMULREAL(s->ext_link_fat[dir], 1/temp, s->ext_link_fat[dir]);
      temp = (double)cabs(&(s->ext_link_lng[dir]));
      CMULREAL(s->ext_link_lng[dir], 1/temp, s->ext_link_lng[dir]);

    }//dir
  }//FORALLSITES

}// create_ext_field()
#endif

/* calculates the U(1) phase from the real potential */
void exp_links() {

  register int i, dir;
  register site *s;
  double temp;
  complex link,temp1,temp2,ptemp,temp_mul;
  register Real t2,t3,t4,t5,t6;

  t2 = 1/2.0;
  t3 = 1/3.0;
  t4 = 1/4.0;
  t5 = 1/5.0;
  t6 = 1/6.0;
  
  
  FORALLSITES(i,s){
    for (dir=XUP; dir<=TUP; dir++){
      ptemp.real = 0.0;
      ptemp.imag = s->potential[dir];
      //link = &(s->link[dir]);
      link.real = 1.0; link.imag = 0.0;
      
      //CMUL( ptemp, *link, temp1);
      
      CMULREAL(ptemp, t6, temp_mul);
      CADD(link, temp_mul, temp2);
      CMUL(ptemp, temp2, temp1);
      
      CMULREAL(temp1, t5, temp_mul);
      CADD(temp_mul, link, temp2);
      CMUL(ptemp, temp2, temp1);
      
      CMULREAL(temp1, t4, temp_mul);
      CADD(temp_mul, link, temp2);
      CMUL(ptemp, temp2, temp1);
      
      CMULREAL(temp1, t3, temp_mul);
      CADD(temp_mul, link, temp2);
      CMUL(ptemp, temp2, temp1);
      
      CMULREAL(temp1, t2, temp_mul);
      CADD(temp_mul, link, temp2);
      CMUL(ptemp, temp2, temp1);
      
      CMULREAL(temp1, 1.0, temp_mul);
      CADD(temp_mul, link, temp2); 
      s->link[dir].real = temp2.real;
      s->link[dir].imag = temp2.imag;
      temp = (double)cabs(&(s->link[dir]));
      CMULREAL(s->link[dir], 1/temp, s->link[dir]);
    }//dir
  }//FORALLSITES

}// exp_links()


void load_longlinks_u1(ferm_links_u1_t *fn, ks_action_paths *ap) {
  complex **t_ll = &fn->lng;
  register int i;
  register site *s;
  int ipath,dir;
  int disp[4];
  int num_q_paths = ap->num_q_paths;
  Q_path *q_paths = ap->q_paths;
  register complex *long1;
  complex *staple = NULL, *tempmat1 = NULL;
  char myname[] = "load_longlinks";
  complex temp_mul;
  Real temp_coeff;

#ifdef LLTIME
  int nflop = 1804;
  double dtime;
  dtime=-dclock();
#endif

  if( phases_in != 1){
    node0_printf("BOTCH: %s needs phases in\n",myname);
    terminate(0);
  }
  /* Allocate space for t_ll if NULL */
  if(*t_ll == NULL){
    *t_ll = (complex *)special_alloc(sites_on_node*4*sizeof(complex));
    if(*t_ll==NULL){
      printf("%s(%d): no room for t_ll\n",myname, this_node);
      terminate(1);
    }
  }
  
  staple = (complex *)special_alloc(sites_on_node*sizeof(complex));
  if(staple == NULL){
    printf("%s(%d): Can't malloc temporary\n",myname,this_node);
    terminate(1);
  }

  tempmat1 = (complex *)special_alloc(sites_on_node*sizeof(complex));
  if(tempmat1 == NULL){
    printf("%s(%d): Can't malloc temporary\n",myname,this_node);
    terminate(1);
  }

  for (dir=XUP; dir<=TUP; dir++){ /* loop over longlink directions */
    /* set longlink to zero */
    FORALLSITES(i,s){
      long1 = *t_ll + 4*i +dir;
      long1->real = 0.0; 
      long1->imag = 0.0;
    }
    /* loop over paths, checking for ones with total displacement 3*dir */
    for( ipath=0; ipath<num_q_paths; ipath++ ){  /* loop over paths */
	/* compute total displacement of path */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	for( i=0; i<q_paths[ipath].length; i++){
	  if( GOES_FORWARDS(q_paths[ipath].dir[i]) )
	    disp[        q_paths[ipath].dir[i]  ]++;
	  else
	    disp[OPP_DIR(q_paths[ipath].dir[i]) ]--;
	}
	for( disp[dir]+=3,i=XUP; i<=TUP; i++)if(disp[i]!=0)break;
	if( i<=TUP )continue;  /* skip if path doesn't go to right place */
/**printf("ipath = %d, found a path:  ",ipath);
for(j=0;j<q_paths[ipath].length;j++)printf("\t%d", q_paths[ipath].dir[j]);
printf("\n");**/

	path_product_u1( q_paths[ipath].dir, q_paths[ipath].length, tempmat1 );
	FORALLSITES(i,s){
	  CONJG( tempmat1[i], staple[i] );
	  long1 = *t_ll + 4*i + dir;
	  //tadpole factors
	  if(dir==TUP)
	    temp_coeff = q_paths[ipath].coeff/(u0_t*u0_t*u0_t);
	  else
	    temp_coeff = q_paths[ipath].coeff/(u0_s*u0_s*u0_s);
	  CMULREAL( staple[i], -temp_coeff, temp_mul);
	  CADD( *long1, temp_mul, *long1 );
#ifdef EXT_FIELD
	  //CMUL( *long1, s->ext_link_lng[dir], *long1);
#endif
	  /* minus sign in coeff. because we used backward path*/
	}
    } /* ipath */
  } /* loop over directions */

  special_free(staple); staple = NULL;
  special_free(tempmat1); tempmat1 = NULL;

#ifdef LLTIME
dtime += dclock();
node0_printf("LLTIME(long): time =  %e (Naik) mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
}  /* load_longlinks() */

/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED_FATTENING - C. DeTar
 * Fatlinks:       61632 for load_fatlinks
 */

/* KS phases and APBC must be in the links. See long comment at 
   end of fermion_force_general.c */
void load_fatlinks_u1(ferm_links_u1_t *fn, ks_action_paths *ap){
  complex **t_fl = &fn->fat;
  register int i, j;
  register site *s;
  complex temp_mul;
  int dir;
  register complex *fat1;
  complex *staple = NULL, *tempmat1 = NULL;
  char myname[] = "load_fatlinks";

#ifdef ASQ_OPTIMIZED_FATTENING
  int  nu,rho,sig ;
  Real one_link;
  Real temp;
  Real temp_coeff[6];
  Real *act_path_coeff = ap->act_path_coeff;
#ifdef LLTIME
  char method[] = "Asqtad opt";
#endif
#else
  int ipath;
  int disp[4];
  int num_q_paths = ap->num_q_paths;
  Q_path *q_paths = ap->q_paths;
#ifdef LLTIME
  char method[] = "FN nonopt";
#endif
#endif

#ifdef LLTIME
  int nflop = 61632;
double dtime;
dtime=-dclock();
#endif

  if( phases_in != 1){
    node0_printf("BOTCH: %s needs phases in\n",myname);
    terminate(1);
  }

  /* Allocate space for t_fl if NULL */
  if(*t_fl == NULL){
    *t_fl = (complex *)special_alloc(sites_on_node*4*sizeof(complex));
    if(*t_fl==NULL){
      printf("NODE %d: no room for t_fl\n",this_node);
      terminate(1);
    }
  }
  
  staple = (complex *)special_alloc(sites_on_node*sizeof(complex));
  if(staple == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

  tempmat1 = (complex *)special_alloc(sites_on_node*sizeof(complex));
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

#ifndef  ASQ_OPTIMIZED_FATTENING   /* general case code */
  //rephase(OFF);
  /*for (dir=XUP; dir<=TUP; dir++){
    FORALLSITES(i,s){
      node0_printf("site %d dir %d %e, %e\n", i, dir, 
		   lattice[i].link[dir].real, 
		   lattice[i].link[dir].imag);
    }
    }*/
  //rephase(ON);
  for (dir=XUP; dir<=TUP; dir++){ /* loop over fatlink directions */
    /* set fatlink to zero */
    FORALLSITES(i,s){
      fat1 = (*t_fl) + 4*i + dir;
      fat1->real = 0.0; 
      fat1->imag = 0.0;
    }
    
    /* loop over paths, checking for ones with total displacement 1*dir */
    for( ipath=0; ipath<num_q_paths; ipath++ ){  /* loop over paths */
	/* compute total displacement of path */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	for( i=0; i<q_paths[ipath].length; i++){
	  if( GOES_FORWARDS(q_paths[ipath].dir[i]) )
	    disp[        q_paths[ipath].dir[i]  ]++;
	  else
	    disp[OPP_DIR(q_paths[ipath].dir[i]) ]--;
	}
	for( disp[dir]+=1,i=XUP; i<=TUP; i++)if(disp[i]!=0)break;
	if( i<=TUP )continue;  /* skip if path doesn't go to right place */
/**printf("dir = %d, found a path:  ",dir);
for(j=0;j<q_paths.[ipath].length;j++)printf("\t%d", q_paths[ipath].dir[j]);
printf("\n");**/

	path_product_u1( q_paths[ipath].dir, q_paths[ipath].length, tempmat1 );
	//DEBUGGING TADPOLE FACTORS FOR ANISOTROPIC LATTICE 9/23/2011 C.W.
	/*node0_printf("PATH #%d of length %d with coeff %e\n", ipath, 
	  q_paths[ipath].length, q_paths[ipath].coeff);
	for(j=0; j<q_paths[ipath].length;j++) {
	  if(q_paths[ipath].dir[j]==TUP || q_paths[ipath].dir[j]==TDOWN) {
	    q_paths[ipath].coeff /= u0_t;
	    //node0_printf("dir[%d] = %d (TIME)\n", j, q_paths[ipath].dir[j]);
	  }
	  else {
	    q_paths[ipath].coeff /= u0_s;
	    //node0_printf("dir[%d] = %d (SPACE)\n", j, q_paths[ipath].dir[j]);
	  }
	}
	//node0_printf("TADPOLE_IMPROVED coeff %e\n", q_paths[ipath].coeff);
	*/
	FORALLSITES(i,s){
	  CONJG( tempmat1[i], staple[i] );
	  fat1 = (*t_fl) +  4*i + dir;
	  CMULREAL( staple[i], -q_paths[ipath].coeff, temp_mul);
	  CADD( *fat1, temp_mul, *fat1 );
	  /* minus sign in coeff. because we used backward path*/
	}
    } /* ipath */
  } /* loop over directions */
#else	/* ASQ_OPTIMIZED_FATTENING, for Asq and Asqtad actions */
/*  Optimized fattening code for the Asq and Asqtad actions.           *
 * I assume that path 0 is the one link path 2 the 3-staple            *
 * path 3 the 5-staple path 4 the 7-staple and path 5 the Lapage term. *
 * Path 1 is the Naik term.                                            */
 
 /* to fix up the Lepage term, included by a trick below */
 //one_link = (act_path_coeff[0] - 6.0*act_path_coeff[5]);
  //node0_printf("ASQ_OPTIMIZED_FATTENING ROUTINE!!!!\n");
  one_link = act_path_coeff[0];
  
  for (dir=XUP; dir<=TUP; dir++){
    //tadpole factor for one-link
    if(dir == TUP)
      temp_coeff[0] = one_link/u0_t;
    else
      temp_coeff[0] = one_link/u0_s;
    FORALLSITES(i,s) /* Intialize fat links with c_1*U_\mu(x) */
      {
	fat1 = (*t_fl) +  4*i + dir;
	CMULREAL( s->link[dir], temp_coeff[0], *fat1 );
      }
    for(nu=XUP; nu<=TUP; nu++) if(nu!=dir) {
	//tapdole factors for simple staple amd Lepage term
	/* also add in appropriate tadpole improved coeff. to fix up the Lepage term */
	if(nu==TUP) {
	  temp_coeff[2] = act_path_coeff[2]/(u0_t*u0_t*u0_s);
	  temp_coeff[5] = act_path_coeff[5]/(u0_t*u0_t*u0_t*u0_t*u0_s);
	  temp = -2.*act_path_coeff[5]/(u0_t*u0_t*u0_t*u0_t*u0_s);
	}
	else {
	  if(dir==TUP) {
	    temp_coeff[2] = act_path_coeff[2]/(u0_s*u0_s*u0_t);
	    temp_coeff[5] = act_path_coeff[5]/(u0_s*u0_s*u0_s*u0_s*u0_t);
	    temp = -2.*act_path_coeff[5]/(u0_s*u0_s*u0_s*u0_s*u0_t);
	  }
	  else {
	    temp_coeff[2] = act_path_coeff[2]/(u0_s*u0_s*u0_s);
	    temp_coeff[5] = act_path_coeff[5]/(u0_s*u0_s*u0_s*u0_s*u0_s);
	    temp = -2.*act_path_coeff[5]/(u0_s*u0_s*u0_s*u0_s*u0_s);
	  }
	}
	
	FORALLSITES(i,s) {
	  fat1 = (*t_fl) +  4*i + dir;
	  CMULREAL( s->link[dir], temp, temp_mul );
	  CADD(*fat1, temp_mul, *fat1);
	}
	compute_gen_staple_site_u1(staple,dir,nu,F_OFFSET(link[dir]),
				   *t_fl, temp_coeff[2]);
	/* The Lepage term */
	/* Note this also involves modifying c_1 (above) */
	compute_gen_staple_field_u1(NULL,dir,nu,staple,
				    *t_fl, temp_coeff[5]);
	
	for(rho=XUP; rho<=TUP; rho++) if((rho!=dir)&&(rho!=nu)) {
	    //tadpole factors for 5-link staple 
	    if(rho == TUP || nu == TUP) {
	      temp_coeff[3] = act_path_coeff[3]/(u0_t*u0_t*u0_s*u0_s*u0_s);
	    }
	    else {
	      if(dir == TUP)
		temp_coeff[3] = act_path_coeff[3]/(u0_s*u0_s*u0_s*u0_s*u0_t);
	      else
		temp_coeff[3] = act_path_coeff[3]/(u0_s*u0_s*u0_s*u0_s*u0_s);
	    }
	    
	    compute_gen_staple_field_u1( tempmat1, dir, rho, staple,
					 *t_fl, temp_coeff[3]);
	    for(sig=XUP; sig<=TUP; sig++)
	      if((sig!=dir)&&(sig!=nu)&&(sig!=rho))
		{
		  //tadpole factors for 7-link staple
		 if(sig==TUP || rho==TUP || nu==TUP) 
		   temp_coeff[4] = act_path_coeff[4]/(u0_t*u0_t*u0_s*u0_s*u0_s*u0_s*u0_s);
		 else {
		   if(dir==TUP)
		     temp_coeff[4] = act_path_coeff[4]/(u0_s*u0_s*u0_s*u0_s*u0_s*u0_s*u0_t);
		   else
		     temp_coeff[4] = act_path_coeff[4]/(u0_s*u0_s*u0_s*u0_s*u0_s*u0_s*u0_s);
		 }
		 
		 compute_gen_staple_field_u1(NULL,dir,sig,tempmat1,
					  *t_fl, temp_coeff[4]);

		 
		} /* sig */
	  } /* rho */
      } /* nu */
  }/* dir */  
#endif

#ifdef EXT_FIELD
/*for(dir=XUP; dir<=TUP; dir++){
   FORALLSITES(i,s) {
     fat1 = (*t_fl) +  4*i + dir;
     CMUL(*fat1, s->ext_link_fat[dir], *fat1);
     node0_printf("fat link at  site[%d] in dir=%d: %e, %e\n", i, dir,
       fat1->real, fat1->imag); 
   }
   } */
#endif
  /*for(dir=XUP; dir<=TUP; dir++){
    FORALLSITES(i,s) {
      fat1 = (*t_fl) +  4*i + dir;
      node0_printf("fat link at  site[%d] in dir=%d: %e, %e\n", i, dir,
		  fat1->real, fat1->imag); 
    }
  }
  exit(1); //DEBUG LINK FATTENING
*/
  special_free(staple);  staple = NULL;
  special_free(tempmat1); tempmat1 = NULL;
#ifdef LLTIME
 dtime += dclock();
 node0_printf("LLTIME(Fat): time = %e (%s) mflops = %e\n",dtime,method,
	      (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
}  /* load_fatlinks() */ 


#ifdef  ASQ_OPTIMIZED_FATTENING   /* Asqtad action only, "_fn" executables */

/*#ifndef FN
BOMB THE COMPILE
#endif*/

static void 
compute_gen_staple_site_u1(complex *staple, int mu, int nu, 
			field_offset link, complex* fatlink, Real coef) {
  complex tmat1,tmat2;
  msg_tag *mtag0,*mtag1;
  complex *tempmat = NULL;
  complex temp_mul;
  register site *s ;
  register int i ;
  register complex *fat1;

  /* Computes the staple :
                   mu
                +-------+
          nu	|	|
		|	|
		X	X
    Where the mu link can be any complex. The result is saved in staple.
    if staple==NULL then the result is not saved.
    It also adds the computed staple to the fatlink[mu] with weight coef.
  */

  /* Upper staple */
  mtag0 = start_gather_site( link, sizeof(complex), nu, EVENANDODD, gen_pt[0] );
  mtag1 = start_gather_site( F_OFFSET(link[nu]), sizeof(complex), mu, 
			EVENANDODD, gen_pt[1] );
  wait_gather(mtag0);
  wait_gather(mtag1);
  
  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      CMUL_J( *(complex *)gen_pt[0][i], *(complex *)gen_pt[1][i], tmat1 );
      CMUL( s->link[nu], tmat1, staple[i] );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      CMUL_J( *(complex *)gen_pt[0][i],  *(complex *)gen_pt[1][i], tmat1);
      CMUL( s->link[nu], tmat1, tmat2 );
      fat1 = &(fatlink[4*i+mu]);
      CMULREAL( tmat2, coef, temp_mul );
      CADD( *fat1, temp_mul, *fat1 );
    }
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);

  /* lower staple */
  tempmat = (complex *)special_alloc( sites_on_node*sizeof(complex) );
  mtag0 = start_gather_site( F_OFFSET(link[nu]),
			sizeof(complex), mu, EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);
  FORALLSITES(i,s){
    CMULJ_( s->link[nu],*(complex *)F_PT(s, link), tmat1 );
    CMUL( tmat1,  *(complex *)gen_pt[0][i], tempmat[i] );
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_field( tempmat, sizeof(complex),
				  OPP_DIR(nu), EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);

  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      CADD( staple[i],  *(complex *)gen_pt[0][i], staple[i] );
      fat1 = &(fatlink[4*i+mu]);
      CMULREAL( staple[i], coef, temp_mul );
      CADD ( *fat1, temp_mul, *fat1 );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      fat1 = &(fatlink[4*i+mu]);
      CMULREAL(  *(complex *)gen_pt[0][i], coef, temp_mul );
      CADD ( *fat1, temp_mul, *fat1 );
    }
  }

  special_free(tempmat); tempmat = NULL;
  cleanup_gather(mtag0);
} /* compute_gen_staple_site_u1 */

/* Asqtad action only, "_fn" executables */

/*#ifndef FN
BOMB THE COMPILE
#endif*/

static void 
compute_gen_staple_field_u1(complex *staple, int mu, int nu, 
			 complex *link, complex *fatlink, Real coef) {
  complex tmat1,tmat2, temp_mul;
  msg_tag *mtag0,*mtag1;
  complex *tempmat = NULL;
  register site *s ;
  register int i ;
  register complex *fat1;

  /* Computes the staple :
                   mu
                +-------+
          nu	|	|
		|	|
		X	X
    Where the mu link can be any complex. The result is saved in staple.
    if staple==NULL then the result is not saved.
    It also adds the computed staple to fatlink[mu] with weight coef.
  */

  /* Upper staple */
  mtag0 = start_gather_field( link, sizeof(complex), nu, EVENANDODD, gen_pt[0] );
  mtag1 = start_gather_site( F_OFFSET(link[nu]), sizeof(complex), mu, 
			EVENANDODD, gen_pt[1] );
  wait_gather(mtag0);
  wait_gather(mtag1);
  
  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      CMUL_J( *(complex *)gen_pt[0][i], *(complex *)gen_pt[1][i], tmat1 );
      CMUL( s->link[nu], tmat1, staple[i] );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      CMUL_J( *(complex *)gen_pt[0][i], *(complex *)gen_pt[1][i], tmat1);
      CMUL( s->link[nu], tmat1, tmat2 );
      fat1 = &(fatlink[4*i+mu]);
      CMULREAL( tmat2, coef, temp_mul );
      CADD( *fat1, temp_mul, *fat1 );
    }
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);

  /* lower staple */
  tempmat = (complex *)special_alloc( sites_on_node*sizeof(complex) );
  mtag0 = start_gather_site( F_OFFSET(link[nu]),
			sizeof(complex), mu, EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);
  FORALLSITES(i,s){
    CMULJ_( s->link[nu], link[i], tmat1 );
    CMUL( tmat1, *(complex *)gen_pt[0][i], tempmat[i] );
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_field( tempmat, sizeof(complex),
				  OPP_DIR(nu), EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);

  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      CADD( staple[i], *(complex *)gen_pt[0][i], staple[i] );
      fat1 = &(fatlink[4*i+mu]);
      CMULREAL( staple[i], coef, temp_mul );
      CADD ( *fat1, temp_mul, *fat1 );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      fat1 = &(fatlink[4*i+mu]);
      CMULREAL( *(complex *)gen_pt[0][i], coef, temp_mul );
      CADD ( *fat1, temp_mul, *fat1 );
    }
  }

  special_free(tempmat); tempmat = NULL;
  cleanup_gather(mtag0);
} /* compute_gen_staple_field */
#endif  /* ASQ_OPTIMIZED_FATTENING   */

/* Move up the backward longlinks.  Result in t_lbl */
void 
load_longbacklinks_u1(ferm_links_u1_t *fn){
  complex **t_lbl = &fn->lngback;
  complex *t_ll = fn->lng;
  register int i;
  register site *s;
  int dir;
  complex *tempmat1 = NULL;
  msg_tag *tag[4];
  char myname[] = "load_longbacklinks";

  /* Allocate space for t_lbl if NULL */
  if(*t_lbl == NULL){
    *t_lbl = (complex *)special_alloc(sites_on_node*4*sizeof(complex));
    if(*t_lbl==NULL){
      printf("%s(%d): no room for t_lbl\n",myname,this_node);
      terminate(1);
    }
  }

  tempmat1 = (complex *)special_alloc(sites_on_node*sizeof(complex));
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

  /* gather backwards longlinks */
  for( dir=XUP; dir<=TUP; dir ++){
    FORALLSITES(i,s){
      tempmat1[i] = t_ll[dir+4*i];
    }
    tag[dir] = start_gather_field( tempmat1,
      sizeof(complex), OPP_3_DIR(DIR3(dir)), EVENANDODD, gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLSITES(i,s){
      CONJG( *(complex *)gen_pt[dir][i],
      (*t_lbl)[dir + 4*i] ); 
    }
    cleanup_gather( tag[dir] );
  }

  special_free(tempmat1); tempmat1 = NULL;
}

/* Move up the backward fatlinks.  Result in t_fbl */
void 
load_fatbacklinks_u1(ferm_links_u1_t *fn){
  complex **t_fbl = &fn->fatback;
  complex *t_fl = fn->fat;
  register int i;
  register site *s;
  int dir;
  complex *tempmat1 = NULL;
  msg_tag *tag[4];
  char myname[] = "load_fatbacklinks";

  /* Allocate space for t_fbl if NULL */
  if(*t_fbl == NULL){
    *t_fbl = (complex *)special_alloc(sites_on_node*4*sizeof(complex));
    if(*t_fbl==NULL){
      printf("%s(%d): no room for t_fbl\n",myname,this_node);
      terminate(1);
    }
  }

  tempmat1 = (complex *)special_alloc(sites_on_node*sizeof(complex));
  if(tempmat1 == NULL){
    printf("%s: Gan't malloc temporary\n",myname);
    terminate(1);
  }

  /* gather backwards fatlinks */
  for( dir=XUP; dir<=TUP; dir ++){
    FORALLSITES(i,s){
      tempmat1[i] = t_fl[dir+4*i];
    }
    tag[dir] = start_gather_field( tempmat1,
      sizeof(complex), OPP_DIR(dir), EVENANDODD, gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLSITES(i,s){
      CONJG( *(complex *)gen_pt[dir][i],
      (*t_fbl)[dir + 4*i] );
    }
    cleanup_gather( tag[dir] );
  }

  special_free(tempmat1); tempmat1 = NULL;
}

static void 
free_t_links(complex **t_l){
  if(*t_l != NULL) special_free(*t_l);
  *t_l = NULL;
}

/* Wrappers for MILC call to QOP */
void 
free_fn_links_u1(ferm_links_u1_t *fn){
  free_t_links(&fn->fat);
  free_t_links(&fn->lng);
#ifdef DBLSTORE_FN
  free_t_links(&fn->fatback);
  free_t_links(&fn->lngback);
#endif
  invalidate_all_ferm_links_u1(fn);
}

#ifdef DM_DU0
/* Routines for dDslash/du0 */

void free_fn_links_dmdu0_u1(ferm_links_u1_t *fn){
  free_t_links_u1(&fn->fat);
  invalidate_all_ferm_links_u1(fn);
}
#endif

