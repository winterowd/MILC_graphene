/******* d_congrad5_fn_u1.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- this version for "fat plus Naik" quark
   actions.  

   Previously called d_congrad5_fn_fewsums.c and d_congrad5_fn_tmp.c
*/

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* 10/02/01 C. DeTar Consolidated with tmp version */
/* 4/18/03  D.T.     Consolidate global sums (fewer reductions) */
/* 5/6/07   C. DeTar True residual stopping criterion now. */

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "src", and the initial guess and answer
   in "dest".  "resid" is the residual vector, and "cg_p" and "ttt" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations before restarting.
   max_restarts = max number of restarts
   rsqmin = desired rsq, quit when we reach rsq <= rsqmin*source_norm.

   reinitialize after niters iterations and try once more.
*/

/* U(1) modification of conjugate gradient routine for SU(3) 4/13/2010*/

#include "generic_ks_includes_u1.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

/*#define CG_DEBUG*/

#define LOOPEND
#include "../include/loopend.h"

/* The Fermilab relative residue */

static Real 
my_relative_residue(complex *p, complex *q, int parity)
{
  double residue, num, den;
  int i;
  site *s;
  
  residue = 0;
  FORMYSITESANDPARITY(i,s,parity){
    num = (double) cabs_sq( &(p[i]) );
    den = (double) cabs_sq( &(q[i]) );
    residue += (den==0) ? 1.0 : (num/den);
  } END_LOOP

  g_doublesum(&residue);

  if(parity == EVENANDODD)
    return sqrt(residue/volume);
  else
    return sqrt(2*residue/volume);
}

static int
ks_congrad_parity_u1( complex *t_src, complex *t_dest, 
		   quark_invert_control *qic, Real mass,
		   ferm_links_t *fn){
  register int i;
  register site *s;
  int iteration;	/* counter for iterations */
  Real a,b;           	/* Sugar's a,b */
#ifdef FEWSUMS
  double actual_rsq = 999.;      /* rsq from actual summation of resid */
  double c_tr,c_tt,tempsum[4];	/* Re<resid|ttt>, <ttt|ttt> */
#endif
  double rsq = 0,relrsq = 0; /* resid**2, rel resid*2 */
  double oldrsq,pkp;	/*last resid*2,pkp = cg_p.K.cg_p */
  Real msq_x4;	/* 4*mass*mass */
  double source_norm;	/* squared magnitude of source vector */
  int otherparity = 0; /* the other parity */
  msg_tag * tags1[16], *tags2[16];	/* tags for gathers to parity and opposite */
  int special_started = 0; /* 1 if dslash_fn_field_special has been called */
  int nrestart;  /* Restart counter */
  complex *ttt, *cg_p, *resid;
  complex temp_mul, temp_conjg; //added for hold result of complex multiplication  and taking a complex conjugate 4/13/2010
  char myname[] = "ks_congrad_parity";

  /* Unpack structure */
  int niter        = qic->max;      /* maximum number of iters per restart */
  int max_restarts = qic->nrestart; /* maximum restarts */
  Real rsqmin      = qic->resid;    /* desired residual - 
			 normalized as sqrt(r*r)/sqrt(src_e*src_e) */
  Real relrsqmin   = qic->relresid; /* desired relative residual (FNAL)*/
  int parity       = qic->parity;   /* EVEN, ODD */

  int max_cg = max_restarts*niter; /* Maximum number of iterations */

  msq_x4 = 4.0*mass*mass;

  switch(parity){
  case(EVEN): otherparity=ODD; break;
  case(ODD):  otherparity=EVEN; break;
  }

  /* Allocate temporary variables */
  /* PAD may be used to avoid cache trashing */
#define PAD 0
  ttt = (complex *) malloc((sites_on_node+PAD)*sizeof(complex));
  cg_p = (complex *) malloc((sites_on_node+PAD)*sizeof(complex));
  resid = (complex *) malloc((sites_on_node+PAD)*sizeof(complex));

  if(ttt == NULL || cg_p == NULL || resid == NULL){
    printf("%s(%d): No room for temporaries\n",myname,this_node);
  }

  /* Source norm */
  source_norm = 0.0;
  FORMYSITESANDPARITY(i,s,parity){
    source_norm += (double)cabs_sq( &(t_src[i]) );
  } END_LOOP
  g_doublesum( &source_norm );
#ifdef CG_DEBUG
  node0_printf("congrad: source_norm = %e\n", (double)source_norm);
#endif

  /* Start CG iterations */
  
  nrestart = 0;
  iteration = 0;
  qic->size_r = 0;
  qic->size_relr = 0;

  while(1) {
    /* Check for completion */
    if( ( iteration % niter == 0 ) || 
	( ( rsqmin    <= 0 || rsqmin    > qic->size_r   ) &&
	  ( relrsqmin <= 0 || relrsqmin > qic->size_relr) ) ) 
      {
	
	/* (re)initialization process */
	
	/* Compute true residual and relative residual */

	/* ttt <-  (-1)*M_adjoint*M*dest
	   resid,cg_p <- src + ttt
	   rsq = |resid|^2
	   source_norm = |src|^2
	*/
	if(special_started==1) {	/* clean up gathers */
	    cleanup_gathers(tags1,tags2);
	    special_started=0;
	}
	rsq = 0.0;
	dslash_fn_field(t_dest, ttt, otherparity, fn);
	dslash_fn_field(ttt, ttt, parity, fn);
	/* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
	FORMYPARITYDOMAIN(i,s,parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( &ttt[i+FETCH_UP], 
			   &t_dest[i+FETCH_UP],
			   &t_src[i+FETCH_UP],
			   &resid[i+FETCH_UP]);
	  }
	  CMULREAL( t_dest[i], -msq_x4, temp_mul);
	  CADD( temp_mul , ttt[i], ttt[i]);
	  CADD( t_src[i], ttt[i], resid[i]);
	  /* remember ttt contains -M_adjoint*M*src */
	  cg_p[i] = resid[i];
	  rsq += (double)cabs_sq( &(resid[i]) );
	} END_LOOP
#ifdef FEWSUMS
	actual_rsq = rsq; /* not yet summed over nodes */
#endif
        g_doublesum( &rsq );

	if(relrsqmin > 0)
	  relrsq = my_relative_residue(resid, t_dest, parity);

	qic->final_rsq    = (Real)rsq/source_norm;
	qic->final_relrsq = (Real)relrsq;


	iteration++ ;  /* iteration counts number of multiplications
			  by M_adjoint*M */
	total_iters++;

#ifdef CG_DEBUG
	if(this_node==0)printf("CONGRAD: (re)start %d rsq = %.10e relrsq %.10e\n",
			       nrestart, qic->final_rsq, qic->final_relrsq);
#endif
	/* Quit when true residual and true relative residual are within
	   tolerance or when we exhaust iterations or restarts */
	
	if( iteration >= max_cg || 
	    nrestart  >= max_restarts ||
	    ( ( rsqmin    <= 0 || rsqmin    > qic->final_rsq   ) &&
	      ( relrsqmin <= 0 || relrsqmin > qic->final_relrsq) ) ) break;
	
	nrestart++;
      }

    /*
      oldrsq <- rsq
      ttt <- (-1)*M_adjoint*M*cg_p
      pkp <- (-1)*cg_p.M_adjoint*M.cg_p
      a <- -rsq/pkp
      dest <- dest + a*cg_p
      resid <- resid + a*ttt
      rsq <- |resid|^2
      b <- rsq/oldrsq
      cg_p <- resid + b*cg_p
    */

#ifdef FEWSUMS
    oldrsq = actual_rsq;	/* not yet summed over nodes */
#else
    oldrsq = rsq;
#endif
    /* sum of neighbors */
    
    if(special_started==0){
      dslash_fn_field_special( cg_p, ttt, otherparity, tags2, 1, fn );
      dslash_fn_field_special( ttt, ttt, parity, tags1, 1, fn);
      special_started=1;
    }
    else {
      dslash_fn_field_special( cg_p, ttt, otherparity, tags2, 0, fn );
      dslash_fn_field_special( ttt, ttt, parity, tags1, 0, fn);
    }
    
    /* finish computation of M_adjoint*M*p and p*M_adjoint*M*Kp */
    /* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
    /* pkp  <- cg_p.(ttt - msq*cg_p) */
    pkp = 0.0;
#ifdef FEWSUMS
    c_tr=0.0; c_tt=0.0;
#endif
    FORMYPARITYDOMAIN(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_VV( &ttt[i+FETCH_UP], &cg_p[i+FETCH_UP] );
      }
      CMULREAL( cg_p[i], -msq_x4, temp_mul);
      CADD( ttt[i], temp_mul, ttt[i]);
      CONJG( cg_p[i], temp_conjg);
      CMUL( temp_conjg, ttt[i] , temp_mul);
      pkp += (double)temp_mul.real;
#ifdef FEWSUMS
      CONJG( ttt[i], temp_conjg);
      CMUL( temp_conjg, resid[i], temp_mul);
      c_tr += (double)temp_mul.real;
      CMUL( temp_conjg, ttt[i], temp_mul );
      c_tr += (double)temp_mul.real;
#endif
    } END_LOOP
#ifdef FEWSUMS
    /* finally sum oldrsq over nodes, also other sums */
    tempsum[0] = pkp; tempsum[1] = c_tr; 
    tempsum[2] = c_tt; tempsum[3] = oldrsq;
    g_vecdoublesum( tempsum, 4 );
    pkp = tempsum[0]; c_tr = tempsum[1]; 
    c_tt = tempsum[2]; oldrsq = tempsum[3];
#else
    g_doublesum( &pkp );
#endif
    iteration++;
    total_iters++;
    
    a = (Real) (-rsq/pkp);
    
    /* dest <- dest + a*cg_p */
    /* resid <- resid + a*ttt */
#ifdef FEWSUMS
    actual_rsq=0.0;
#else
    rsq=0.0;
#endif
    FORMYPARITYDOMAIN(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_VVVV( &t_dest[i+FETCH_UP], 
		       &cg_p[i+FETCH_UP], 
		       &resid[i+FETCH_UP], 
		       &ttt[i+FETCH_UP] );
      }
      CMULREAL( cg_p[i], a, temp_mul);
      CADD( t_dest[i], temp_mul, t_dest[i] );
      CMULREAL( ttt[i], a, temp_mul);
      CADD( resid[i], temp_mul, resid[i] );
#ifdef FEWSUMS
      actual_rsq += (double)cabs_sq( &resid[i] );
#else
      rsq += (double)cabs_sq( &resid[i] );
#endif
    } END_LOOP
#ifdef FEWSUMS
    /**printf("XXX:  node %d\t%e\t%e\t%e\n",this_node,oldrsq,c_tr,c_tt);**/
    rsq = oldrsq + 2.0*a*c_tr + a*a*c_tt; /*TEST - should equal actual_rsq */
    /**c_tt = actual_rsq;**/ /* TEMP for test */
    /**g_doublesum(&c_tt);**/ /* TEMP true value for rsq */
    /**node0_printf("RSQTEST: %e\t%e\t%e\n",rsq,c_tt,rsq-c_tt);**/
#else
    g_doublesum(&rsq);
#endif	

    if(relrsqmin > 0)
      relrsq = my_relative_residue(resid, t_dest, parity);
    
    qic->size_r    = (Real)rsq/source_norm;
    qic->size_relr = (Real)relrsq;

#ifdef CG_DEBUG
    if(mynode()==0){printf("iter=%d, rsq/src= %e, relrsq= %e, pkp=%e\n",
			   iteration,(double)qic->size_r,
			   (double)qic->size_relr,
			   (double)pkp);fflush(stdout);}
#endif
    
    b = (Real)rsq/oldrsq;
    /* cg_p  <- resid + b*cg_p */
    FORMYSITESANDPARITY(i,s,parity){
      CMULREAL( cg_p[i], b, temp_mul );
      CADD( resid[i], temp_mul, cg_p[i] );
    } END_LOOP
  }

  if(nrestart == max_restarts || iteration == max_cg){
    node0_printf("ks_congrad: CG not converged after %d iterations and %d restarts, \n",
		 iteration, nrestart);
    node0_printf("rsq. = %e wanted %e relrsq = %e wanted %e\n",
		 qic->final_rsq,rsqmin,qic->final_relrsq,relrsqmin);
    fflush(stdout);
  }

  if(special_started==1) {
    cleanup_gathers(tags1,tags2);
    special_started = 0;
  }
  cleanup_dslash_temps();

  free(ttt); free(cg_p); free(resid);
  return iteration;
}

/* API for field arguments */

int ks_congrad_field_u1( complex *src, complex *dest, 
		      quark_invert_control *qic, Real mass,
		      ferm_links_t *fn)
{
  int iters = 0;
  double dtimec;
  double nflop = 1187;
  int parity = qic->parity;

  if(parity==EVENANDODD)nflop *=2;

  dtimec = -dclock(); 

  if(parity == EVEN || parity == EVENANDODD){
    qic->parity = EVEN;
    iters += ks_congrad_parity(src, dest, qic, mass, fn);
  }
  if(parity == ODD || parity == EVENANDODD){
    qic->parity = ODD;
    iters += ks_congrad_parity(src, dest, qic, mass, fn);
  }

  qic->parity = parity;

  dtimec += dclock();
#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (fn %s) masses = 1 iters = %d mflops = %e\n",
	   dtimec, prec_label[PRECISION-1], iters, 
	   (double)(nflop*volume*iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

  return iters;
}

//Site routines have been removed
