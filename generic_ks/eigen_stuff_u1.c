/****** eigen_stuff_u1.c  ******************/
/* Eigenvalue and Eigevector computation routines.
* K.O. 8/99 Started. 
* UH did some stuff to this I think
* EBG 2/15/02 changed phi -> phi1, mass -> mass1 so could be used 
*    w/ ks_imp_dyn2
* EBG 6/2004 fixed some memory leaks
* CW 9/30/2010 Converted all auxilary routines to U(1)
* MIMD version 7
*
*  These routines are for the computation of the Eigenvalues and Eigevectors
* of the Kogut-Susskind dslash^2. 
*/
//#define DEBUG
/* This Variable controls the tollerence of the JACOBI iteration */
#define JACOBI_TOL 1.110223e-16
#define MINITER 5

/*   If you use strict convergence then you stop when g = |M*v - l*v|     *
 *  becomes less than the eigenvalue tolerance. Other wise you stop       *
 *  when the estimated error of the eigenvalue gets smaller than          *
 *  the eigenvalue tolerance. The later option takes abvantage of         *
 *  the quadratic convergence of the algorithm. (See Kalkreuter's paper    * 
 *  for details). I prefer the strict convergenve which results some      *
 *  overkill.                                                             */
/**#define STRICT_CONVERGENCE	(taken out by UMH)**/


/* If you do not define the USE_DSLASH_SPECIAL then calls to the standard   *
 * dslash are made. If you do define USE_DSLASH_SPECIAL then DSLASH_SPECIAL *
 * is used                                                                  */
#define USE_DSLASH_SPECIAL


/* Include files */
#include "generic_ks_includes_u1.h"
#include "../include/jacobi.h"
#include "../include/dslash_ks_redefine.h"
#include <string.h>
#define LOOPEND
#include "../include/loopend.h"

static struct {
  int d[4];
  Real sign ;
} eps[24]={
  {{0,1,2,3},+1.0},
  {{0,3,1,2},+1.0},
  {{0,2,3,1},+1.0},
  {{0,3,2,1},-1.0},
  {{0,1,3,2},-1.0},
  {{0,2,1,3},-1.0},

  {{1,0,2,3},-1.0},
  {{1,3,0,2},-1.0},
  {{1,2,3,0},-1.0},
  {{1,3,2,0},+1.0},
  {{1,0,3,2},+1.0},
  {{1,2,0,3},+1.0},

  {{2,1,0,3},-1.0},
  {{2,3,1,0},-1.0},
  {{2,0,3,1},-1.0},
  {{2,3,0,1},+1.0},
  {{2,1,3,0},+1.0},
  {{2,0,1,3},+1.0},

  {{3,1,2,0},-1.0},
  {{3,0,1,2},-1.0},
  {{3,2,0,1},-1.0},
  {{3,0,2,1},+1.0},
  {{3,1,0,2},+1.0},
  {{3,2,1,0},+1.0}
};

void Matrix_Vec_mult_u1(complex *src, complex *res, int parity,
		     ferm_links_u1_t *fn );
void cleanup_Matrix() ;
void measure_chirality_u1(complex *src, double *chirality, int parity) ;
void print_densities_u1(complex *src, char *tag, int y,int z,int t,int parity);

void GramSchmidt_u1(complex **vector, int Num, int parity) ;/* Not used */
void copy_Vector_u1(complex *src, complex *res ) ; 
void norm2_u1(complex *vec, double *norm, int parity); 
void dot_product_u1(complex *vec1, complex *vec2, 
		 double_complex *dot, int parity) ;
void complex_vec_mult_sub_u1(double_complex *cc, complex *vec1, 
			  complex *vec2, int parity) ;
void complex_vec_mult_add_u1(double_complex *cc, complex *vec1, 
			  complex *vec2, int parity) ;
void double_vec_mult_u1(double *a, complex *vec1, 
		     complex *vec2, int parity) ;
void double_vec_mult_sub_u1(double *rr, complex *vec1,  
			 complex *vec2, int parity) ;
void double_vec_mult_add_u1(double *rr, complex *vec1,  
			 complex *vec2, int parity) ;
void dax_p_by_u1(double *a, complex *vec1, double *b, complex *vec2, 
	      int parity) ; 
void vec_plus_double_vec_mult_u1(complex *vec1, double *a, complex *vec2, 
			      int parity) ;
void normalize_u1(complex *vec,int parity) ;
void project_out_u1(complex *vec, complex **vector, int Num, int parity);
void constructArray_u1(complex **eigVec, complex **MeigVec, Matrix *A,
		    double *err, int *converged, int parity,
		    ferm_links_u1_t *fn);
void RotateBasis_u1(complex **eigVec, Matrix *V, int parity) ;

/*Can move these U(1) "flavor_ops" routines to new file if 
  need the rest later */
void mult_spin_pseudoscalar_u1(field_offset src, field_offset dest ) ;
void eta_shift_u1(int n, int *d, field_offset src, field_offset dest );
void sym_shift_u1(int dir, field_offset src,field_offset dest);

/************************************************************************/


/* Message tags to be used for the Matrix Vector multiplication */
static msg_tag *tags1[16],*tags2[16];
/* Temporary su3_vectors used for the squaring */ 
static complex *temp = NULL ;
/* flag indicating wether to start dslash. */
static int dslash_start = 1 ; /* 1 means start dslash */
/* This routine is the routine that applies the Matrix, whose eigenvalues    *
 * we want to compute, to a vector. For this specific application it is the  *
 * -D_slash^2 of the KS fermions. We only compute on the "parity" sites.     *
 * Where parity can be EVEN, ODD, or ENENANDODD                              */


#ifndef USE_DSLASH_SPECIAL
/* The Matrix_Vec_mult_u1 and cleanup_Matrix() WITHOUT using dslash_special updated 9/29/2010 CW */
/************************************************************************/

void Matrix_Vec_mult_u1(complex *src, complex *res, int parity,
		     ferm_links_u1_t *fn )
{  
  register site *s;
  register  int i;
  int otherparity = EVENANDODD;

  if(temp == NULL ){
    temp = (complex *)malloc(sites_on_node*sizeof(complex));
  }
  
  switch(parity){
  case EVEN:
    otherparity = ODD ;
    break ;
  case ODD:
    otherparity = EVEN ;
    break ;
  case EVENANDODD:
    otherparity = EVENANDODD ;
    break ;
  default:
    node0_printf("ERROR: wrong parity in eigen_stuff::Matrix_Vec_mult\n") ;
    terminate(1) ;
  }

  dslash_fn_field_u1(src , temp, otherparity, fn) ; 
  dslash_fn_field_u1(temp, res , parity     , fn) ;
  FORMYSITESANDPARITY(i,s,parity){ 
    CMULREAL( res[i], -1.0, res[i] );
   } END_LOOP
}
/*****************************************************************************/
/* Deallocates the tags and the temporaries the Matrix_Vec_mult needs */
void cleanup_Matrix(){
  if(temp != NULL ){
    free(temp) ;
    temp = NULL ;
  }
}
#else
/*****************************************************************************/
/* The Matrix_Vec_mult_1 and cleanup_Matrix() using dslash_special  Updated 9/29/2010 CW */
void Matrix_Vec_mult_u1(complex *src, complex *res, int parity,
		     ferm_links_u1_t *fn ){
  
  register site *s;
  register  int i;
  int otherparity = EVENANDODD;
  /* store last source so that we know when to reinitialize the message tags */
  static complex *last_src=NULL ;

  if(dslash_start){
    temp = (complex *)malloc(sites_on_node*sizeof(complex));
  }

  /*reinitialize the tags if we have a new source */
  if(last_src != src){
    if(!dslash_start) cleanup_gathers(tags1,tags2); 
    dslash_start = 1 ;
    last_src = src ;
  }
  switch(parity){
  case EVEN:
    otherparity = ODD ;
    break ;
  case ODD:
    otherparity = EVEN ;
    break ;
  case EVENANDODD:
    otherparity = EVENANDODD ;
    break ;
  default:
    node0_printf("ERROR: wrong parity in eigen_stuff::Matrix_Vec_mult\n") ;
    terminate(1) ;
  }
  
  dslash_fn_field_special_u1(src , temp, otherparity, tags1, dslash_start, fn) ;
  dslash_fn_field_special_u1(temp, res , parity     , tags2, dslash_start, fn) ;
  
  FORMYSITESANDPARITY(i,s,parity){ 
    CMULREAL( res[i], -1.0, res[i] );
  } END_LOOP
  dslash_start = 0 ;
}
/* Deallocates the tags and the temporaries the Matrix_Vec_mult needs */

/************************************************************************/
void cleanup_Matrix(){
  if(!dslash_start) {
    cleanup_gathers(tags1,tags2); 
    cleanup_dslash_temps() ;
    free(temp) ;
  }
  dslash_start = 1 ;
#ifdef DEBUG
  node0_printf("cleanup_Matrix(): done!\n") ; fflush(stdout) ;
#endif
}
#endif /* USE_DSLASH_SPECIAL */

/************************************************************************/

/* Modified Gram-Schmidt orthonormalization pg. 219 Golub & Van Loan     *
 * Num is the number of vectors. vector are the vectors. They get        *
 * overwritten by the orthonormal vectors.                               *
 * parity is the parity on which we work on (EVEN,ODD,ENENANDODD).       */
void GramSchmidt_u1(complex **vector, int Num, int parity)/* NEVER USED */{
  register int i,j ;
  double norm ;
  double_complex cc ;

  for(i=0;i<Num;i++){
    norm2_u1(vector[i], &norm, parity);
    norm = 1.0/norm ;
    double_vec_mult_u1(&norm, vector[i], vector[i], parity);
    for(j=i+1;j<Num;j++){
      dot_product_u1(vector[i], vector[j], &cc, parity) ;
      complex_vec_mult_sub_u1(&cc, vector[i], vector[j], parity);
    }
  }
}
/************************************************************************/

/*  Projects out the *vectors from the  vec. Num is the Number of vectors  *
 * and parity is the parity on which we work on.                           *
 * The vectors are assumed to be orthonormal.                              */
void project_out_u1(complex *vec, complex **vector, int Num, int parity){
  register int i ;
  double_complex cc ;
  for(i=Num-1;i>-1;i--){
    dot_product_u1(vector[i], vec, &cc, parity) ;
    complex_vec_mult_sub_u1(&cc, vector[i], vec, parity);
  }
}

/************************************************************************/
/* Copies scr to res */
void copy_Vector_u1(complex *src, complex *res){
  memcpy((void *)res, (void *)src, sites_on_node*sizeof(complex)) ;
}

/************************************************************************/
/* Returns the 2-norm of a fermion vector */
void norm2_u1(complex *vec, double *norm, int parity){
  register double n ;
  register site *s;
  register  int i;
  
  n=0 ; 
  FORMYSITESANDPARITY(i,s,parity){
    n += (double)cabs_sq(&(vec[i]));
  } END_LOOP
  *norm = n ;
  g_doublesum(norm);
  *norm = sqrt(*norm) ;
}
 
/*****************************************************************************/
/* Returns the dot product of two fermion vectors */
void dot_product_u1(complex *vec1, complex *vec2, 
		   double_complex *dot, int parity) {
  register double re,im ;
  register site *s;
  register  int i;
  complex cc ;
  
  re=im=0.0;
  FORMYSITESANDPARITY(i,s,parity){
    CMULJ_( vec1[i], vec2[i], cc ); 
    re += cc.real ;
    im += cc.imag ;
  } END_LOOP
  dot->real = re ; 
  dot->imag = im ;
  g_dcomplexsum(dot);
}

/*****************************************************************************/
/* Returns vec2 = vec2 - cc*vec1   cc is a double complex   */
void complex_vec_mult_sub_u1(double_complex *cc, complex *vec1, 
			  complex *vec2, int parity){

  register site *s;
  register  int i;
  complex sc, temp_mul;
  
  sc.real= (Real)(cc->real) ; 
  sc.imag= (Real)(cc->imag) ;

  FORMYSITESANDPARITY(i,s,parity){
    CMUL( vec1[i], sc, temp_mul );
    CSUB( vec2[i], temp_mul, vec2[i] );
  } END_LOOP
}

/*****************************************************************************/
/* Returns vec2 = vec2 + cc*vec1   cc is a double complex   */
void complex_vec_mult_add_u1(double_complex *cc, complex *vec1, 
			  complex *vec2, int parity){
  register site *s;
  register  int i;
  complex sc, temp_mul;
  
  sc.real= (Real)(cc->real) ; sc.imag= (Real)(cc->imag) ;
  FORMYSITESANDPARITY(i,s,parity){
    CMUL( vec1[i], sc, temp_mul );
    CADD( vec2[i], temp_mul, vec2[i] );
  } END_LOOP
}

/*****************************************************************************/
/* Returns vec2 = vec2 - rr*vec1   rr is a double    */
void double_vec_mult_sub_u1(double *rr, complex *vec1,  
			 complex *vec2, int parity){
  register site *s;
  register  int i;
  complex temp_mul;
  
  FORMYSITESANDPARITY(i,s,parity){
    CMULREAL( vec1[i], *rr, temp_mul );
    CSUB( vec2[i], temp_mul, vec2[i] );
  } END_LOOP
}

/*****************************************************************************/
/* Returns vec2 = vec2 + rr*vec1   rr is a double    */
void double_vec_mult_add_u1(double *rr, complex *vec1,  
			 complex *vec2, int parity){
  register site *s;
  register  int i;
  complex temp_mul;

  FORMYSITESANDPARITY(i,s,parity){
    CMULREAL( vec1[i], *rr, temp_mul );
    CADD( vec2[i], temp_mul, vec2[i] );
  } END_LOOP
}

/*****************************************************************************/
/* Returns vec2 = a*vec1   a is a double vec2 can be vec1*/
void double_vec_mult_u1(double *a, complex *vec1, 
		     complex *vec2, int parity){
  
  register site *s;
  register  int i;
  
  FORMYSITESANDPARITY(i,s,parity){ 
    CMULREAL( vec1[i], *a, vec2[i] );
  } END_LOOP
}

/*****************************************************************************/
/* Returns vec2 = vec1 + a*vec2 */
void vec_plus_double_vec_mult_u1(complex *vec1, double *a, complex *vec2, 
			      int parity){
  register site *s;
  register  int i;

  FORMYSITESANDPARITY(i,s,parity){ 
    CMULREAL( vec2[i], *a, vec2[i] );
    CADD( vec2[i], vec1[i], vec2[i] );
  } END_LOOP
}

/*****************************************************************************/
/* Returns vec1 = a*vec1 + b*vec2   a,b are double    */
void dax_p_by_u1(double *a, complex *vec1, double *b, complex *vec2, 
	      int parity) {

  register site *s;
  register  int i;
  complex temp_mul;
  
  FORMYSITESANDPARITY(i,s,parity){
    CMULREAL( vec1[i], *a, vec1[i] );
    CMULREAL( vec2[i], *b, temp_mul );
    CADD( vec1[i], temp_mul, vec1[i] );
  } END_LOOP
}

/*****************************************************************************/
/* normalizes the vecror vec. Work only on parity. */
void normalize_u1(complex *vec,int parity){

  double norm ;
  norm2_u1(vec,&norm,parity) ;
  norm = 1.0/norm ;
  double_vec_mult_u1(&norm,vec,vec,parity) ;
}

/*****************************************************************************/
/* Run a CG iteration to minimize the Rayleigh-Ritz quotient over the
   span "eigVec".  The initial vector is "vec".  The CG iteration stops
   when the tolerance is reached or the norm of the residual drops
   faster than RelTol, or when the MaxIter and Restart are exhausted. */

int Rayleigh_min_u1(complex *vec, complex **eigVec, Real Tolerance, 
		 Real RelTol, int Nvecs, int MaxIter, int Restart, 
		 int parity, ferm_links_u1_t *fn){

  int iter ;
  double beta, cos_theta, sin_theta ;
  double quot, old_quot, P_norm, theta, real_vecMp, pMp ;
#ifdef DEBUG
  double vec_norm ;
#endif
  double g_norm, old_g_norm, start_g_norm ;
  double_complex cc ;
  complex *Mvec, *grad, *P, *MP ;
  
  Mvec     = (complex *)malloc(sites_on_node*sizeof(complex)) ;
  grad     = (complex *)malloc(sites_on_node*sizeof(complex)) ;
  P        = (complex *)malloc(sites_on_node*sizeof(complex)) ;
  MP       = (complex *)malloc(sites_on_node*sizeof(complex)) ;


  old_quot = 1.0e+16 ;
  project_out_u1(vec, eigVec, Nvecs, parity);
  normalize_u1(vec,parity) ; 
  Matrix_Vec_mult_u1(vec,Mvec,parity, fn) ;
  project_out_u1(Mvec, eigVec, Nvecs, parity);
  
  /* Compute the quotient quot=vev*M*vec */
  dot_product_u1(vec, Mvec, &cc, parity) ;
  /* quot is real since M is hermitian. quot = vec*M*vec */
  quot = cc.real ;
#ifdef DEBUG
  node0_printf("Rayleigh_min: Start -- quot=%g,%g\n", quot,cc.imag) ;
#endif  
  /* Compute the grad=M*vec - quot*vec */
  copy_Vector_u1(Mvec,grad) ;
  double_vec_mult_sub_u1(&quot, vec, grad, parity) ;
  /* set P (the search direction) equal to grad */
  copy_Vector_u1(grad,P) ;
  /* compute the norms of P and grad */
  norm2_u1(P   , &P_norm, parity) ;
  g_norm = P_norm ;
  start_g_norm = g_norm ;
#ifdef DEBUG
  node0_printf("Rayleigh_min: Start -- g_norm=%g\n", g_norm) ;
  node0_printf("Tolerance: %e, MaxIter: %d, MINITER: %d\n", Tolerance, MaxIter,
	       MINITER);
#endif  

  iter = 0 ;
  while( (g_norm>Tolerance)&&
	 ( ((iter<MaxIter)&&(g_norm/start_g_norm>RelTol)) || (iter<MINITER) ) 
	 ){
    iter++ ;
    Matrix_Vec_mult_u1(P,MP,parity, fn) ;
    dot_product_u1(vec, MP, &cc, parity) ;
    real_vecMp = cc.real ;
    dot_product_u1(P, MP, &cc, parity) ;
    pMp = cc.real ; /*pMp is real */
    theta = 0.5*atan(2.0*real_vecMp/(quot*P_norm - pMp/P_norm)) ;
    sin_theta = sin(theta) ;
    cos_theta = cos(theta) ;
    if(sin_theta*cos_theta*real_vecMp>0){
      theta = theta - 0.5*PI ; /* chose the minimum not the maximum */
      sin_theta = sin(theta) ; /* the sin,cos calls can be avoided */
      cos_theta = cos(theta) ;
    }
    sin_theta = sin_theta/P_norm ;
    /* vec = cos(theta)*vec +sin(theta)*P/p_norm */
    dax_p_by_u1(&cos_theta, vec,&sin_theta, P,parity) ;
    /* Mvec = cos(theta)*Mvec +sin(theta)*MP/p_norm */
    dax_p_by_u1(&cos_theta,Mvec,&sin_theta,MP,parity) ;
    /* renormalize vec ... */
    if(iter%Restart == 0 ) {
#ifdef DEBUG
      node0_printf("Renormalizing...");
      norm2_u1(vec,&vec_norm,parity) ;
      node0_printf("  norm: %g\n",1.0/vec_norm);
#endif
      /* Project vec on the orthogonal complement of eigVec */
      project_out_u1(vec, eigVec, Nvecs, parity);
      normalize_u1(vec,parity) ;
      Matrix_Vec_mult_u1(vec,Mvec,parity, fn);
      /* Recompute the quotient */
      dot_product_u1(vec, Mvec, &cc, parity) ;
      /* quot is real since M is hermitian. quot = vec*M*vec */
      quot = cc.real ;
      /* Recompute the grad */
      copy_Vector_u1(Mvec,grad) ;
      double_vec_mult_sub_u1(&quot, vec, grad, parity) ;
      norm2_u1(grad,&g_norm,parity) ; /* recompute the g_norm */
      /* Project P on the orthogonal complement of eigVec */
      project_out_u1(P, eigVec, Nvecs, parity);
      /* make P orthogonal to vec */
      dot_product_u1(vec, P, &cc, parity) ;
      complex_vec_mult_sub_u1(&cc, vec, P, parity);
      /* make P orthogonal to grad */
      dot_product_u1(grad, P, &cc, parity) ;
      complex_vec_mult_sub_u1(&cc, grad, P, parity);
      norm2_u1(P,&P_norm,parity) ; /* recompute the P_norm */
    }
    dot_product_u1(vec, Mvec, &cc, parity) ;
    /* quot is real since M is hermitian. quot = vec*M*vec */
    quot = cc.real ;
#ifdef DEBUG
    node0_printf("Rayleigh_min: %i, quot=%8g g=%8g b=%6g P:%6g\n",
		 iter,cc.real,g_norm,beta,P_norm) ;
#endif      
    old_g_norm = g_norm ;
       
    copy_Vector_u1(Mvec,grad) ;
    double_vec_mult_sub_u1(&quot, vec, grad, parity) ;

    norm2_u1(grad,&g_norm,parity) ;
    beta = cos_theta*g_norm*g_norm/(old_g_norm*old_g_norm) ;
    /* Cut off beta */
    if( beta>2.0 ) {
      beta = 2.0 ;
    }
    dot_product_u1(vec, P, &cc, parity) ;
    cc.real *= beta ; cc.imag *= beta ;
    /*      P_norm = 1.0/P_norm ;
	    double_vec_mult(&P_norm,  P,  P, parity) ; */
    vec_plus_double_vec_mult_u1(grad, &beta, P, parity) ;/* P = grad + beta*P */
    complex_vec_mult_sub_u1(&cc, vec, P, parity) ; /* P = P - cc*vec */
    norm2_u1(P   , &P_norm, parity) ;

    if(fabs(old_quot -quot)< Tolerance/100.0){
      /*break*/
      g_norm = Tolerance/10.0; 
    }
#ifdef DEBUG
    node0_printf("At end of while loop:\n");
    node0_printf("iter: %d, g_norm: %e, quot: %e, old_quot: %e\n", iter, g_norm, 
		 quot, old_quot);
#endif
    old_quot = quot ;

  }

  project_out_u1(vec, eigVec, Nvecs, parity);
  normalize_u1(vec,parity) ;
  cleanup_Matrix() ;
  free(MP) ;
  free(P) ;
  free(grad) ;
  free(Mvec) ;
  iter++ ;
  return iter ;
}

/*****************************************************************************/
/* Returns the projected matrix A and the error of each eigenvector */
void constructArray_u1(complex **eigVec, complex **MeigVec, Matrix *A,
		    double *err, int *converged, int parity,
		    ferm_links_u1_t *fn){
  int i,j,Nvecs ;
  complex *grad ;
  double_complex cc,Aij,Aji ;

  Nvecs = A->N ;
  grad = (complex *)malloc(sites_on_node*sizeof(complex));
  for(i=0;i<Nvecs;i++){
    if(converged[i]==0) Matrix_Vec_mult_u1(eigVec[i],MeigVec[i],parity, fn) ;
    dot_product_u1(MeigVec[i], eigVec[i], &cc, parity) ;
    A->M[i][i].real = cc.real ; A->M[i][i].imag = 0.0 ;
    copy_Vector_u1(MeigVec[i], grad) ;
    double_vec_mult_sub_u1(&cc.real, eigVec[i], grad, parity) ;
    norm2_u1(grad, &err[i], parity) ;
    for(j=i+1;j<Nvecs;j++){
      dot_product_u1(MeigVec[i], eigVec[j], &cc, parity) ;
      Aij=cc ;
      CONJG(cc, Aji) ;
      dot_product_u1(eigVec[j], MeigVec[i], &cc, parity) ;
      CSUM(Aji,cc) ;
      CONJG(cc, cc) ;
      CSUM(Aij,cc) ;
      CMULREAL(Aij,0.5,A->M[i][j]) ;
      CMULREAL(Aji,0.5,A->M[j][i]) ;
    }
  }
  free(grad) ;
}


/*****************************************************************************/
void RotateBasis_u1(complex **eigVec, Matrix *V, int parity){

  complex **Tmp ;
  register site *s;
  register  int i,N ;
  int j,k ;
  
  N = V->N ;
  /* Allocate the temporary vectors needed */
  Tmp = (complex **)malloc(N*sizeof(complex *));
  for(j=0;j<N;j++){
    Tmp[j]=(complex *)malloc(sites_on_node*sizeof(complex));
  }

  for(j=0;j<N;j++){
    FORMYSITESANDPARITY(i,s,parity){
      //clearvec(&Tmp[j][i]) ;
      Tmp[j][i].real = Tmp[j][i].imag = 0.0;
    } END_LOOP
    for(k=0;k<N;k++){
      complex_vec_mult_add_u1(&V->M[k][j],eigVec[k],Tmp[j],parity) ;
    }
  }
  
  /* Copy rotated basis to the eigVec and free temporaries */
  for(j=0;j<N;j++) {
    /* Copy only wanted parity (UMH) */
    /** copy_Vector(Tmp[j], eigVec[j]) ; **/
    FORMYSITESANDPARITY(i,s,parity){
      eigVec[j][i].real = Tmp[j][i].real;
      eigVec[j][i].imag = Tmp[j][i].imag;
    } END_LOOP
    free(Tmp[j]) ;
  }
  free(Tmp) ;
}

/*****************************************************************************/
int Kalkreuter_u1(complex **eigVec, double *eigVal, Real Tolerance, 
	       Real RelTol, int Nvecs, int MaxIter, 
	       int Restart, int Kiters, int parity,
	       ferm_links_u1_t *fn ){

  int total_iters=0 ;
  int j;
  Matrix Array,V ;
  register site *s ;
  register  int i ;
  complex *vec ;
  complex **MeigVec ;
  double max_error = 1.0e+10 ;
  double *grad, *err ;
  int iter = 0 ;
  int *converged ;
  Real ToleranceG ;
#ifdef EIGTIME
  double dtimec;
#endif

  ToleranceG = 10.0*Tolerance ;

  /** Allocate the array **/
  Array = AllocateMatrix(Nvecs) ;  
  /** Allocate the Eigenvector matrix **/
  V = AllocateMatrix(Nvecs) ;  

  vec = (complex *)malloc(sites_on_node*sizeof(complex));
  grad = (double *) malloc(Nvecs*sizeof(double)) ;
  err = (double *) malloc(Nvecs*sizeof(double)) ;
  converged = (int *)malloc(Nvecs*sizeof(int));
  MeigVec = (complex **)malloc(Nvecs*sizeof(complex*));
  for(i=0;i<Nvecs;i++){
    MeigVec[i] = (complex *)malloc(sites_on_node*sizeof(complex));
  }

  /* Initiallize all the eigenvectors to a random vector */
  for(j=0;j<Nvecs;j++) {
    grad[j] = 1.0e+10 ;
    //grsource_plain(F_OFFSET(g_rand), parity);  
    grsource_plain_u1(F_OFFSET(g_rand), parity);
    FORMYSITESANDPARITY(i,s,parity){
      eigVec[j][i].real = s->g_rand.real;
      eigVec[j][i].imag = s->g_rand.imag;
    } END_LOOP
    eigVal[j] = 1.0e+16 ;
  }

#ifdef EIGTIME
  dtimec = -dclock();
#endif

  while((max_error>Tolerance)&&(iter<Kiters)) {
    iter++ ;
    /* Run through all current eigenvectors, minimizing the
       Rayleigh-Ritz quotient in the space spanned by the current
       eigenvectors.  Replace the eigenvectors with each
       improvement. Stop improving the eigenvector when it has
       converged. */
    for(j=0;j<Nvecs;j++){
      if(grad[j]>(ToleranceG)) {
	converged[j] = 0 ;
	copy_Vector_u1(eigVec[j],vec) ;
	total_iters += Rayleigh_min_u1(vec, eigVec, ToleranceG, RelTol,
				    j, MaxIter , Restart, parity, fn) ;
	/* Copy only wanted parity (UMH) */
	/** copy_Vector(vec,eigVec[j]) ; **/
	FORMYSITESANDPARITY(i,s,parity){
	  eigVec[j][i].real = vec[i].real;
	  eigVec[j][i].imag = vec[i].imag;
	} END_LOOP
      }else{ 
	converged[j] = 1;
      }
    }

    /* Diagonalize the operator on the subspace of improved
       eigenvectors and rotate the eigenvectors to the new basis. */

    /* if you didn't act on eigVec[i] last time, converged[i]=1,
       and  MeigVec hasn't changed, so don't compute it */
    constructArray_u1(eigVec, MeigVec, &Array, grad, converged, parity, fn) ;

#ifdef DEBUG
    node0_printf("Eigenvalues before diagonalization\n");
    for(i=0;i<Nvecs;i++){
      node0_printf("quot(%i) = %g |grad|=%g\n",i,Array.M[i][i].real,grad[i]);
    }
#endif

    Jacobi(&Array, &V, JACOBI_TOL) ;
    sort_eigenvectors(&Array,&V) ;
    RotateBasis_u1(eigVec,&V,parity) ;
    RotateBasis_u1(MeigVec,&V,parity) ;

    /* find the maximum error */
    max_error = 0.0 ;
    for(i=0;i<Nvecs;i++){
      /* First recompute gradient--vec=MeigVec - eigVal*eigVec */
      /* recall we are recycling MeigVec's */
      Array.M[i][i].imag=0.0;
      copy_Vector_u1(MeigVec[i], vec) ;
      double_vec_mult_sub_u1(&Array.M[i][i].real, eigVec[i], vec, parity) ;
      norm2_u1(vec, &grad[i], parity) ;
      err[i] = eigVal[i] ;
      eigVal[i] = Array.M[i][i].real ;
      err[i] = fabs(err[i] - eigVal[i])/(1.0 - RelTol*RelTol) ;
#ifndef STRICT_CONVERGENCE
      if(err[i]>max_error){
	max_error = err[i] ;
      }
#else	
      if(grad[i]>max_error){
	max_error = grad[i] ;
      }
#endif
    }

    /* this makes a LOT of output! Im  commenting it out for a bit. -EBG */

    /*
    node0_printf("\nEigenvalues after diagonalization at iteration %i\n",
		 iter);
    for(i=0;i<Nvecs;i++)
      node0_printf("quot(%i) = %g +/- %8e |grad|=%g\n",
		   i,eigVal[i],err[i],grad[i]);

    */
  }
  
#ifdef EIGTIME
  dtimec += dclock();
  node0_printf("KAULKRITER: time = %e iters = %d iters/vec = %e\n",
	   dtimec,total_iters, (double)(total_iters)/Nvecs);
#endif

  node0_printf("Eigenvalues after diagonalization at iteration %i\n",iter);

  node0_printf("BEGIN RESULTS\n");
  for(i=0;i<Nvecs;i++){
    node0_printf("Eigenvalue(%i) = %g +/- %8e\t cvg? %d  \n",
		 i,eigVal[i],err[i],converged[i]);
  }

  /** Deallocate the arrays **/
  deAllocate(&V) ;
  deAllocate(&Array) ;
  free(err) ;
  free(grad) ;
  free(vec) ;
  for(i=0;i<Nvecs;i++){
    free(MeigVec[i]);
  }
  free(MeigVec);
  cleanup_Matrix();

  return total_iters ;
}



/*****************************************************************************/
/* measures the chiraliry of a normalized fermion state */
void measure_chirality_u1(complex *src, double *chirality, int parity){
  register int i;
  register site *s;
  register double cc ;
  complex tmp ;

  FORMYSITESANDPARITY(i,s,parity){
    s->tempvec[3].real = src[i].real;
    s->tempvec[3].imag = src[i].imag;
  } END_LOOP

  mult_spin_pseudoscalar_u1(F_OFFSET(tempvec[3]),F_OFFSET(ttt)) ;

  cc = 0.0 ; 
  FORMYSITESANDPARITY(i,s,parity){ 
    CMULJ_( s->tempvec[3], s->ttt, tmp );
    //tmp = su3_dot( &(s->tempvec[3]), &(s->ttt) ) ;
    cc +=  tmp.real ; /* chirality is real since Gamma_5 is hermitian */
  } END_LOOP
  *chirality = cc ;
  g_doublesum(chirality);
}


/*****************************************************************************/
/* prints the density and chiral density of a normalized fermion state */
void print_densities_u1(complex *src, char *tag, int y,int z,int t, 
		     int parity){

  register int i;
  register site *s;
  complex tmp1,tmp2 ;

  FORMYSITESANDPARITY(i,s,parity){
    s->tempvec[3].real = src[i].real;
    s->tempvec[3].imag = src[i].imag;
  } END_LOOP
      
  mult_spin_pseudoscalar_u1(F_OFFSET(tempvec[3]),F_OFFSET(ttt)) ;
  
  FORMYSITESANDPARITY(i,s,parity){ 
    if((s->y==y)&&(s->z==z)&&(s->t==t)){
      CMULJ_( s->tempvec[3], s->ttt, tmp1 );
      CMULJ_( s->tempvec[3], s->tempvec[3], tmp2 );
      node0_printf("%s: %i %e %e %e\n",tag,
		   s->x,tmp2.real,tmp1.real,tmp1.imag);
    }
  } END_LOOP
  
}/* print_densities_u1 */

/* Apply the symmetric shift opperator in direction "dir" *
 * This is the explicit version                           *
 * Covariant shifts are used                              */
void sym_shift_u1(int dir, field_offset src,field_offset dest)
{
  register int i ;
  register site *s ;
  msg_tag *tag[2];
  complex *tvec;
  
  tvec = (complex *)malloc( sites_on_node*sizeof(complex) );

  tag[0] = start_gather_site( src, sizeof(complex), dir, EVENANDODD ,gen_pt[0] );
  FORALLMYSITES(i,s)
    {
      CMULJ_( s->link[dir], (*(complex *)F_PT(s,src)), tvec[i]);
      /*mult_adj_su3_mat_vec( &(s->link[dir]), (complex *)F_PT(s,src), 
	&(tvec[i]) ) ;*/
    }
  tag[1] = start_gather_field(tvec, sizeof(complex), OPP_DIR(dir), 
				  EVENANDODD ,gen_pt[1] );
  wait_gather(tag[0]);
  FORALLMYSITES(i,s)
    {
      CMUL( s->link[dir], (*(complex *)gen_pt[0][i]), 
	    (*(complex *)F_PT(s,dest) ) );
      /*mult_su3_mat_vec( &(s->link[dir]), (complex *)gen_pt[0][i], 
	(complex *)F_PT(s,dest) ) ;  */  
    }
  wait_gather(tag[1]);
  FORALLMYSITES(i,s)
    {
      CADD( (*(complex *)F_PT(s,dest)),(*(complex *)gen_pt[1][i]), 
	    (*(complex *)F_PT(s,dest)) ); 
      /*add_su3_vector( (su3_vector *)F_PT(s,dest), (su3_vector *)gen_pt[1][i], 
	(su3_vector *)F_PT(s,dest) ) ;    */
    }
  /* Now devide by 2 eq. (4.2b) of Golderman's Meson paper*/
 FORALLMYSITES(i,s)
   {
     CMULREAL( (*(complex *)F_PT(s,dest)), .5, (*(complex *)F_PT(s,dest)) );
     /*scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), .5,
       (su3_vector *)F_PT(s,dest) );*/
   }
  for(i=0;i<2;i++) cleanup_gather(tag[i]) ;
  free(tvec);
} /* sym_shift_u1 */

/* It applies the symmetric shift with directions                       *
 * stored in the array d. Each shift is multiplied by \eta_k            *
 * In fact since \eta_k are already absorbed into the U matrices do     *
 * no action is needed. Just do the symmetric shift.                    *
 * n is the number of shifts                                            */
void eta_shift_u1(int n, int *d, field_offset src, field_offset dest )
{
  int c ;
  field_offset temp0,temp1,tmp ;
  
  if(n==1)
    {
      sym_shift_u1(d[0], src, dest);
    }
  else
    {
      temp0 = F_OFFSET(tempvec[0]) ;
      temp1 = F_OFFSET(tempvec[1]) ;
      if(dest==temp0||dest==temp1)
	{ 
	  node0_printf("eta_shift(): ERROR! dest should not be ");
	  node0_printf("tempvec[0] or tempvec[1]\n" );
	  terminate(102) ;
	}
      if(src==temp0||src==temp1)
	{ 
	  node0_printf("eta_shift(): ERROR! src should not be ");
	  node0_printf("tempvec[0] or tempvec[1]\n" );
	  terminate(101) ;
	}
      sym_shift_u1(d[0], src, temp0 );
      for(c=1;c<n-1;c++)
	{  
	  sym_shift_u1(d[c], temp0, temp1);
	  /* switch the pointers */
	  tmp = temp0   ;
	  temp0 = temp1 ;
	  temp1 = tmp   ;
	}
      /* do the last shift */
      sym_shift_u1(d[n-1], temp0, dest );
    }
} /* eta_shift_u1 */

void mult_spin_pseudoscalar_u1(field_offset src, field_offset dest ) {

 register int i;
 register site *s;
 complex temp;
 int p ; 
  
  if(dest==F_OFFSET(tempvec[2])||src==F_OFFSET(tempvec[2]))
    { 
      node0_printf("mult_spin_pseudoscalar_u1(): ERROR! dest or src should ");
      node0_printf("not be tempvec[2]\n" );
      terminate(101) ;
    }
  /*clean up dest */
  FORALLMYSITES(i,s){
    //clearvec((su3_vector *)F_PT(s,dest));
    (*(complex *)F_PT(s, dest)).real = (*(complex *)F_PT(s, dest)).imag = 0.0;
  }   
  for(p=0;p<24;p++)
     {
       eta_shift_u1(4,eps[p].d,src,F_OFFSET(tempvec[2])) ;
       /*  Multiply the the extra 1/24 needed by the            *
	* definition of the operator (number of permutations)   */
       FORALLMYSITES(i,s){
	 CMULREAL( s->tempvec[2], eps[p].sign/24.0, temp);
	 CSUM( (*(complex *)F_PT(s,dest)), temp);
	 /*scalar_mult_sum_su3_vector((su3_vector *)F_PT(s,dest), 
	   &(s->tempvec[2]), eps[p].sign/24.0);*/
       }
     } 

} /* mult_spin_pseudoscalar_u1 */
