/******* dslash_fn_u1.c - dslash for improved KS fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- improved.  This version for "fat plus
   Naik" quark action.  Connection to nearest neighbors stored in
   fatlink and to third nearest neighbors in longlink */

/* This version overlaps computation and gathers from negative
   directions, and has an extra lattice loop devoted to exclusively to
   sub_four_vectors (traditional algorithm) */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* C. DeTar 9/29/01 Standardized prefetching and synced versions */

#include "generic_ks_includes_u1.h"	/* definitions files and prototypes */
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

// 4/11/2010 added "_u1" to the end of function names

#define INDEX_3RD(dir) (dir - 8)      /* this gives the 'normal' direction */


/* Temporary work space for dslash_fn_field_special */ 
static complex *temp[9] ;
/* Flag indicating if temp is allocated               */
static int temp_not_allocated=1 ;

static void 
cleanup_one_gather_set(msg_tag *tags[])
{
  int i;

  FORALLMYUPDIR(i){
    cleanup_gather( tags[i] );
    cleanup_gather( tags[OPP_DIR(i)] );
  }

  FORALLMYUP_3_DIR(i){
    cleanup_gather( tags[i] );
    cleanup_gather( tags[OPP_3_DIR(i)] );
  }
}

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[])
{
  cleanup_one_gather_set(tags1);
  cleanup_one_gather_set(tags2);
}

void cleanup_dslash_temps(){
  register int i ;
  if(!temp_not_allocated)
    for(i=0;i<9;i++) {
      free(temp[i]) ; 
    }
  temp_not_allocated=1 ;
}



void dslash_fn_field_u1(complex *src, complex *dest, int parity,
		      ferm_links_t_u1 *fn) {
  msg_tag *tag[16];
    
   dslash_fn_field_special(src, dest, parity, tag, 1, fn);
   cleanup_one_gather_set(tag);
}

/* Special dslash for use by congrad.  Uses restart_gather_field() when
  possible. Next to last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_field, otherwise use restart_gather_field. 
  The calling program must clean up the gathers and temps! */
void dslash_fn_field_special_u1(complex *src, complex *dest,
			     int parity, msg_tag **tag, int start,
			     ferm_links_t_u1 *fn){
  register int i;
  register site *s;
  register int dir,otherparity=0;
  register complex *fat4, *long4;
  complex add_temp[4];
  complex *t_fatlink;
  complex *t_longlink;
  
  /* allocate temporary work space only if not already allocated */

  if(temp_not_allocated)
    {
      FORALLMYUPDIR(dir){
	temp[dir]  =(complex *)malloc(sites_on_node*sizeof(complex));
	temp[dir+4]=(complex *)malloc(sites_on_node*sizeof(complex));
      }
      temp[8]=(complex *)malloc(sites_on_node*sizeof(complex));
      temp_not_allocated = 0 ;
    }
  
  /* load fatlinks and longlinks */
  if(!fn->valid){
    printf("dslash_fn_field_special: invalid fn links!\n");
    terminate(1);
  }
  t_longlink = fn->lng;
  t_fatlink = fn->fat;

  switch(parity)
    {
    case EVEN:	otherparity=ODD; break;
    case ODD:	otherparity=EVEN; break;
    case EVENANDODD:	otherparity=EVENANDODD; break;
    }
  
  /* Start gathers from positive directions */
  /* And start the 3-step gather too */
  FORALLMYUPDIR(dir){
    if(start==1)
      {
	tag[dir] = start_gather_field( src, sizeof(complex), 
					   dir, parity,gen_pt[dir] );
	tag[DIR3(dir)] = start_gather_field(src, sizeof(complex),
						DIR3(dir),parity, 
						gen_pt[DIR3(dir)] );
      }
    else
      {
	restart_gather_field( src, sizeof(complex), 
				  dir, parity,gen_pt[dir], tag[dir]);
	restart_gather_field(src, sizeof(complex), DIR3(dir), parity, 
				 gen_pt[DIR3(dir)], tag[DIR3(dir)]);
      }
  }
  
  /* Multiply by adjoint matrix at other sites */
  /* Use fat link for single link transport */
  FORMYSITESANDPARITY( i, s, otherparity ){
    if( i < loopend-FETCH_UP ){
       fat4 = &(t_fatlink[4*(i+FETCH_UP)]);
       long4 = &(t_longlink[4*(i+FETCH_UP)]);
       prefetch_V(&(src[i+FETCH_UP]));
       prefetch_4MVVVV( 
		       fat4,
		       &(temp[0][i+FETCH_UP]),
		       &(temp[1][i+FETCH_UP]),
		       &(temp[2][i+FETCH_UP]),
		       &(temp[3][i+FETCH_UP]) );
       prefetch_4MVVVV( 
		       long4,
		       &(temp[4][i+FETCH_UP]),
		       &(temp[5][i+FETCH_UP]),
		       &(temp[6][i+FETCH_UP]),
		       &(temp[7][i+FETCH_UP]) );
    }

    fat4 = &(t_fatlink[4*i]);
    long4 = &(t_longlink[4*i]);
    //insert a routine for complex numbers that multiplies conjugate of fat(long) (NOTE: IS THERE SUPPOSED TO BE A POINTER HERE)
    CMULJ_(fat4[0], src[i], temp[0][i]);
    CMULJ_(fat4[1], src[i], temp[1][i]);
    CMULJ_(fat4[2], src[i], temp[2][i]); //z-direction
    CMULJ_(fat4[3], src[i], temp[3][i]);

    /* multiply by 3-link matrices too */
    CMULJ_(long4[0], src[i], temp[4][i]);
    CMULJ_(long4[1], src[i], temp[5][i]);
    CMULJ_(long4[2], src[i], temp[6][i]); //z-direction
    CMULJ_(long4[3], src[i], temp[7][i]);

  } END_LOOP
      
      /* Start gathers from negative directions */
   FORALLMYUPDIR(dir) {
    if (start==1) tag[OPP_DIR(dir)] = start_gather_field( temp[dir],
							  sizeof(complex), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
    else restart_gather_field( temp[dir], sizeof(complex), 
			       OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)], tag[OPP_DIR(dir)] );
  }
  
  /* Start 3-neighbour gathers from negative directions */
  FORALLMYUP_3_DIR(dir){
    if (start==1) tag[OPP_3_DIR(dir)]=start_gather_field(
		 temp[INDEX_3RD(dir)+4], sizeof(complex), 
		 OPP_3_DIR( dir), parity, gen_pt[OPP_3_DIR(dir)] );
    else restart_gather_field(temp[INDEX_3RD(dir)+4], 
	      sizeof(complex), OPP_3_DIR( dir),parity, 
	      gen_pt[OPP_3_DIR(dir)], tag[OPP_3_DIR(dir)] );
  }

  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  /* wait for the 3-neighbours from positive directions, multiply */
  FORALLMYUPDIR(dir) {
    wait_gather(tag[dir]);
    wait_gather(tag[DIR3(dir)]);
  }
  
  FORMYSITESANDPARITY(i,s,parity){
    if( i < loopend-FETCH_UP ){
      fat4 = &(t_fatlink[4*(i+FETCH_UP)]);
      long4 = &(t_longlink[4*(i+FETCH_UP)]);
      prefetch_4MVVVV( 
		      fat4,
		      (complex *)gen_pt[XUP][i+FETCH_UP],
		      (complex *)gen_pt[YUP][i+FETCH_UP],
		      (complex *)gen_pt[ZUP][i+FETCH_UP],
		      (complex *)gen_pt[TUP][i+FETCH_UP] );
      prefetch_4MVVVV( 
		      long4,
		      (complex *)gen_pt[X3UP][i+FETCH_UP],
		      (complex *)gen_pt[Y3UP][i+FETCH_UP],
		      (complex *)gen_pt[Z3UP][i+FETCH_UP],
		      (complex *)gen_pt[T3UP][i+FETCH_UP] );
      prefetch_VVVV( 
		    (complex *)gen_pt[XDOWN][i+FETCH_UP],
		    (complex *)gen_pt[YDOWN][i+FETCH_UP],
		    (complex *)gen_pt[ZDOWN][i+FETCH_UP],
		    (complex *)gen_pt[TDOWN][i+FETCH_UP] );
      prefetch_VVVV( 
		    (complex *)gen_pt[X3DOWN][i+FETCH_UP],
		    (complex *)gen_pt[Y3DOWN][i+FETCH_UP],
		    (complex *)gen_pt[Z3DOWN][i+FETCH_UP],
		    (complex *)gen_pt[T3DOWN][i+FETCH_UP] );
    }
    
    fat4 = &(t_fatlink[4*i]);
    long4 = &(t_longlink[4*i]);
    //multiply each fat link by field and store in temporary array then add and store in dest[i]
    CMUL(fat4[0],  (complex *)gen_pt[XUP][i], add_temp[0]);
    CMUL(fat4[1],  (complex *)gen_pt[YUP][i], add_temp[1]);
    CMUL(fat4[2],  (complex *)gen_pt[ZUP][i], add_temp[2]); //z-direction
    CMUL(fat4[3],  (complex *)gen_pt[TUP][i], add_temp[3]);
    for(int j=0; j<4; j++) if(j!=2) dest[i] += add_temp[j]; //do not add in z-direction
  
    //now do the same for the long links and store in temp[8][i]
    CMUL(long4[0],  (complex *)gen_pt[X3UP][i], add_temp[0]);
    CMUL(long4[1],  (complex *)gen_pt[Y3UP][i], add_temp[1]);
    CMUL(long4[2],  (complex *)gen_pt[Z3UP][i], add_temp[2]); //z-direction
    CMUL(long4[3],  (complex *)gen_pt[T3UP][i], add_temp[3]);
    for(int j=0; j<4; j++) if(j!=2) temp[8][i] += add_temp[j];

  } END_LOOP
   
  /* Wait gathers from negative directions, accumulate (negative) */
  /* and the same for the negative 3-rd neighbours */
  for(dir=XUP; dir<=TUP; dir++){
    wait_gather(tag[OPP_DIR(dir)]);
  }
  for(dir=X3UP; dir<=T3UP; dir++){
    wait_gather(tag[OPP_3_DIR(dir)]);
  }
  
  FORMYSITESANDPARITY(i,s,parity){
    if( i < loopend-FETCH_UP ){
      prefetch_VVVVV( 
		     &(dest[i+FETCH_UP]),
		     (complex *)gen_pt[XDOWN][i+FETCH_UP],
		     (complex *)gen_pt[YDOWN][i+FETCH_UP],
		     (complex *)gen_pt[ZDOWN][i+FETCH_UP],
		     (complex *)gen_pt[TDOWN][i+FETCH_UP] );
      prefetch_VVVVV( 
		     &(temp[8][i+FETCH_UP]), 
		     (complex *)gen_pt[X3DOWN][i+FETCH_UP],
		     (complex *)gen_pt[Y3DOWN][i+FETCH_UP],
		     (complex *)gen_pt[Z3DOWN][i+FETCH_UP],
		     (complex *)gen_pt[T3DOWN][i+FETCH_UP] );
    }

    //subract all the 1st down link matrices and store in dest[i]
    CSUB(dest[i], gen_pt[XDOWN][i], dest[i]);
    CSUB(dest[i], gen_pt[YDOWN][i], dest[i]);
    CSUB(dest[i], gen_pt[ZDOWN][i], dest[i]); //z-direction
    CSUB(dest[i], gen_pt[TDOWN][i], dest[i]);

    //subract all the 3rd down link matrices and store in temp[8][i]
    CSUB(temp[8][i], gen_pt[XDOWN][i], temp[8][i]);
    CSUB(temp[8][i], gen_pt[YDOWN][i], temp[8][i]);
    CSUB(temp[8][i], gen_pt[ZDOWN][i], temp[8][i]); //z-direction
    CSUB(temp[8][i], gen_pt[TDOWN][i], temp[8][i]);

    /* Now need to add these things together */
    CADD(dest[i], temp[8][i], dest[i]);

  } END_LOOP 
      
}



