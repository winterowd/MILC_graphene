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
#include "../graphene/lattice.h"
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

void dslash_fn_site_u1( field_offset src, field_offset dest, int parity,
		     ferm_links_u1_t *fn )
{
  msg_tag *tag[16];
  
  dslash_fn_site_special_u1(src, dest, parity, tag, 1, fn );
  cleanup_one_gather_set(tag);
}

void dslash_fn_site_special_u1( field_offset src, field_offset dest,
			     int parity, msg_tag **tag, int start,
			     ferm_links_u1_t *fn){
  register int i, j;
  register site *s;
  register int dir,otherparity=0;
  register complex *fat4, *long4;
  complex *temp;
  complex add_temp[4];
  complex *t_fatlink;
  complex *t_longlink;

  if(!fn->valid){
    printf("dslash_fn_site_special: invalid fn links!\n");
    terminate(1);
  }
  t_longlink = fn->lng;
  t_fatlink = fn->fat;

  switch(parity){
  case EVEN:	otherparity=ODD; break;
  case ODD:	otherparity=EVEN; break;
  case EVENANDODD:	otherparity=EVENANDODD; break;
  }
  
  /* Start gathers from positive directions */
  FORALLMYUPDIR(dir) {
    if(start==1) tag[dir] = start_gather_site( src, sizeof(complex),
					       dir, parity, gen_pt[dir] );
    else restart_gather_site( src, sizeof(complex),
			      dir, parity, gen_pt[dir] , tag[dir] );
  }
  
  /* and start the 3rd neighbor gather */
  FORALLMYUP_3_DIR(dir) {
    if(start==1) tag[dir] = start_gather_site( src, sizeof(complex),
					       dir, parity, gen_pt[dir] );
    else restart_gather_site( src, sizeof(complex),
			      dir, parity, gen_pt[dir] , tag[dir] ); 
  }
  
  /* Multiply by adjoint matrix at other sites */
  FORMYSITESANDPARITY(i,s,otherparity){
    if( i < loopend-FETCH_UP ){
      fat4 = &(t_fatlink[4*(i+FETCH_UP)]);
      long4 = &(t_longlink[4*(i+FETCH_UP)]);
      prefetch_4MV4V( 
		     fat4,
		     (complex *)F_PT(s+FETCH_UP,src),
		     (s+FETCH_UP)->tempvec );
      prefetch_4MV4V(
		     long4,
		     (complex *)F_PT(s+FETCH_UP,src),
		     (s+FETCH_UP)->templongvec );
    }
    
    fat4 = &(t_fatlink[4*i]);
    long4 = &(t_longlink[4*i]);
    /*mult_adj_su3_mat_vec_4dir( fat4,
      (su3_vector *)F_PT(s,src), s->tempvec ); */
    CMULJ_( fat4[0], (*(complex *)F_PT(s,src)), s->tempvec[0]);
    CMULJ_( fat4[1], (*(complex *)F_PT(s,src)), s->tempvec[1]);
#ifdef FOUR_DIM
    CMULJ_( fat4[2], (*(complex *)F_PT(s,src)), s->tempvec[2]);
#endif
    CMULJ_( fat4[3], (*(complex *)F_PT(s,src)), s->tempvec[3]);
    /* multiply by 3-link matrices too */
    /*mult_adj_su3_mat_vec_4dir( long4,
      (su3_vector *)F_PT(s,src), s->templongvec );*/
    CMULJ_( long4[0], (*(complex *)F_PT(s,src)), s->templongvec[0]);
    CMULJ_( long4[1], (*(complex *)F_PT(s,src)), s->templongvec[1]);
#ifdef FOUR_DIM
    CMULJ_( long4[2], (*(complex *)F_PT(s,src)), s->templongvec[2]);
#endif
    CMULJ_( long4[3], (*(complex *)F_PT(s,src)), s->templongvec[3]);
  } END_LOOP
      
      /* Start gathers from negative directions */
      FORALLMYUPDIR(dir) {
    if (start==1) tag[OPP_DIR(dir)] = start_gather_site( F_OFFSET(tempvec[dir]),
							 sizeof(complex), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
    else restart_gather_site( F_OFFSET(tempvec[dir]), sizeof(complex),
			      OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] , tag[OPP_DIR(dir)] );
  }
  
  /* and 3rd neighbours */
  FORALLMYUP_3_DIR(dir) {
    /**printf("dslash_fn_site_special: down gathers, start=%d\n",start);**/
    if (start==1) tag[OPP_3_DIR(dir)] = 
		    start_gather_site( F_OFFSET(templongvec[INDEX_3RD(dir)]),
				       sizeof(complex), OPP_3_DIR(dir), parity, gen_pt[OPP_3_DIR(dir)] );
    else restart_gather_site( F_OFFSET(templongvec[INDEX_3RD(dir)]),
			      sizeof(complex), OPP_3_DIR( dir), parity, gen_pt[OPP_3_DIR(dir)],
			      tag[OPP_3_DIR(dir)] );
  }
  
  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  FORALLMYUPDIR(dir){
    wait_gather(tag[dir]);
  }
  
  /* wait for the 3-neighbours from positive directions, multiply */
  FORALLMYUP_3_DIR(dir) {
    wait_gather(tag[dir]);
  }
  FORMYSITESANDPARITY(i,s,parity){
    if( i < loopend-FETCH_UP ){
      fat4 = &(t_fatlink[4*(i+FETCH_UP)]);
      long4 = &(t_longlink[4*(i+FETCH_UP)]);
      prefetch_VV(
		  (complex *)F_PT(s+FETCH_UP,dest),
		  (complex *) &((s+FETCH_UP)->templongv1));
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
    }
    fat4 = &(t_fatlink[4*i]);
    long4 = &(t_longlink[4*i]);
    /*mult_su3_mat_vec_sum_4dir( fat4,
	       (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
	       (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
	       (su3_vector *)F_PT(s,dest)); */
    CMUL( fat4[0], (*(complex *)(gen_pt[XUP][i])), add_temp[0] );
    CMULREAL( add_temp[0], v_Fermi, add_temp[0] );
    CMUL( fat4[1], (*(complex *)(gen_pt[YUP][i])), add_temp[1] );
    CMULREAL( add_temp[1], v_Fermi, add_temp[1] );
#ifdef FOUR_DIM
    CMUL( fat4[2], (*(complex *)(gen_pt[ZUP][i])), add_temp[2] );
#endif
    CMUL( fat4[3], (*(complex *)(gen_pt[TUP][i])), add_temp[3] );

    /* clear destination vector then sum and put result in dest[i] */
    temp = (complex *)F_PT(s,dest);
    temp->real = 0.0;
    temp->imag = 0.0;

    for(j=0; j<4; j++) { 
#ifndef FOUR_DIM
      if(j != 2) 
#endif
	CADD( add_temp[j], (*(complex *)F_PT(s,dest)), (*(complex *)F_PT(s,dest)) );
    }
    /*mult_su3_mat_vec_sum_4dir( long4,
	       (su3_vector *)gen_pt[X3UP][i], (su3_vector *)gen_pt[Y3UP][i],
	       (su3_vector *)gen_pt[Z3UP][i], (su3_vector *)gen_pt[T3UP][i],
	       (su3_vector *) &(s->templongv1)); */
    /* do the same with templongv1 */
    s->templongv1.real = 0.0;
    s->templongv1.imag = 0.0;
    CMUL( long4[0], (*(complex *)(gen_pt[X3UP][i])), add_temp[0] );
    CMULREAL( add_temp[0], v_Fermi, add_temp[0] );
    CMUL( long4[1], (*(complex *)(gen_pt[Y3UP][i])), add_temp[1] );
    CMULREAL( add_temp[1], v_Fermi, add_temp[1] );
#ifdef FOUR_DIM
    CMUL( long4[2], (*(complex *)(gen_pt[Z3UP][i])), add_temp[2] );
#endif
    CMUL( long4[3], (*(complex *)(gen_pt[T3UP][i])), add_temp[3] );
    for(j=0; j<4; j++) { 
#ifndef FOUR_DIM
      if(j != 2)
#endif
	CADD( add_temp[j], s->templongv1, s->templongv1 );
    }
  } END_LOOP
      

  /* Wait gathers from negative directions, accumulate (negative) */
  FORALLMYUPDIR(dir) {
    wait_gather(tag[OPP_DIR(dir)]);
  }
  
  /* and the same for the negative 3-rd neighbours */
  
  FORALLMYUP_3_DIR(dir) {
    wait_gather(tag[OPP_3_DIR(dir)]);
  }
  
  FORMYSITESANDPARITY(i,s,parity){
    if( i < loopend-FETCH_UP ){
      prefetch_VV(
		  (complex *)F_PT(s+FETCH_UP,dest),
		  (complex *) &((s+FETCH_UP)->templongv1));
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
    CMULREAL(*(complex *)(gen_pt[XDOWN][i]), v_Fermi, add_temp[0]);
    CSUB( (*(complex *)F_PT(s,dest)), add_temp[0], 
	  (*(complex *)F_PT(s,dest)) );
    CMULREAL(*(complex *)(gen_pt[YDOWN][i]), v_Fermi, add_temp[1]);
    CSUB( (*(complex *)F_PT(s,dest)), add_temp[1], 
	  (*(complex *)F_PT(s,dest)) );
#ifdef FOUR_DIM
    CSUB( (*(complex *)F_PT(s,dest)), *(complex *)(gen_pt[ZDOWN][i]), 
      (*(complex *)F_PT(s,dest)) );
#endif
    CSUB( (*(complex *)F_PT(s,dest)), *(complex *)(gen_pt[TDOWN][i]), 
	  (*(complex *)F_PT(s,dest)) );
    /* sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
		       (su3_vector *)(gen_pt[XDOWN][i]),
		       (su3_vector *)(gen_pt[YDOWN][i]),
		       (su3_vector *)(gen_pt[ZDOWN][i]),
		       (su3_vector *)(gen_pt[TDOWN][i]) ); */
    CMULREAL(*(complex *)(gen_pt[X3DOWN][i]), v_Fermi, add_temp[0]);
    CSUB( s->templongv1, add_temp[0], s->templongv1 );
    CMULREAL(*(complex *)(gen_pt[Y3DOWN][i]), v_Fermi, add_temp[1]);
    CSUB( s->templongv1, add_temp[1], s->templongv1 );
#ifdef FOUR_DIM
    CSUB( s->templongv1, *(complex *)(gen_pt[Z3DOWN][i]), s->templongv1 ); 
#endif
    CSUB( s->templongv1, *(complex *)(gen_pt[T3DOWN][i]), s->templongv1 );
    /* sub_four_su3_vecs( & (s->templongv1), 
		       (su3_vector *)(gen_pt[X3DOWN][i]),
		       (su3_vector *)(gen_pt[Y3DOWN][i]),
		       (su3_vector *)(gen_pt[Z3DOWN][i]),
		       (su3_vector *)(gen_pt[T3DOWN][i]) ); */
    /*** Now need to add these things together ***/
    CADD( (*(complex *)F_PT(s,dest)), s->templongv1, 
	  (*(complex *)F_PT(s,dest)) );
    /* add_su3_vector((su3_vector *)F_PT(s,dest), &(s->templongv1),
       (su3_vector *)F_PT(s,dest)); */
  } END_LOOP
      

}


void dslash_fn_field_u1(complex *src, complex *dest, int parity,
		      ferm_links_u1_t *fn) {
  msg_tag *tag[16];
    
  dslash_fn_field_special_u1(src, dest, parity, tag, 1, fn);
  cleanup_one_gather_set(tag);
}

/* Special dslash for use by congrad.  Uses restart_gather_field() when
  possible. Next to last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_field, otherwise use restart_gather_field. 
  The calling program must clean up the gathers and temps! */
void dslash_fn_field_special_u1(complex *src, complex *dest,
			     int parity, msg_tag **tag, int start,
			     ferm_links_u1_t *fn){
  register int i;
  register site *s;
  register int dir,otherparity=0;
  register complex *fat4, *long4;
  complex add_temp[4];
  complex *t_fatlink;
  complex *t_longlink;
  int j;
  
  for(j=0; j<4; j++) {
    add_temp[j].real = 0.0;
    add_temp[j].imag = 0.0;
  }
  /* allocate temporary work space only if not already allocated */
  if(temp_not_allocated)
    {
      FORALLMYUPDIR(dir){
	//printf("temp1 dslash\n");
	temp[dir]  =(complex *)malloc(sites_on_node*sizeof(complex));
	temp[dir+4]=(complex *)malloc(sites_on_node*sizeof(complex));
	//printf("temp2 dslash\n");
      }
      //printf("out of updir dslash\n");
      temp[8]=(complex *)malloc(sites_on_node*sizeof(complex));
      temp_not_allocated = 0 ;
      //printf("out of updir dslash2\n");
    }
  //printf("out of not allocated dslash2\n");
  /* load fatlinks and longlinks */
  if(!fn->valid){
    printf("dslash_fn_field_special: invalid fn links!\n");
    terminate(1);
  }
  t_longlink = fn->lng;
  t_fatlink = fn->fat;
  //printf("loaded links dslash\n");
  switch(parity)
    {
    case EVEN:	otherparity=ODD; break;
    case ODD:	otherparity=EVEN; break;
    case EVENANDODD:	otherparity=EVENANDODD; break;
    }
  
  /* Start gathers from positive directions */
  /* And start the 3-step gather too */
  FORALLMYUPDIR(dir){
    //printf("inside gathers with dir %d\n", dir);
    if(start==1)
      {
	//printf("inside gathers 2\n");
	tag[dir] = start_gather_field( src, sizeof(complex), 
					   dir, parity,gen_pt[dir] );
	//printf("inside gathers 3\n");
	tag[DIR3(dir)] = start_gather_field(src, sizeof(complex),
						DIR3(dir),parity, 
						gen_pt[DIR3(dir)] );
	//printf("inside gathers 4\n");
      }
    else
      {
	//printf("inside gathers 5\n");
	restart_gather_field( src, sizeof(complex), 
				  dir, parity,gen_pt[dir], tag[dir]);
	restart_gather_field(src, sizeof(complex), DIR3(dir), parity, 
				 gen_pt[DIR3(dir)], tag[DIR3(dir)]);
      }
  }
  //printf("first gathers dslash\n");
  /* Multiply by adjoint matrix at other sites */
  /* Use fat link for single link transport */
  FORMYSITESANDPARITY( i, s, otherparity ){
    if( i < loopend-FETCH_UP ){
       fat4 = &(t_fatlink[4*(i+FETCH_UP)]);
       long4 = &(t_longlink[4*(i+FETCH_UP)]);
       //printf("fat4 and long4 initialized\n");
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
       //printf("prefetch accomplished\n");
    }

    fat4 = &(t_fatlink[4*i]);
    long4 = &(t_longlink[4*i]);
    //printf("fat4 and long4 initialized again\n");
    //insert a routine for complex numbers that multiplies conjugate of fat(long) (NOTE: IS THERE SUPPOSED TO BE A POINTER HERE)
    CMULJ_(fat4[0], src[i], temp[0][i]);
    CMULJ_(fat4[1], src[i], temp[1][i]);
#ifdef FOUR_DIM
    CMULJ_(fat4[2], src[i], temp[2][i]); //z-direction
#endif
    CMULJ_(fat4[3], src[i], temp[3][i]);
    //printf("multiplied source by fat links\n");
    /* multiply by 3-link matrices too */
    CMULJ_(long4[0], src[i], temp[4][i]);
    CMULJ_(long4[1], src[i], temp[5][i]);
#ifdef FOUR_DIM
    CMULJ_(long4[2], src[i], temp[6][i]); //z-direction
#endif
    CMULJ_(long4[3], src[i], temp[7][i]);
    //printf("multiplied source by long links\n");
  } END_LOOP
      //printf("out of loop sites and parity\n");
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
    CMUL(fat4[0],  *(complex *)gen_pt[XUP][i], add_temp[0]);
    CMULREAL(add_temp[0], v_Fermi, add_temp[0]); 
    CMUL(fat4[1],  *(complex *)gen_pt[YUP][i], add_temp[1]);
    CMULREAL(add_temp[1], v_Fermi, add_temp[1]); 
#ifdef FOUR_DIM
    CMUL(fat4[2],  *(complex *)gen_pt[ZUP][i], add_temp[2]); 
#endif
    CMUL(fat4[3],  *(complex *)gen_pt[TUP][i], add_temp[3]);
    
    //reset before adding to dest
    dest[i].real = 0.0; dest[i].imag = 0.0;
    for(j=0; j<4; j++) { 
#ifndef FOUR_DIM
      if(j!=2) //do not add in z-direction for graphene
#endif
      CADD(dest[i], add_temp[j], dest[i]); } //do not add in z-direction
    
    //now do the same for the long links and store in temp[8][i]
    CMUL(long4[0],  *(complex *)gen_pt[X3UP][i], add_temp[0]);
    CMULREAL(add_temp[0], v_Fermi, add_temp[0])
    CMUL(long4[1],  *(complex *)gen_pt[Y3UP][i], add_temp[1]);
    CMULREAL(add_temp[1], v_Fermi, add_temp[1]);
#ifdef FOUR_DIM
    CMUL(long4[2],  *(complex *)gen_pt[Z3UP][i], add_temp[2]); //z-direction
#endif
    CMUL(long4[3],  *(complex *)gen_pt[T3UP][i], add_temp[3]);
    temp[8][i].real = 0.0; temp[8][i].imag = 0.0;
    for(j=0; j<4; j++) { 
#ifndef FOUR_DIM
      if(j!=2) 
#endif
	CADD(temp[8][i], add_temp[j], temp[8][i]); }

  } END_LOOP
  
  /* Wait gathers from negative directions, accumulate (negative) */
  /* and the same for the negative 3-rd neighbours */
  FORALLMYUPDIR(dir){
    wait_gather(tag[OPP_DIR(dir)]);
  }
  FORALLMYUP_3_DIR(dir) {
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
    CMULREAL(*(complex *)gen_pt[XDOWN][i], v_Fermi, add_temp[0]);
    CSUB(dest[i], add_temp[0], dest[i]);
    //CSUB(dest[i], *(complex *)gen_pt[XDOWN][i], dest[i]);
    CMULREAL(*(complex *)gen_pt[YDOWN][i], v_Fermi, add_temp[1]);
    CSUB(dest[i], add_temp[1], dest[i]);
    //CSUB(dest[i], *(complex *)gen_pt[YDOWN][i], dest[i]);
#ifdef FOUR_DIM
    CSUB(dest[i], *(complex *)gen_pt[ZDOWN][i], dest[i]); //z-direction
#endif
    CSUB(dest[i], *(complex *)gen_pt[TDOWN][i], dest[i]);

    //subract all the 3rd down link matrices and store in temp[8][i]
    CMULREAL(*(complex *)gen_pt[X3DOWN][i], v_Fermi, add_temp[0]);
    CSUB(temp[8][i], add_temp[0], temp[8][i]);
    //CSUB(temp[8][i], *(complex *)gen_pt[X3DOWN][i], temp[8][i]);
    CMULREAL(*(complex *)gen_pt[Y3DOWN][i], v_Fermi, add_temp[1]);
    CSUB(temp[8][i], add_temp[1], temp[8][i]);
    //CSUB(temp[8][i], *(complex *)gen_pt[Y3DOWN][i], temp[8][i]);
#ifdef FOUR_DIM
    CSUB(temp[8][i], *(complex *)gen_pt[Z3DOWN][i], temp[8][i]); //z-direction
#endif
    CSUB(temp[8][i], *(complex *)gen_pt[T3DOWN][i], temp[8][i]);

    /* Now need to add these things together */
    CADD(dest[i], temp[8][i], dest[i]);
   
  } END_LOOP 
      /*printf("end of dslash_fn_field_special_u1\n");
  FORALLMYSITES(i,s){
    printf("At site (%d, %d, %d, %d), index %d\n", s->x, s->y, s->z, s->t, i);
    printf("Destination vector: Real %f, Imag %f\n", dest[i].real, 
	   dest[i].imag);
    printf("Source vector: Real %f, Imag %f\n", src[i].real, 
	   src[i].imag);
  }
      */
}



