/****** gauge_stuff_u1.c  -- ******************/
/* MIMD version 7 */
/* gauge action stuff for improved action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file 
* C.W. modified 10/10, U(1) improved gauge action calculation */

/**#define GFTIME**/ /* For timing gauge force calculation */
#include "generic_includes.h"	/* definitions files and prototypes */

#ifdef LOOPEND
#undef FORALLSITES
#define FORALLSITES(i,s) \
{ register int loopend; loopend=sites_on_node; \
for( i=0,  s=lattice ; i<loopend; i++,s++ )
#define END_LOOP }
#else
#define END_LOOP        /* define it to be nothing */
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)
void printpath( int *path, int length );

#define GAUGE_ACTION_PART1
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#include <gauge_action.h>
#undef GAUGE_ACTION_PART1

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

char gauge_action_description[128];
int  gauge_action_nloops=NLOOP;
int  gauge_action_nreps=NREPS;
static int loop_length[NLOOP];	/* lengths of various kinds of loops */
static int loop_num[NLOOP];	/* number of rotations/reflections  for each kind */

    /* table of directions, for each rotation and reflection of each kind of
	loop.  tabulated with "canonical" starting point and direction. */
static int*** loop_table;
    /* table of coefficients in action, for various "representations" (actually,
	powers of the trace) */
static Real **loop_coeff;
    /* for each rotation/reflection, an integer distinct for each starting
	point, or each cyclic permutation of the links */
int loop_char[MAX_NUM];

static void char_num( int *dig, int *chr, int length);

/* Make table of loops in action */
void make_loop_table() {

    int perm[8],pp[8],ir[4];
    int length,iloop,i,j,chr;
    int vec[MAX_LENGTH];
    int count,flag;
    int total_dyn_flavors;
    char myname[] = "make_loop_table";

    total_dyn_flavors = 0;
    for(i = 0; i < n_dyn_masses; i++){
      total_dyn_flavors += dyn_flavors[i];
    }

    /* Allocate as if loop_table[NLOOP][MAX_NUM][MAX_LENGTH] */

    loop_table = (int ***)malloc(sizeof(int **)*NLOOP);
    if(loop_table == NULL){
      printf("%s(%d): No room for loop_table\n",myname,this_node);
      terminate(1);
    }

    for(iloop = 0; iloop < NLOOP; iloop++){
      loop_table[iloop] = (int **)malloc(sizeof(int *)*MAX_NUM);
      if(loop_table[iloop] == NULL){
	printf("%s(%d): No room for loop_table\n",myname,this_node);
	terminate(1);
      }

      for(count = 0; count < MAX_NUM; count++){
	loop_table[iloop][count] = (int *)malloc(sizeof(int)*MAX_LENGTH);
	if(loop_table[iloop][count] == NULL){
	  printf("%s(%d): No room for loop_table\n",myname,this_node);
	  terminate(1);
	}
      }
    }

    /* Allocate as if loop_coeff[NLOOP][NREPS] */

    loop_coeff = (Real **)malloc(sizeof(Real *)*NLOOP);
    if(loop_coeff == NULL){
      printf("%s(%d): No room for loop_coeff\n",myname,this_node);
      terminate(1);
    }

    for(iloop = 0; iloop < NLOOP; iloop++){
      loop_coeff[iloop] = (Real *)malloc(sizeof(Real)*NREPS);
      if(loop_coeff[iloop] == NULL){
	printf("%s(%d): No room for loop_coeff\n",myname,this_node);
	terminate(1);
      }
    }

#define GAUGE_ACTION_PART2
/* defines all loops and their coefficients */
#include <gauge_action.h>
#undef GAUGE_ACTION_PART2

    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	count=0;
	/* permutations */
	for(perm[0]=0;perm[0]<4;perm[0]++)
	for(perm[1]=0;perm[1]<4;perm[1]++)
	for(perm[2]=0;perm[2]<4;perm[2]++)
	for(perm[3]=0;perm[3]<4;perm[3]++){
	    if(perm[0] != perm[1] && perm[0] != perm[2] 
		&& perm[0] != perm[3] && perm[1] != perm[2]
	    	&& perm[1] != perm[3] && perm[2] != perm[3] ) {
	        /* reflections*/
	 	for(ir[0]=0;ir[0]<2;ir[0]++)
		for(ir[1]=0;ir[1]<2;ir[1]++)
		for(ir[2]=0;ir[2]<2;ir[2]++)
		for(ir[3]=0;ir[3]<2;ir[3]++){
		    for(j=0;j<4;j++){
			pp[j]=perm[j];

			if(ir[j] == 1) pp[j]=7-pp[j];
			pp[7-j]=7-pp[j];
		    }
		    /* create new vector*/
		    for(j=0;j<length;j++) vec[j]=pp[loop_ind[iloop][j]];

         	    char_num(vec,&chr,length);
         	    flag=0;
		    /* check if it's a new set: */
		    for(j=0;j<count;j++) if(chr == loop_char[j])flag=1;
		    if(flag == 0 ){
			loop_char[count]=chr;
			for(j=0;j<length;j++)
			    loop_table[iloop][count][j]=vec[j];
			count++;
/**node0_printf("ADD LOOP: "); printpath( vec, length );**/
		    }
		    if(count>MAX_NUM){
			node0_printf("OOPS: MAX_NUM too small\n");
			exit(0);
		    }
		    loop_num[iloop]=count;

		} /* end reflection*/
	    } /* end permutation if block */
	} /* end permutation */
    } /* end iloop */

    /* print out the loop coefficients */
    node0_printf("gauge_action: total_dyn_flavors = %d\n",total_dyn_flavors);
    node0_printf("loop coefficients: nloop rep loop_coeff  multiplicity\n");
    for(i=0;i<NREPS;i++) for(j=0;j<NLOOP;j++) {
	node0_printf("                    %d %d      %e     %d\n",
	    j,i,loop_coeff[j][i],loop_num[j]);
    }

} /* make_loop_table */

int get_max_length(){
  return MAX_LENGTH;
}

int get_nloop(){
  return NLOOP;
}

int get_nreps(){
  return NREPS;
}

int *get_loop_length(){
  return loop_length;
}

int *get_loop_num(){
  return loop_num;
}

int ***get_loop_table(){
  return loop_table;
}

Real **get_loop_coeff(){
  return loop_coeff;
}

/* find a number uniquely identifying the cyclic permutation of a path,
   or the starting point on the path.  Backwards paths are considered
   equivalent here, so scan those too. */
static void char_num( int *dig, int *chr, int length){
    int j;
    int bdig[MAX_LENGTH],tenl,newv,old;
    /* "dig" is array of directions.  "bdig" is array of directions for
	backwards path. */

  tenl=1;
  for(j=0;j<length-1;j++) tenl=tenl*10;

  *chr=dig[length-1];
  for(j=length-2;j>=0;j--) *chr= *chr*10+dig[j];

  /* forward*/
  old=*chr;
  for(j=length-1;j>=1;j--){
       newv=old-tenl*dig[j];
       newv=newv*10+dig[j];
       if(newv < *chr) *chr=newv;
       old=newv;           }

   /* backward*/
   for(j=0;j<length;j++)bdig[j]=7-dig[length-j-1];
   old=bdig[length-1];
   for(j=length-2;j>=0;j--) old=old*10+bdig[j];
   if(old < *chr ) *chr=old;
   for(j=length-1;j>=1;j--){
       newv=old-tenl*bdig[j];
       newv=newv*10+bdig[j];
       if(newv < *chr) *chr=newv;
       old=newv;           }

} /* char_num */

double imp_gauge_action_u1() {
    register int i;
    int rep;
    register site *s;
    complex trace;
    double g_action;
    double action,act2,total_action;
    complex *tempmat1;
    int length;

    /* these are for loop_table  */
    int ln,iloop;

    g_action=0.0;

    tempmat1 = (complex *)special_alloc(sites_on_node*sizeof(complex));
    if(tempmat1 == NULL){
      printf("imp_gauge_action: Can't malloc temporary\n");
      terminate(1);
    }

    /* gauge action */
    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	/* loop over rotations and reflections */
	for(ln=0;ln<loop_num[iloop];ln++){

	    path_product_u1( loop_table[iloop][ln] , length, tempmat1 );

	    FORALLSITES(i,s){
	      //trace=trace_su3( &tempmat1[i] );
	      action =  1.0 - (double)tempmat1[i].real;
	      /* need the "3 -" for higher characters */
	      total_action= (double)loop_coeff[iloop][0]*action;
	      act2=action;
	      for(rep=1;rep<NREPS;rep++){
		act2 *= action;
		total_action += (double)loop_coeff[iloop][rep]*act2;
	      }
	      
	      g_action  += total_action;

	    } END_LOOP /* sites */
	} /* ln */
    } /* iloop */

    g_doublesum( &g_action );
    special_free(tempmat1);
    return( g_action );
} /* imp_gauge_action_u1 */

double gauge_action_nc() {
Real *temp_real;
register int i,dir1,dir2;
register site *s;
register Real *m1,*m4;
Real mtmp, temp_sum;
double ss_sum,st_sum;
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;

    temp_real = (Real *)malloc(sizeof(Real)*sites_on_node);
    if(temp_real == NULL)
      {
	printf("gauge_action_nc: can't malloc temp_complex\n");
	fflush(stdout); terminate(1);
      }
    //#ifdef FOUR_DIM
    for(dir1=YUP;dir1<=TUP;dir1++){
      /*#else
	for(dir1=TUP;dir1<=TUP;dir1++){
	#endif*/
      for(dir2=XUP;dir2<dir1;dir2++){

	    mtag0 = start_gather_site( F_OFFSET(potential[dir2]), sizeof(Real),
		dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather_site( F_OFFSET(potential[dir1]), sizeof(Real),
		dir2, EVENANDODD, gen_pt[1] );

	    FORALLSITES(i,s){
		m1 = &(s->potential[dir1]);
		m4 = &(s->potential[dir2]);
		temp_real[i] = *m1 - *m4;
		//CMULJ_( *m4, *m1, temp_complex[i] );
		//mult_su3_an(m4,m1,&su3mat[i]);
	    }

	    wait_gather(mtag0);
	    wait_gather(mtag1);

	    FORALLSITES(i,s){
	      mtmp = temp_real[i] + *(Real *)(gen_pt[0][i]);
	      //CMUL( temp_complex[i], *(complex *)(gen_pt[0][i]), mtmp );
	      /* mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
		 &mtmp); */
	      
	      if(dir1==TUP ) {
		temp_sum = mtmp - *(Real *)(gen_pt[1][i]);
		//CMULJ_( *(complex *)(gen_pt[1][i]), mtmp, temp_mul );
		temp_sum = temp_sum*temp_sum;
		st_sum += (double)temp_sum;
		/*st_sum += (double)
		  realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp); */
	      }
	      else {          
		temp_sum = mtmp - *(Real *)(gen_pt[1][i]);
		//CMULJ_( *(complex *)(gen_pt[1][i]), mtmp, temp_mul );
		temp_sum = temp_sum*temp_sum;
		ss_sum += (double)temp_sum;
		/* ss_sum += (double)
		   realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp); */
	      }
	    }
	    
	    cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	}
    }
    g_doublesum( &ss_sum );
    g_doublesum( &st_sum );

    free(temp_real);

    return(ss_sum + st_sum);
}/* gauge_action_nc() */


/* Measure gauge observables:
    Loops in action (time and space directions treated differently)
    Polyakov loop

*/

void g_measure_u1( ){
    double ss_plaquette, st_plaquette;
    complex p_loop;
    register int i;
    register site *s;
    complex trace;
    double average[NREPS],action,act2,total_action;
    int length;
    complex *tempmat1;
    /* these are for loop_table  */
    int ln,iloop,rep;

    tempmat1 = (complex *)special_alloc(sites_on_node*sizeof(complex));
    if(tempmat1 == NULL){
      printf("g_measure: Can't malloc temporary\n");
      terminate(1);
    }

    /* KS and BC minus signs should be out for this routine */
    d_plaquette_u1( &ss_plaquette, &st_plaquette );
    if(this_node==0)printf("PLAQ:\t%f\t%f\n", ss_plaquette, st_plaquette );

    p_loop = ploop_u1();
    if(this_node==0)printf("P_LOOP:\t%e\t%e\n", p_loop.real, p_loop.imag );

    /* gauge action, all loops that contribute */
    total_action=0.0;
    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	/* loop over rotations and reflections */
	for(ln=0;ln<loop_num[iloop];ln++){

	    path_product_u1( loop_table[iloop][ln] , length, tempmat1 );

	    for(rep=0;rep<NREPS;rep++)average[rep] = 0.0;
	    FORALLSITES(i,s){
	      average[0] += (double)tempmat1[i].real;
	      action = 1.0 - (double)tempmat1[i].real;
	      total_action += (double)loop_coeff[iloop][0]*action;
	      /* need the "3 -" for higher characters */
	      act2=action;
		for(rep=1;rep<NREPS;rep++){
		    act2 *= action;
		    average[rep] += act2;
		    total_action += (double)loop_coeff[iloop][rep]*act2;
		} /* reps */
	    } END_LOOP /* sites */
	    g_vecdoublesum( average, NREPS );
	    /* dump the loop */
	    node0_printf("G_LOOP:  %d  %d  %d   ",iloop,ln,length);
	    for(rep=0;rep<NREPS;rep++)node0_printf("\t%e",average[rep]/volume);
	    node0_printf("\t( ");
	    for(i=0;i<length;i++)node0_printf("%d ",loop_table[iloop][ln][i]);
	    node0_printf(" )\n");
	} /* ln */
    } /* iloop */
    g_doublesum( &total_action );
    node0_printf("GACTION: %e\n",total_action/volume);
    /**node0_printf("CHECK:   %e   %e\n",total_action,imp_gauge_action() );**/

    if(this_node==0)fflush(stdout);
    special_free(tempmat1);
} /* g_measure_u1() */

void printpath( int *path, int length ){
    register int i;
    node0_printf("\t( ");
    for(i=0;i<length;i++)node0_printf("%d ",path[i]);
    node0_printf(",  L = %d )\n", length );
}

#ifdef N_SUBL32
/*** code from symanzik_sl32/dsdu_qhb.c  -- compute the staple ***/
/* This is a version for extended actions where 32 sublattices are
   needed to make the links independent. */
/* U.M. Heller August 1997 */

/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */
#include <assert.h>

/* NOTE: the staple is returned in the site strcture for use in monte
   and relax in symanzik_sl32! */

void dsdu_qhb_subl(int dir, int subl)
{
register site *st;
register int i;
int iloop, ln, k, j;
int dirs[MAX_LENGTH], length;
int path_dir[MAX_LENGTH], path_length;
su3_matrix tmat1, *tempmat1;
int fsubl;

 assert(NREPS==1);   /* This procedure designed only for NREPS = 1 */

 tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
 if(tempmat1 == NULL){
   printf("dsdu_qhb_subl: Can't malloc temporary\n");
   terminate(1);
 }

    FORSOMESUBLATTICE(i,st,subl) {
	clear_su3mat(&(st->staple));
    }

    for(iloop=0;iloop<NLOOP;iloop++){
	length=loop_length[iloop];
	for(ln=0;ln<loop_num[iloop];ln++){
	    /* set up dirs.  we are looking at loop starting in "XUP"
	       direction, rotate so it starts in "dir" direction. */
	    for(k=0;k<length;k++){
		if( GOES_FORWARDS(loop_table[iloop][ln][k]) ){
		    dirs[k]=(dir+loop_table[iloop][ln][k] )% 4;
		}
		else {
		    dirs[k]=OPP_DIR(
			(dir+OPP_DIR(loop_table[iloop][ln][k]))%4 );
		}
	    }

	    path_length = length-1;	/* generalized "staple" */
	    /* The path starts at the forward end of the link */
	    fsubl = neighsubl[subl][dir];

	    /* check for links in direction of link to be updated.
	       Note the direction of the path - opposite the link. */
	    for(k=0;k<length;k++)if( dirs[k]==dir||dirs[k]==OPP_DIR(dir)) {
		if( GOES_FORWARDS(dirs[k]) ) for(j=0;j<path_length;j++) {
		    path_dir[j] = dirs[(k+j+1)%length];
		}
		if( GOES_BACKWARDS(dirs[k]) ) for(j=0;j<path_length;j++) {
		    path_dir[path_length-1-j] =
			OPP_DIR(dirs[(k+j+1)%length]);
		}
		path_prod_subl(path_dir, path_length, fsubl, tempmat1);

		/* We took the path in the other direction from our old
		   convention in order to get it to end up "at our site".
		   So now take adjoint */
		FORSOMESUBLATTICE(i,st,subl) {
		    su3_adjoint( &tempmat1[i], &tmat1 );
		    scalar_mult_add_su3_matrix(&(st->staple), &tmat1,
			loop_coeff[iloop][0], &(st->staple) );
		}
	    } /* k (location in path) */
	} /* ln */
    } /* iloop */

    special_free(tempmat1);
    g_sync();

} /* dsdu_qhb */

#endif /* N_SUBL32 */
