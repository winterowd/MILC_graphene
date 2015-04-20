/*********************** gauge_force_imp.c  -- ****************************/
/* MIMD version 7 */
/* gauge action stuff for improved action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file
* C.D. split from gauge_stuff.c 10/06 
* C.W. modified to U(1) 10/10 */

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

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

/* update the momenta with the gauge force */
void imp_gauge_force( Real eps, field_offset mom_off ){
    register int i,dir;
    register site *st;
    complex tmat1,tmat2;
    register Real eb;
    register complex* momentum;
    complex *staple, *tempmat1, temp_mul;
    
    /* lengths of various kinds of loops */
    int *loop_length = get_loop_length();
    /* number of rotations/reflections  for each kind */
    int *loop_num = get_loop_num();
    /* table of directions, 1 for each kind of loop */
    int ***loop_table = get_loop_table();
    /* table of coefficients in action, for various "representations"
	(actually, powers of the trace) */
    Real **loop_coeff = get_loop_coeff();
    int max_length = get_max_length();
    int nloop = get_nloop();
    int nreps = get_nreps();

#ifdef GFTIME
    int nflop = 153004;  /* For Symanzik1 action */
    double dtime;
#endif
    int j,k;
    int *dirs,length;
    int *path_dir,path_length;

    int ln,iloop;
    Real action,act2,new_term;

    int ncount;
    char myname[] = "imp_gauge_force";

#ifdef GFTIME
    dtime=-dclock();
#endif

    dirs = (int *)malloc(max_length*sizeof(int));
    if(dirs == NULL){
      printf("%s(%d): Can't malloc dirs\n",myname,this_node);
      terminate(1);
    }
    path_dir = (int *)malloc(max_length*sizeof(int));
    if(path_dir == NULL){
      printf("%s(%d): Can't malloc path_dir\n",myname,this_node);
      terminate(1);
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

    eb = eps*beta; //coefficient for U(1) action

    /* Loop over directions, update mom[dir] */
    for(dir=XUP; dir<=TUP; dir++){

      FORALLSITES(i,st) {
	staple[i]=cmplx(0.0,0.0);
      } END_LOOP

	ncount=0;
	for(iloop=0;iloop<nloop;iloop++){
	    length=loop_length[iloop];
	    for(ln=0;ln<loop_num[iloop];ln++){
/**printf("UPD:  "); printpath( loop_table[iloop][ln], length );**/
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

		path_length= length-1;  /* generalized "staple" */

		/* check for links in direction of momentum to be
		   updated, each such link gives a contribution. Note
		   the direction of the path - opposite the link. */
		for(k=0;k<length;k++)if( dirs[k]==dir||dirs[k]==OPP_DIR(dir)) {
		    if( GOES_FORWARDS(dirs[k]) ) for(j=0;j<path_length;j++) {
			path_dir[j] = dirs[(k+j+1)%length];
		    }
		    if( GOES_BACKWARDS(dirs[k]) ) for(j=0;j<path_length;j++) {
			path_dir[path_length-1-j] =
			    OPP_DIR(dirs[(k+j+1)%length]);
		    }
/**if(dir==XUP)printf("X_UPDATE PATH: "); printpath( path_dir, path_length );**/
		    path_product_u1(path_dir,path_length, tempmat1);

		    /* We took the path in the other direction from our
			old convention in order to get it to end up
			"at our site", so now take adjoint */
		    /* then compute "single_action" contribution to
			staple */
		    FORALLSITES(i,st){
		      CONJG( tempmat1[i], tmat1 );
		      new_term = loop_coeff[iloop][0];

			/* now we add in the higher representations */
			if(nreps > 1){
node0_printf("WARNING: THIS CODE IS NOT TESTED\n"); exit(0);
			    act2=1.0;
			    CMULJ_( st->link[dir], tmat1, temp_mul );
			    action = 1.0 - temp_mul.real;
			    /* action = 3.0 - realtrace_su3(&(st->link[dir]),
			       &tmat1 ); */

			    for(j=1;j<nreps;j++){
				act2 *= action;
				new_term +=
				    loop_coeff[iloop][j]*act2*(Real)(j+1);
			    }
			}  /* end if nreps > 1 */

			CMULREAL( tmat1, new_term, temp_mul );
			CADD( temp_mul, staple[i], staple[i] );
			/* scalar_mult_add_su3_matrix( &(staple[i]), &tmat1,
			   new_term, &(staple[i]) ); */

		    } END_LOOP

		    ncount++;

		} /* k (location in path) */
	    } /* ln */
	} /* iloop */

	/* Now multiply the staple sum by the link, then update momentum */
	FORALLSITES(i,st){
	  CMUL_J( st->link[dir], staple[i], tmat1 );
	  //mult_su3_na( &(st->link[dir]), &(staple[i]), &tmat1 );
	  momentum = (complex *)F_PT(st,mom_off);
	  tmat2.real = 0.0;
	  tmat2.imag = momentum[dir].imag;
	  //uncompress_anti_hermitian( &momentum[dir], &tmat2 );
	  CMULREAL( tmat1, eb, temp_mul );
	  CSUB( tmat2, temp_mul, staple[i] );
	  /*scalar_mult_sub_su3_matrix( &tmat2, &tmat1,
	    eb, &(staple[i]) ); */
	  momentum[dir].real = 0.0;
	  momentum[dir].imag = staple[i].imag;
	  //make_anti_hermitian( &(staple[i]), &momentum[dir] );
	} END_LOOP
    } /* dir loop */
#ifdef GFTIME
dtime+=dclock();
node0_printf("GFTIME:   time = %e (Symanzik1) mflops = %e\n",dtime,
	     nflop*(double)volume/(1e6*dtime*numnodes()) );
#endif
 free(path_dir);
 free(dirs);
 special_free(staple); 
 special_free(tempmat1); 
} /* imp_gauge_force.c */

/* update the momenta with the non-compact gauge force */
void gauge_force( Real eps ) {
register int i,dir1,dir2;
register site *st;
msg_tag *tag0,*tag1,*tag2;
int start;
Real t1,t2;
register Real eb;
/**TEMP**
Real gf_x,gf_av,gf_max;
int gf_i,gf_j;
gf_av=gf_max=0.0;
**END TEMP**/

#ifdef GFTIME
 double dtime;
 dtime = -dclock();
#endif
 
    eb = eps*beta;
    /* Loop over directions, update mom[dir1] */
    //#ifdef FOUR_DIM
    for(dir1=XUP; dir1<=TUP; dir1++){
      /*#else
    for(dir1=TUP; dir1<=TUP; dir1++){
    #endif*/
	/* Loop over other directions, computing force from plaquettes in
	   the dir1,dir2 plane */
	start=1; /* indicates staple sum not initialized */
	for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1){

	    /* get potential[dir2] from direction dir1 */
	    tag0 = start_gather_site( F_OFFSET(potential[dir2]), 
	       sizeof(Real),dir1, EVENANDODD, gen_pt[0] );

	    /* Start gather for the "upper staple" */
	    tag2 = start_gather_site( F_OFFSET(potential[dir1]), sizeof(Real),
		dir2, EVENANDODD, gen_pt[2] );

	    /* begin the computation "at the dir2DOWN point", we will
		later gather the intermediate result "to the home point" */

	    wait_gather(tag0);
	    FORALLSITES(i,st){
	      t1 = st->potential[dir1] - st->potential[dir2];
	      st->temp1 = t1 + (*(Real *)gen_pt[0][i]);
	      /*mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1 );
	        mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[0][i],
		  &(st->tempmat1) ); */
	    }

	    /* Gather this partial result "up to home site" */
	    tag1 = start_gather_site( F_OFFSET(temp1), sizeof(Real),
		OPP_DIR(dir2), EVENANDODD, gen_pt[1] );

	    /* begin the computation of the "upper" staple.  Note that
		one of the links has already been gathered, since it
		was used in computing the "lower" staple of the site
		above us (in dir2) */
	    wait_gather(tag2);
	    if(start){	/* this is the first contribution to staple */
	        FORALLSITES(i,st){
		  t1 = st->potential[dir2] + (*(Real *)gen_pt[2][i]);
		  st->stple = -t1 + (*(Real *)gen_pt[0][i]);
		  //stple collects the "staple" part of phi_{\mu, \nu}
		  /*mult_su3_nn( &(st->link[dir2]), 
			(su3_matrix *)gen_pt[2][i], &tmat1);
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
		    &(st->staple) ); */
		}
		start=0;
	    }
	    else{
	        FORALLSITES(i,st){
		  t1 = st->potential[dir2] + (*(Real *)gen_pt[2][i]);
		  t2 = -t1 + (*(Real *)gen_pt[0][i]);
		  st->stple = t2 + st->stple;
		  /*mult_su3_nn( &(st->link[dir2]), 
			(su3_matrix *)gen_pt[2][i], &tmat1);
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
		    add_su3_matrix( &(st->staple),&tmat2,&(st->staple));*/
	        }
	    }

	    wait_gather(tag1);
	    /* subtract lower "staple" from the upper "staple */
	    FORALLSITES(i,st){
	      st->stple = st->stple - (*(Real *)gen_pt[1][i]);
	      /*add_su3_matrix( &(st->staple), (su3_matrix *)gen_pt[1][i],
		&(st->staple));*/
	    }
	    cleanup_gather(tag0);
	    cleanup_gather(tag1);
	    cleanup_gather(tag2);
	}
	/* Now add 2*potential[dir1] to "staple" and update 
	   (purely complex) momenta */
	FORALLSITES(i,st){
	  t1 =  6*st->potential[dir1] + st->stple;
	  st->mom[dir1].real = 0.0;
	  st->mom[dir1].imag = st->mom[dir1].imag - eb*t1;
	  /*mult_su3_na( &(st->link[dir1]), &(st->staple), &tmat1 );
	    uncompress_anti_hermitian( &(st->mom[dir1]), &tmat2 );
	    scalar_mult_add_su3_matrix( &tmat2, &tmat1,
		eb3, &(st->staple) );
		make_anti_hermitian( &(st->staple), &(st->mom[dir1]) );*/
/** FIND AVERAGE AND MAXIMUM GAUGE FORCE **/
/**TEMP **
for(gf_x=0,gf_i=0;gf_i<3;gf_i++)for(gf_j=0;gf_j<3;gf_j++)
gf_x+=cabs_sq(&(tmat1.e[gf_i][gf_j]));
gf_x *= beta/3.0;
gf_av += gf_x;
if(gf_x > gf_max) gf_max = gf_x;
** END TEMP **/
	}
    }
/**TEMP**
g_floatsum( &gf_av );
g_floatmax( &gf_max );
gf_av /= (4*volume);
if(this_node==0)printf("GF: %e %e\n",gf_av,gf_max);
**END TEMP **/
#ifdef GFTIME
    dtime += dclock();
    if(this_node==0)printf("G_FORCE: time = %e mflops = %e\n",
			   dtime, (double)(10848.0*volume/(1.0e6*dtime*numnodes())) );
#endif

}// gauge_force()
