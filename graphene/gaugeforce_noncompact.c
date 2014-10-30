/****** gaugeforce_noncompact.c  ******************/
/* MIMD version 6 */

#include "generic_includes.h"

/* update the momenta with the non-compact gauge force */
void gauge_force( Real eps ) {
register int i,dir1,dir2;
register site *st;
msg_tag *tag0,*tag1,*tag2;
int start;
Real t1,t2;
register Real eb2;
/**TEMP**
Real gf_x,gf_av,gf_max;
int gf_i,gf_j;
gf_av=gf_max=0.0;
**END TEMP**/

/**double dtime,dclock();
dtime = -dclock();**/

    eb3 = eps*beta*2.0;
    /* Loop over directions, update mom[dir1] */
    for(dir1=XUP; dir1<=TUP; dir1++){
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
		  /*mult_su3_nn( &(st->link[dir2]), 
			(su3_matrix *)gen_pt[2][i], &tmat1);
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
		    &(st->staple) ); */
		}
		start=0;
	    }
	    else{
	        FORALLSITES(i,st){
		    mult_su3_nn( &(st->link[dir2]), 
			(su3_matrix *)gen_pt[2][i], &tmat1);
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
		    add_su3_matrix( &(st->staple),&tmat2,&(st->staple));
	        }
	    }

	    wait_gather(tag1);
	    FORALLSITES(i,st){
		add_su3_matrix( &(st->staple), (su3_matrix *)gen_pt[1][i],
		    &(st->staple));
	    }
	    cleanup_gather(tag0);
	    cleanup_gather(tag1);
	    cleanup_gather(tag2);
	}
	/* Now multiply the staple sum by the link, then update momentum */
	FORALLSITES(i,st){
	    mult_su3_na( &(st->link[dir1]), &(st->staple), &tmat1 );
	    uncompress_anti_hermitian( &(st->mom[dir1]), &tmat2 );
	    scalar_mult_add_su3_matrix( &tmat2, &tmat1,
		eb3, &(st->staple) );
	    make_anti_hermitian( &(st->staple), &(st->mom[dir1]) );
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
/**dtime += dclock();
if(this_node==0)printf("G_FORCE: time = %e mflops = %e\n",
dtime, (double)(10848.0*volume/(1.0e6*dtime*numnodes())) );**/
}// gauge_force()
