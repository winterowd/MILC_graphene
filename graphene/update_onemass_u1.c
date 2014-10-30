/********** update.c ****************************************************/
/* MIMD version 7 */

/*
 Update lattice.
 Improved method for 1-4 flavors:
	update U by (epsilon/2)*(1-Nf/4)
	compute PHI
	update U to epsilon/2
	compute X
	update H, full step
	update U to next time needed

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
*/
#include "graphene_includes.h"	/* definitions files and prototypes */

int update_u1()  {
int step, iters=0;
Real final_rsq;
void predict_next_xxx(Real *oldtime,Real *newtime,Real *nexttime);
Real cg_time;	/* simulation time for last two CG's */
double temp_ss, temp_st;
#ifdef PHI_ALGORITHM
Real old_cg_time,next_cg_time;	/* simulation time for last two CG's */
#endif
#ifdef HMC_ALGORITHM
double startaction,endaction,d_action_u1();
Real xrandom;
#endif
//rephase(OFF);
 //gauge_force(epsilon); 
//rephase(ON);
 /* refresh the momenta */
 ranmom_u1();
 /* register int i, dir;
 register site *s;
 FORALLSITES(i,s) {
   node0_printf("(%d, %d, %d, %d): \n",s->x,s->y,s->z,s->t);
   for(dir=XUP;dir<=TUP;dir++)
     node0_printf("mom[%d] = %e ,", dir, s->mom[dir].imag);
   node0_printf("\n");		     
   }*/

/*DEBUG*/
/**grsource_imp(F_OFFSET(phi), mass, EVENANDODD);
ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,niter,nrestart,rsqmin,PRECISION,EVENANDODD,&final_rsq);
checkmul();**/
/*ENDDEBUG*/

    /* do "steps" microcanonical steps  */
    for(step=1; step <= steps; step++){
 
#ifdef PHI_ALGORITHM
        /* generate a pseudofermion configuration only at start*/
      if(step==1){
	load_ferm_links_u1(&fn_links, &ks_act_paths);
	/*f_meas_imp_u1(F_OFFSET(phi),F_OFFSET(xxx),mass, &fn_links, 
	  &fn_links_dmdu0);*/
	//d_plaquette_u1(&temp_ss, &temp_st);
	clear_latvec_u1(F_OFFSET(xxx), EVENANDODD);
	grsource_imp_u1(F_OFFSET(phi), mass, EVEN, &fn_links); 
	old_cg_time = cg_time = -1.0e6;
      }

#ifdef HMC_ALGORITHM
        /* find action */
        /* do conjugate gradient to get (Madj M)inverse * phi */
        if(step==1){
	  /* do conjugate gradient to get (Madj M)inverse * phi */
	  load_ferm_links_u1(&fn_links, &ks_act_paths);
	  iters += ks_congrad_u1( F_OFFSET(phi), F_OFFSET(xxx), mass, 
				  niter, nrestart, rsqmin, PRECISION, 
				  EVEN, &final_rsq, &fn_links );
	  cg_time = 0.0;

     	    startaction=d_action_u1();
            /* copy link field to old_link */
	    gauge_field_copy_u1( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
#ifdef NON_COMPACT
	    gauge_field_copy_nc( F_OFFSET(potential[0]), F_OFFSET(old_potential[0]));
#endif
        }
#endif

	/* update U's to middle of interval */
	test_flag = 1;
     	update_u_u1(0.5*epsilon);

	/*node0_printf("AFTER FIRST UPDATE_U!!!\n");
	FORALLSITES(i,s) {
	node0_printf("(%d, %d, %d, %d): \n",s->x,s->y,s->z,s->t);
	for(dir=XUP;dir<=TUP;dir++) {
	node0_printf("link[%d] = %e %e, potential[%d] = %e",dir, 
	s->link[dir].real, s->link[dir].imag, dir, s->potential[dir]);
	}
	node0_printf("\n");		     
	}*/
	
	/* save conjugate gradient solution, predict next one */
	next_cg_time = ((Real)step-0.5)*epsilon;
	predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);

#else /* "R" algorithm */
       	/* first update the U's to special time interval */
       	update_u_u1(epsilon*(0.5-nflavors/8.0));

        /* generate a pseudofermion configuration */
	load_ferm_links_u1(&fn_links, &ks_act_paths);
     	grsource_imp_u1(F_OFFSET(phi), mass, EVEN, &fn_links); 
	cg_time = -1.0e6;

	/* update U's to middle of interval */
     	update_u_u1(epsilon*nflavors/8.0);
#endif

        /* do conjugate gradient to get (Madj M)inverse * phi */
	load_ferm_links_u1(&fn_links, &ks_act_paths);
     	iters += ks_congrad_u1( F_OFFSET(phi), F_OFFSET(xxx), mass, 
			     niter, nrestart, rsqmin, PRECISION, 
			     EVEN, &final_rsq, &fn_links );
	dslash_fn_site_u1( F_OFFSET(xxx), F_OFFSET(xxx), ODD, &fn_links);
	cg_time = ((Real)step - 0.5)*epsilon;
	/* now update H by full time interval */
    	update_h_u1(epsilon);
	
	/*node0_printf("AFTER UPDATE_H!!!\n");
	FORALLSITES(i,s) {
	node0_printf("(%d, %d, %d, %d): \n",s->x,s->y,s->z,s->t);
	for(dir=XUP;dir<=TUP;dir++) {
	node0_printf("mom[%d] = %e, ", dir, s->mom[dir].imag);
	}
	node0_printf("\n");		     
	}*/

    	/* update U's by half time step to get to even time */
    	update_u_u1(epsilon*0.5); 

	/*node0_printf("AFTER SECOND UPDATE_U!!!\n");
	FORALLSITES(i,s) {
	node0_printf("(%d, %d, %d, %d): \n",s->x,s->y,s->z,s->t);
	for(dir=XUP;dir<=TUP;dir++) {
	node0_printf("link[%d] = %e %e, potential[%d] = %e",dir, s->link[dir].real,
	s->link[dir].imag, dir, s->potential[dir]);
	}
	node0_printf("\n");		     
	}*/

        /* reunitarize the gauge field */
	rephase( OFF );
	//node0_printf("before reunitarize step %d\n", step);
        reunitarize_u1();
	rephase( ON );

	//d_plaquette_u1(&temp_ss, &temp_st);
	//printf("temp_ss = %e, temp_st = %e \n", temp_ss, temp_st);
    }	/* end loop over microcanonical steps */

#ifdef HMC_ALGORITHM
    /* find action */
    /* do conjugate gradient to get (Madj M)inverse * phi */
    next_cg_time = steps*epsilon;
    predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);
    load_ferm_links_u1(&fn_links, &ks_act_paths);
    iters += ks_congrad_u1( F_OFFSET(phi), F_OFFSET(xxx), mass,
			 niter, nrestart, rsqmin, PRECISION, 
			 EVEN, &final_rsq, &fn_links );
    cg_time = steps*epsilon;
    endaction=d_action_u1();
    /* decide whether to accept, if not, copy old link field back */
    /* careful - must generate only one random number for whole lattice */
    if(this_node==0)xrandom = myrand(&node_prn);
    broadcast_float(&xrandom);
    if( exp( (double)(startaction-endaction) ) < xrandom ){
      if(steps > 0)  {
	gauge_field_copy_u1( F_OFFSET(old_link[0]), F_OFFSET(link[0]) );
#ifdef NON_COMPACT
	gauge_field_copy_nc( F_OFFSET(old_potential[0]), F_OFFSET(potential[0]));
#endif
      }
#ifdef FN
	free_fn_links_u1(&fn_links);
	free_fn_links_u1(&fn_links_dmdu0);
#endif
	node0_printf("REJECT: delta S = %e\n", (double)(endaction-startaction));
	num_rejects++;
    }
    else {
	node0_printf("ACCEPT: delta S = %e\n", (double)(endaction-startaction));
	num_accepts++;
    }
#endif

    if(steps > 0)return (iters/steps);
    else return(-99);
}

#ifdef PHI_ALGORITHM
/* use linear extrapolation to predict next conjugate gradient solution */
/* only need even sites */
void predict_next_xxx(Real *oldtime,Real *newtime,Real *nexttime)
{
register int i;
register site *s;
register Real x;
 complex tvec, temp_mul;
    if( *newtime != *oldtime ) x = (*nexttime-*newtime)/(*newtime-*oldtime);
    else x = 0.0;
    if( *oldtime < 0.0 ){
        FORMYEVENSITES(i,s){
	    s->old_xxx = s->xxx;
        }
    }
    else  {
        FORMYEVENSITES(i,s){
            CSUB( s->xxx, s->old_xxx, tvec);
	    s->old_xxx = s->xxx;
	    CMULREAL( tvec, x, temp_mul );
	    CADD( temp_mul, s->xxx, s->xxx );
	}
    }
    *oldtime = *newtime;
    *newtime = *nexttime;
}
#endif
