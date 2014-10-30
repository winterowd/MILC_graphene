/****************** ploop2_u1.c ************************************/
/* MIMD version 7 */
/* evaluate the Polyakov loops for U(1).  This version uses gathers. */
/* Updated 10/10, C.W. */

#include "generic_includes.h"

complex ploop_u1() {
  register int i,t;
  register site *st;
  msg_tag *tag;
  complex sum;
  complex plp;
  complex *tempmat1, *tempmat2, *staple;

  staple = (complex *)malloc(sites_on_node*sizeof(complex));
  if(staple == NULL){
    printf("ploop: Can't allocate temporary\n");
    terminate(1);
  }

  tempmat1 = (complex *)malloc(sites_on_node*sizeof(complex));
  if(tempmat1 == NULL){
    printf("ploop: Can't allocate temporary\n");
    terminate(1);
  }

  tempmat2 = (complex *)malloc(sites_on_node*sizeof(complex));
  if(tempmat2 == NULL){
    printf("ploop: Can't allocate temporary\n");
    terminate(1);
  }

  sum = cmplx(0.0,0.0);
  FORALLSITES(i,st){tempmat2[i] = lattice[i].link[TUP];}
  for(t=1;t<nt;t++){
    tag=start_gather_field( tempmat2, sizeof(complex),
			    TUP, EVENANDODD, gen_pt[0] );
    wait_gather(tag);
    FORALLSITES(i,st){
      if( st->t != 0 )continue;
      if(t==1){
	CMUL( st->link[TUP], *(complex *)gen_pt[0][i], staple[i] );
	/*mult_su3_nn( &(st->link[TUP]), (complex *)gen_pt[0][i],
	  &staple[i]); */
      }
      else {
	CMUL( staple[i], *(complex *)gen_pt[0][i], tempmat2[i] );
	/*mult_su3_nn( &staple[i], (complex *)gen_pt[0][i],
	  &(tempmat2[i])); */
	staple[i] = tempmat2[i];
      }
    }
    FORALLSITES(i,st){
      tempmat1[i] = *(complex *)(gen_pt[0][i]);
    }
    FORALLSITES(i,st){
     tempmat2[i] = tempmat1[i];
    }
    cleanup_gather(tag);
  }
  FORALLSITES(i,st){
    if( st->t != 0 )continue;
    plp.real = staple[i].real; plp.imag = staple[i].imag;
    CSUM(sum,plp);
  }
  g_complexsum( &sum );
  plp.real = sum.real /((Real)(nx*ny*nz));
  plp.imag = sum.imag /((Real)(nx*ny*nz));
  free(tempmat1);
  free(tempmat2);
  free(staple);
  return(plp);
}
