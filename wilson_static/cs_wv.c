/**************************** cs_wv.c ****************************/
/* MIMD version 7 */
/*
*  $Header: /u/grad/winterow/cvs/graphene/wilson_static/cs_wv.c,v 1.1.1.1 2013/04/29 23:20:03 winterow Exp $
*
*  Overwrite the a wilson vector with the wilso vector
*  times a complex scalar
*
*  m --> phase * m
*
*/
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_wilsonvec(wilson_vector *m, complex *phase)
{
  register int ic, ispin;
  complex z ;

  for(ic=0;ic<3;ic++)
    for(ispin=0;ispin<4; ispin++)
    {
      z = cmul(&m->d[ispin].c[ic],phase) ;
      m->d[ispin].c[ic] =  z;
    }

}


