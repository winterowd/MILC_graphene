/************************ grsource_imp_u1.c *****************************/
/* MIMD version 7 */
/* Added 7/30/2010, C.W.
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  
*/
#include "generic_ks_includes_u1.h"	/* definitions files and prototypes */
#include "../include/loopend.h"

/* construct a gaussian random vector, g_rand, and phi=M(dagger)*g_rand  */
/* "parity" is EVEN, ODD, or EVENANDODD.  The parity is the parity at
    which phi is computed.  g_rand must always be computed at all sites. */

void grsource_imp_u1( field_offset dest, Real mass, int parity,
		   ferm_links_u1_t *fn ) {
  register int i;
  register site *s;
  FORALLMYSITES(i,s){
#ifdef SITERAND
    s->g_rand.real = gaussian_rand_no(&(s->site_prn));
    s->g_rand.imag = gaussian_rand_no(&(s->site_prn));
#else
    s->g_rand.real = gaussian_rand_no(&node_prn);
    s->g_rand.imag = gaussian_rand_no(&node_prn);
#endif
  }
  
  dslash_fn_site_u1( F_OFFSET(g_rand), dest, parity, fn);
  scalar_mult_latvec_u1( dest, -1.0, dest, parity ); 
  scalar_mult_add_latvec_u1( dest, F_OFFSET(g_rand), 2.0*mass,
			    dest, parity );    
}/* grsource_imp_u1 */

/* construct a plain gaussian random vector in the site structure  */
/* "parity" is EVEN, ODD, or EVENANDODD. */

void grsource_plain_u1( field_offset dest, int parity ) {
  int i,j;
  site *s;
  complex *rand;
  FORMYSITESANDPARITY(i,s,parity){
#ifdef SITERAND
      rand = (complex *)F_PT(s,dest);
            rand->real = gaussian_rand_no(&(s->site_prn));
            rand->imag = gaussian_rand_no(&(s->site_prn));
#else
            rand->real = gaussian_rand_no(&node_prn);
            rand->imag = gaussian_rand_no(&node_prn);
#endif
  } END_LOOP    
}/* grsource */
