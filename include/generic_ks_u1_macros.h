/* generic_ks_u1_macros.h */
/* Some MACROS for the 2+1 U(1) Gauge Theory */
/* Adapted from macros.h */

#define FORSOMEPARITYNOTZ(i,s,choice) \
{ register int loopend;  \
loopend= (choice)==EVEN ? even_sites_on_node : sites_on_node ; \
for( i=((choice)==ODD ? even_sites_on_node : 0 ), s= &(lattice[i]); \
i<loopend; i++,s++)if(s->z==0)
#define END_LOOP }

#define FORALLUP_3_DIR(dir) for(dir=X3UP; dir<=T3UP; dir++)

#define FORALLUP_3_DIRBUT(direction,dir) \
   FORALLUP_3_DIR(dir)if(dir != direction)

#ifndef FOUR_DIM
#define FORALLMYUP_3_DIR(dir) FORALLUP_3_DIRBUT(Z3UP, dir)

#define FORALLMYUPDIR(dir) FORALLUPDIRBUT(ZUP, dir)

#define FORMYSITESANDPARITY(i,s,choice) FORSOMEPARITYNOTZ(i,s,choice)

#define FORALLMYSITES(i,s) \
  for(i=0,s=lattice;i<sites_on_node;i++,s++) if(s->z==0)

#define FORMYODDSITES(i,s) FORODDSITES(i,s) if(s->z==0)

#define FORMYEVENSITES(i,s) FOREVENSITES(i,s) if(s->z==0)
#else
#define FORMYSITESANDPARITY(i,s,choice) FORSOMEPARITY(i,s,choice)
#define FORALLMYUP_3_DIR(dir) FORALLUP_3_DIR(dir)
#define FORALLMYUPDIR(dir) FORALLUPDIR(dir)
#define FORMYEVENSITES(i,s) FOREVENSITES(i,s)
#define FORMYODDSITES(i,s) FORODDSITES(i,s)
#define FORALLMYSITES(i,s) FORALLSITES(i,s)
#endif

#ifdef SCHROED_FUN
#define FORMYEVENSITESDOMAIN(i,s) \
  FORMYEVENSITES(i,s) if(s->t > 0)
#define FORMYODDSITESDOMAIN(i,s) \
  FORMYODDSITES(i,s) if(s->t > 0)
#define FORALLMYSITESDOMAIN(i,s) \
  FORALLMYSITES(i,s) if(s->t > 0)
#define FORMYPARITYDOMAIN(i,s,parity) \
  FORMYSITESANDPARITY(i,s,parity) if(s->t > 0)
#else
#define FORMYPARITYDOMAIN FORMYSITESANDPARITY
#define FORMYEVENSITESDOMAIN FORMYEVENSITES
#define FORMYODDSITESDOMAIN FORMYODDSITES
#define FORALLMYSITESDOMAIN FORALLMYSITES
#endif
