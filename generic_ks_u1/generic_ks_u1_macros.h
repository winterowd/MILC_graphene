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

#define FORALLMYUP_3_DIR(dir) FORALLUP_3_DIRBUT(ZUP, dir)

#define FORALLMYUPDIR(dir) FORALLUPDIRBUT(ZUP, dir)

#define FORMYSITESANDPARITY(i,s,choice) FORSOMEPARITYNOTZ(i,s,choice)

#define FORALLMYSITES(i,s) \
  for(i=0,s=lattice;i<sites_on_node;i++,s++) if(s->z==0)

#define FORMYODDSITES(i,s) FORODDSITES(i,s) if(s->z==0)

#define FORMYPARITYDOMAIN(i,s,parity) \
  FORMYSITESANDPARITY(i,s,choice) if(s->t > 0)
