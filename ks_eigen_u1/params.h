#ifndef _PARAMS_H
#define _PARAMS_H

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;   /* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
	int iseed;	/* for random numbers */
    /*  REPEATING BLOCK */
	Real mass; /* gauge coupling, quark mass */
        Real u0_s; /* space-like tadpole parameter */
	Real u0_t; /* time-like tadpole parameter */
        Real v_Fermi; /* Fermi velocity */
	int niter; 	/* maximum number of c.g. iterations */
	int nrestart; 	/* maximum number of c.g. restarts */
	Real rsqmin,rsqprop;  /* for deciding on convergence */
        int Nvecs ; /* number of eigenvectors */
        Real eigenval_tol ; /* Tolerance for the eigenvalue computation */
        Real error_decr ; /* error decrease per Rayleigh minimization */
        int MaxIter ; /* max  Rayleigh iterations */
        int Restart ; /* Restart  Rayleigh every so many iterations */
        int Kiters ; /* Kalkreuter iterations */
	int startflag;  /* what to do for beginning lattice */
	char startfile[MAXFILENAME];
}  params;

#endif /* _PARAMS_H */
