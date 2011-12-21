/*
 *  qapprox2.c
 *  
 *
 *  Created by Richard Hudson on 8/1/11.
 *
 */
/* gcc -o qapprox2 qapprox2.c  -lm -lgsl */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#   alphas come from the file:   "alphas",    initial frequencies from stdin  


main(int argc, char **argv )
{

int  **genos[2] ;
FILE *fp1, *palpha  ;
int popsize, nloci, ngen , gen  ;
int paroff, matpat, ind, loc , nfix ;
double ran1() ;
double s,  alpha, sum , sumsq,  env, initfreq1, *freqinit  ;
double *effects,  opt, sig , fitmax, mutrate, *freqsp, *freqsoff , *ptemp ;
double var, meanphen, sumpa, suma, delta, pprime ;

if( argc < 7 ) { fprintf( stderr, " qapprox2  popsize nloci optimum  s(strength of selection) ngen  mutrate \n"); exit(1); }
popsize = strtol( argv[1],NULL,10);
nloci = strtol( argv[2],NULL,10);
opt = strtod( argv[3],NULL);
s = strtod( argv[4],NULL);
ngen = strtol( argv[5],NULL,10);
	mutrate = strtod( argv[6],NULL);
/*	alpha = strtod( argv[7],NULL); */
	/*	nfix = strtol( argv[8],NULL,10); */
/*		initfreq1  = strtod( argv[9],NULL); */

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng, 123 );

/* to gen a binomial: 		nwt = gsl_ran_binomial( rng, epwt,(unsigned int)n) ;   */

	env = 0. ;

fp1 = fopen("phenout", "w");
	
/* vectors to store  effect size at each locus,  phenotypes of individuals, and fitnesses of individuals */
	effects = (double *)malloc( (unsigned)nloci*sizeof(double) ) ;
	freqsp = (double *)malloc( (unsigned)nloci*sizeof(double) ) ;
	freqsoff = (double *)malloc( (unsigned)nloci*sizeof(double) ) ;
	freqinit = (double *)malloc( (unsigned)nloci*sizeof(double) ) ;

palpha = fopen( "alphas","r");
for( loc = 0 ; loc<nloci; loc++)   fscanf(palpha," %lf", effects+loc)   ;
for( loc = 0; loc <nloci ; loc++) fscanf(stdin," %lf", freqinit+loc ) ;



for( loc = 0; loc < nloci; loc++) printf(" %lf", freqinit[loc] );
printf("\n");
for( loc = 0; loc < nloci; loc++) freqsp[loc] = freqinit[loc] ;
	
suma = 0.0 ;
for( loc=0; loc <nloci; loc++) suma += effects[loc] ;

/* main loop */
for( gen = 0 ; gen < ngen ; gen++){
/*   if( gen == 3000 ) opt += 2.0 ; */ 

	sumpa = var =  0. ;
	for( loc=0; loc<nloci; loc++) { 
	       sumpa +=  freqsp[loc]*effects[loc] ; 
		   var += 2.*effects[loc]*effects[loc]*freqsp[loc]*(1.0-freqsp[loc] ) ;
		   } 
	meanphen = 2.*sumpa - suma ;
	 delta = meanphen - opt ;
	 fprintf(fp1,"%lf\t%lf\n",meanphen, var ) ;

 /*   pprime =    */
    for( loc = 0; loc< nloci; loc++){
		pprime = freqsp[loc]  + 0.5*s*effects[loc]*effects[loc]*freqsp[loc]*(1.0-freqsp[loc] )*( (2.0*freqsp[loc]-1.0) - 2.0*delta/effects[loc] ) - mutrate*(2.0*freqsp[loc] -1.0 ) ;
	     freqsoff[loc] = gsl_ran_binomial( rng, pprime, (unsigned int)(2*popsize) )/(2.0*popsize) ;
	}
	for( loc = 0; loc < nloci; loc++) printf(" %lf", freqsoff[loc] );
	printf("\n");
	
/*		fprintf(fp1,"%lf\t%lf\n",sum/popsize, sumsq/popsize - (sum/popsize)*(sum/popsize) ) ; */
	ptemp = freqsp;
	freqsp = freqsoff ;
	freqsoff = ptemp ;
  }
/* end of main loop */
}

