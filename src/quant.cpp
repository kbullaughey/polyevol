/*
 *  quant.cpp
 *
 *  Created by Richard Hudson on 8/1/11.
 *  Altered by Kevin Bullaughey on 10/26/11.
 *
 */

#include <math.h>
#include <valarray>
#include <vector>
#include <iostream>

#include "error_handling.h"
#include "command_line.h"
#include "sim_rand.h"
#include "genome.h"
#include "population.h"
#include "common.h"

#define OFFSPRING_POP (1-parent_pop)

using std::valarray;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::cin;

void usage(void);

int
main(int argc, char **argv) { try {
  if (argc == 1) {
    usage();
    return 0;
  }

  /* read in the command line arguments and print them out */
  Args ar(argc, argv);
  cout << ar << endl;

  /* set up simulation-wide genome parameters */
  int nloci = (int)ar.loci_counts.sum();
  Genome::initialize(ar.mu, 2.0/ar.s, ar.opts[0], ar.env);
  Population::initialize(ar.popsize, ar.sites_model);
  /* set the optimum to the first one */
  Genome::new_optimum(ar.opts[0]);

  /* bookkeeping for population simulation. I maintain two population objects, 
   * one for the parent generation and one for the offspring generation */
  Population *pops = new Population[2];

  /* model-specific setup */
  if (ar.sites_model == infinite_sites) {
    GenomeInfiniteSites::setup_effect_probabilities(ar.effect_probabilities, ar.effect_sizes);
  } else {
    /* loop over loci counts for each effect size and make a site with this effect */
    for (int i=0; i < (int)ar.loci_counts.size(); i++) {
      for (int j=0; j < ar.loci_counts[i]; j++) 
        Population::create_site(ar.effect_sizes[i]);
    }
  }

  /* initial frequencies, read in from standard input */
  if (nloci > 0) {
    valarray<int> heterozygotes(nloci);
    valarray<int> derived_homozygotes(nloci);
    double f;
    int loc = 0;
    while (1) {
      cin >> f;
      if (cin.eof() != 0 || loc == nloci) break;
      heterozygotes[loc] = (int)round(2.0*f*(1.0-f)*ar.popsize);
      derived_homozygotes[loc] = (int)round(f*f*ar.popsize);
      loc++;
    }
    if (!(cin.eof() && loc == nloci)) 
      throw SimError(0, "incorrect number of frequencies on stdin. I got more than %d frequencies.", nloci);
  
    /* initialize the first with the initial genotypes */
    pops[0].setup_initial_genotypes(heterozygotes, derived_homozygotes);
  }

  /* I use Dicks' trick of flipping back and forth between populations 
   * There's a macro OFFSPRING_POP which is defined as (1-parent_pop) to
   * make the code as readable as possible */
  int parent_pop = 0;

  /* epochs correspond to periods between which opt is constant and across which it changes */
  int gen = 0;
  for (int epoch=0; epoch < (int)ar.times.size(); epoch++) {
    /* update the optimum, for the first epoch, this has been done above */
    if (epoch > 0) Genome::new_optimum(ar.opts[epoch]);

    while (gen < ar.times[epoch]) {
      /* only print output if we've discarded the burnin */
      if (ar.burnin <= 0) {
        cout << "gen: " << gen << " "; pops[parent_pop].print_frequency_summary();
        cout << "gen: " << gen << " "; pops[parent_pop].print_phenotype_summary();
      } 

      /* advance the population simulation one generation */
      pops[OFFSPRING_POP].populate_from(pops[parent_pop]);

      /* if we're using the infinite sites model, we should do some cleanup */
      if (ar.sites_model == infinite_sites)
        pops[OFFSPRING_POP].purge_lost();
    
      /* generations start counting after the burnin is over */
      if (ar.burnin <= 0) {
        gen++;
      } else {
        ar.burnin--;
      }

      /* swap the parent and offspring in preparation for the next gen */
      parent_pop = OFFSPRING_POP;
    }
  } /* end of main loop */

  /* print the final state */
  cout << "gen: " << gen << " "; pops[parent_pop].print_frequency_summary();
  cout << "gen: " << gen << " "; pops[parent_pop].print_phenotype_summary();

/* catch any errors that were thrown anywhere inside this block */
} catch (SimUsageError e) {
   cerr << endl << "detected usage error: " << e.detail << endl << endl;
   usage();
   return 1;
} catch(SimError &e) {
   cerr << "uncaught exception: " << e.detail << endl;
   return 1;
} return 0; }

/* print a help message */
void
usage(void) {
  cerr << "usage: quant <options>\n"
    << "Required options:\n"
    << "  -m/--model infinite|finite\n"
    << "General options:\n"
    << "  -N/--popsize <int>    population size\n"
    << "  -u/--mu <float>       per locus mutation rate\n"
    << "  -s/--s <float>        Barton's selection parameter\n"
    << "  --seed=<int>          seed for random number generator\n"
    << "  --burnin=<int>        number of generations of burnin (discarded)\n"
    << "  --env=<float>         environmental variance\n"
    << "  --effects=<float vec> vector of effect sizes (comma-separated)\n"
    << "  --opts=<float vec>    vector of optima (comma-separated)\n"
    << "  --times=<int vec>     times (in generations) at which the optima change (comma-separated)\n"
    << "Infinite-sites-specific options:\n"
    << "  --eprobs=<double vec> effect size probabilities (comma-separated)\n"
    << "  --loci=<int>          number of loci to start the simulation (provided on stdin)\n"
    << "Finite-sites-specific options:\n"
    << "  --loci=<int vec>      number of loci of each effect size (comma-separated)\n"
    << "\n";
  return;
}

/* END */
