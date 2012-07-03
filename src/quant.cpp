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
#include "statistic.h"

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

  /* set up the default statistics */
  Statistic::initialize_defaults(); 

  /* read in the command line arguments and print them out */
  Args ar(argc, argv);
  cout << ar << endl;

  /* set up simulation-wide genome parameters */
  Genome::initialize(ar.mu, 2.0/ar.s, ar.opts[0], ar.env);
  Site::ploidy_level = ar.ploidy_level;
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
      for (int j=0; j < ar.loci_counts[i]; j++) {
        Population::create_site(ar.effect_sizes[i]);
        /* unlike the Barton model, I represent genotypes 0,1,2 instead of 
         * -1,0,1. This means I need to store the baseline of the genotype 
         * that is entirely homozygote ancestral, and add the effects of 
         * derived alleles on top of that */
        GenomeFiniteSites::baseline -= ar.effect_sizes[i];
      }
    }
  }

  /* initial frequencies, read in from standard input */
  if (ar.nloci > 0) {
    if (ar.ploidy_level == haploid) 
      throw SimError("haploid version doesn't support frequencies on stdin");
    valarray<int> heterozygotes(ar.nloci);
    valarray<int> derived_homozygotes(ar.nloci);
    double f;
    int loc = 0;
    while (1) {
      /* this function works like a generator, so we can call it repeatedly. 
       * It return -1 when it's done generating */
      f = ar.get_initial_frequency();
      if (f < 0 || loc == ar.nloci) break;
      heterozygotes[loc] = (int)round(2.0*f*(1.0-f)*ar.popsize);
      derived_homozygotes[loc] = (int)round(f*f*ar.popsize);
      loc++;
    }
    if (!(f < 0 && loc == ar.nloci)) 
      throw SimError(0, "incorrect number of frequencies. Expecting %d.", ar.nloci);
  
    /* initialize the first with the initial genotypes */
    pops[0].setup_initial_genotypes(heterozygotes, derived_homozygotes);
  }

  /* I use Dicks' trick of flipping back and forth between populations 
   * There's a macro OFFSPRING_POP which is defined as (1-parent_pop) to
   * make the code as readable as possible */
  int parent_pop = 0;

  /* epochs correspond to periods between which opt is constant and across which it changes */
  Population::generation = 0;
  for (int epoch=0; epoch < (int)ar.times.size(); epoch++) {
    /* update the optimum, for the first epoch, this has been done above */
    if (epoch > 0) Genome::new_optimum(ar.opts[epoch]);

    while (Population::generation < ar.times[epoch]) {
      /* only print output if we've discarded the burnin */
      if (ar.burnin <= 0) {
        pops[parent_pop].stat_frequency_summary();
        pops[parent_pop].stat_phenotype_summary();
        pops[parent_pop].stat_increment_visits();
        pops[parent_pop].stat_fixations();
        pops[parent_pop].stat_segsites();
      } 

      /* advance the population simulation one generation */
      pops[OFFSPRING_POP].populate_from(pops[parent_pop]);

      /* if we're using the infinite sites model, we should do some cleanup */
      if (ar.sites_model == infinite_sites)
        pops[OFFSPRING_POP].purge_lost();
    
      /* generations start counting after the burnin is over */
      if (ar.burnin == 0 && Population::generation == 0 && Statistic::is_activated("burnin"))
        cout << "end burnin" << endl;
      if (ar.burnin <= 0) {
        Population::generation++;
        if (Population::generation == 0) {
          if (Statistic::is_activated("burnin"))
            cout << "burnin mutations: " << Genome::mutation_count << endl;
          Genome::mutation_count = 0;
        }
        pops[OFFSPRING_POP].stat_update_p_moments();
      } else {
        ar.burnin--;
      }

      /* swap the parent and offspring in preparation for the next gen */
      parent_pop = OFFSPRING_POP;
    }
  } /* end of main loop */

  /* print the final state */
  pops[parent_pop].stat_frequency_summary();
  pops[parent_pop].stat_phenotype_summary();
  pops[parent_pop].stat_segsites();
  Population::stat_print_visits();
  Population::stat_print_p_moments();
  if (Statistic::is_activated("mutation")) 
    cout << "mutations: " << Genome::mutation_count << endl;

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
    << "  --effects=<float vec> vector of effect sizes (comma-separated)\n"
    << "  --times=<int vec>     times (in generations) at which the optima change (comma-separated)\n"
    << "  --opts=<float vec>    vector of optima (comma-separated)\n"
    << "  --loci=<int>          number of loci to start the simulation (provided on stdin)\n"
    << "Optional options:\n"
    << "  -N/--popsize <int>    population size\n"
    << "  -u/--mu <float>       per locus mutation rate\n"
    << "  -s/--s <float>        Barton's selection parameter\n"
    << "  --seed=<int>          seed for random number generator\n"
    << "  --burnin=<int>        number of generations of burnin discarded\n"
    << "  --env=<float>         environmental variance\n"
    << "  --haploid             use a haploid population (default is diploid)\n"
    << "Infinite-sites-specific options:\n"
    << "  --eprobs=<double vec> effect size probabilities (comma-separated)\n"
    << "Finite-sites-specific options:\n"
    << "  --loci=<int vec>      number of loci of each effect size (comma-separated)\n"
    << "Initial frequency initialization:\n"
    << "  --freqs stdin | even   method for initializing frequencies\n"
    << "      stdin: read in frequencies from standard input\n"
    << "      even: evenly spaced. ith freq is i/(1+n) where n is the number of loci\n"
    << "Statistics control:\n"
    << "  --enable-stat=<str>   enable a statistic\n"
    << "  --disable-stat=<str>  disable a statistic\n"
    << "  --disable-all-stats   turn off all statistics (must precede enable options)\n"
    << "      Available statistics (default):\n"
    << "        frequencies     print allele IDs and frequencies (on)\n"
    << "        phenotype       mean phenotype and variance (on)\n"
    << "        mutation        ID and generation for each new mutation (on)\n"
    << "        sojourn         sojourn time in gen for infinite sites model (on)\n"
    << "        burnin          notices about the burnin period (on)\n"
    << "        visits          report (final) number of visits to each allele count (off)\n"
    << "        fixations       number of fixations of each effect size (off)\n"
    << "        segsites        number of segregating sites of each effect size (off)\n"
    << "        p_moments       empirical first and second moments of the change in allele frequency (off)\n"
    << "\n";
  return;
}

/* END */
