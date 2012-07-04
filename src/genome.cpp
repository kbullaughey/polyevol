#include "error_handling.h"
#include "common.h"
#include "genome.h"
#include "sim_rand.h"
#include "site.h"
#include "population.h"

using std::vector;
using std::valarray;
using std::ostream;
using std::cerr;
using std::endl;

/************************* 
 * Genome abstract class *
 *************************/

/* static member variables */
double Genome::mu;
double Genome::environmental_noise;
double Genome::optimum;
double Genome::sig;
double Genome::baseline;
int Genome::mutation_count;

/* A new genome is created with no mutant alleles */
Genome::Genome(Population *p, int indiv) { 
  pop = p;
  individual = indiv;
  update_phenotype();
  update_fitness();
}

/* used to set a few class variables */
void 
Genome::initialize(double u, double sg, double opt, double env) {
  mu = u;
  sig = sg;
  optimum = opt;
  environmental_noise = env;
  baseline = 0;
  mutation_count = 0;
  return;
}

/* clear a genome of all derived mutations */
void
Genome::clear(void) {
  for (vector<mutation_loc>::iterator it = mutant_sites.begin(); it != mutant_sites.end(); it++) {
    pop->sites[*it].set_genotype(individual, homozygote_ancestral);
  }
  mutant_sites.clear();
#ifdef EXTRA_CHECKS
  for (int s=0; s < (int)pop->sites.size(); s++) {
    if (pop->sites[s][individual] != 0) 
      throw SimError(0, "site %d should be clear in individual %d", s, individual);
  }
#endif /* EXTRA_CHECKS */
  return;
}

/* compute the combined genotype contribution to phenotype */
double 
Genome::genvalue(void) {
  double sum = baseline;

  for (vector<mutation_loc>::iterator it = mutant_sites.begin(); it != mutant_sites.end(); it++)
    sum += pop->sites[*it][individual] * pop->sites[*it].effect;

  return sum;
}

/* update the phenotype */
double 
Genome::update_phenotype(void) {
  phenotype = genvalue() + ran1()*environmental_noise;

  return phenotype;
}

/* update the fitness, assumes the phenotype is up to date */
double 
Genome::update_fitness(void) {
  fitness = exp( -(phenotype-optimum)*(phenotype-optimum)/sig );
  /* keep track of the maximum fitness */
  if (pop->max_fitness < fitness)
    pop->max_fitness = fitness;
  return fitness;
}

/* set a new optimum, used when the environment changes */
void 
Genome::new_optimum(double opt) { optimum = opt; }

/* replace this genome with a recombined product of two other genomes */
void 
Genome::mate(Genome *mother, Genome *father) {
  /* go through the mother's mutations first */
  for (vector<mutation_loc>::iterator it=mother->mutant_sites.begin(); it != mother->mutant_sites.end(); it++) {
    /* if it's not a heterozygote or if we happen to sample the derived allele,
     * then we start the genotype out as a heterozygote, i.e., it has 
     * inherrited one derived allele from the mother. 
     * 
     * HAPLOID: there are only two genotypes 0 and 1 (het). The 0 genotype doesn't
     * store this site in the maternal genome (because there's no derived allele, 
     * and thus the only time we get here, is if it's genotype 1 (het). So it's 
     * not possible for the genotype to not be heterozygote, and thus the first 
     * logical operatio nbelow is always false, and so we always inherit the maternal 
     * derived allele with probability 1/2. Thus adding the haploid case doesn't 
     * change this expression form the original diploid implementation. */
    if (mother->pop->sites[*it][mother->individual] != heterozygote || ran1() < 0.5) {
      pop->sites[*it].set_genotype(individual, heterozygote); /* set the genotype in the Site object */
      mutant_sites.push_back(*it); /* add this site to the list of ones with derived alleles */
    }
  }

  /* now go through the father's mutations */
  for (vector<mutation_loc>::iterator it=father->mutant_sites.begin(); it != father->mutant_sites.end(); it++) {
    enum genotype child_genotype = pop->sites[*it][individual];
    /* see if the child already has one derived allele */
    if (child_genotype > homozygote_ancestral) {
      /* for the haploid case, if we already have a derived allele, we can't inherit another */
      if (Site::ploidy_level == haploid) continue;
      /* if the father is not a heterozygote or we happen to sample the derived
       * allele, make the childe homozygote-derived */
      if (father->pop->sites[*it][father->individual] != heterozygote || ran1() < 0.5)
        child_genotype = homozygote_derived;
    } else {
      /* if the mother has a derived allele (equivalent to het in the diploid encoding), and 
       * we didn't copy it, then we necessarily copy the father's allele (which must be 
       * derived, otherwise we wouldn't be here, as we're iterating through the father's 
       * derived alleles. If the mother does not have a derived allele, then we copy that 
       * ancestral allele with probability 0.5 (in which case no derived allele needs to be 
       * noted, or we copy the paternal one which is necessarily derived (again, because 
       * we're here) */
      if (Site::ploidy_level == haploid) {
        if (mother->pop->sites[*it][mother->individual] == heterozygote) {
          child_genotype = heterozygote;
        } else {
          if (ran1() < 0.5) child_genotype = heterozygote;
        }
      } else {
        if (father->pop->sites[*it][father->individual] != heterozygote || ran1() < 0.5)
          child_genotype = heterozygote;
      }
      if (child_genotype > homozygote_ancestral)
        mutant_sites.push_back(*it); /* add this site to the list of ones with derived alleles */
    }
    if (child_genotype > homozygote_ancestral)
      pop->sites[*it].set_genotype(individual, child_genotype); /* set the genotype in the Site object */
  }
  
  mutate_genome();
  update_phenotype();
  update_fitness();
  return;
}

/* print a genome's mutations */
ostream& 
operator<<(ostream &s, Genome &g) {
  for (vector<mutation_loc>::iterator it=g.mutant_sites.begin(); it != g.mutant_sites.end(); it++)
    s << " " << g.pop->sites[*it].id << ":" << g.pop->sites[*it][g.individual];
  return s;
}

void 
Genome::mutate_genome(void) {
  /* draw a poisson number of mutations */
  int num_muts;
  if (Site::ploidy_level == diploid) {
    num_muts = poidev(2.0*mu);
  } else {
    num_muts = poidev(1.0*mu);
  }
  for (int i = 0; i < num_muts; i++) 
    mutate_site();
  mutation_count += num_muts;
  return;
}

/* Remove all mutations corresponding to derived alleles at site loc.
 * This function is assuming that someone else is doing the book-keeping 
 * for sites[] */
void
Genome::purge_site(mutation_loc loc) {
  vector<mutation_loc>::iterator new_end = remove(mutant_sites.begin(), mutant_sites.end(), loc);
  if (new_end != mutant_sites.end())
    mutant_sites.erase(new_end, mutant_sites.end());
}

/* Check the integrity of an individual's genome */
void
Genome::check(void) {
  int derived_alleles_1 = 0;
  int derived_alleles_2 = 0;
  for (vector<mutation_loc>::iterator it=mutant_sites.begin(); it != mutant_sites.end(); it++) {
    derived_alleles_1 += pop->sites[*it][individual];
  }
  for (int i=0; i < (int)pop->sites.size(); i++) {
    derived_alleles_2 += pop->sites[i][individual];
  }
  if (derived_alleles_1 != derived_alleles_2)
    throw SimError(0, "derived allele counts don't match: %d != %d", derived_alleles_1, derived_alleles_2);
}

/************************************************** 
 * Genome implementation for infinite sites model *
 **************************************************/

/* class (static) variables */
std::valarray<double> GenomeInfiniteSites::effect_sizes;
std::valarray<double> GenomeInfiniteSites::effect_probabilities;

GenomeInfiniteSites::GenomeInfiniteSites(Population *p, int indiv) : Genome(p, indiv) { 
}

/* Set up the class globals. This allows us to keep track of some constants
 * and do some book keeping regarding mutations */
void 
GenomeInfiniteSites::setup_effect_probabilities(valarray<double> &ep, 
    valarray<double> &es) {
  effect_sizes.resize(es.size());
  effect_sizes = es;
  effect_probabilities.resize(ep.size());
  effect_probabilities = ep;
  return;
}

/* sample an effect size from the effect size probability distribution */
double 
GenomeInfiniteSites::sample_effect_size(void) {
  double max = effect_probabilities.max();
  int r;
  double sign = 1.0;

  /* The way I've implemented the infinite sites model is a bit different than
   * the Barton model. In the Barton model, genotype 0,1,2 corresponds to -a,0,a
   * But for the infinite sites model, this doesn't make sense because the 
   * change should be relative to zero. If this were not the case, then introducing
   * a new mutation changes the fitnesses of the rest of the individuals. 
   * I want to allow the new mutation to either increase or decrease the phenotype 
   * relative to the background on which it arose. So I need to also randomly 
   * pick the sign of the effect. So with probability 0.5, the effects are 0,a,2a 
   * and with probability 0.5 they are 0,-a,-2a. HAPLOID: for haploid, the effect 
   * sizes are 0,a or 0,-a. This difference results from genotypes only having values
   * 0,1 instead of 0,1,2. */
  if (ran1() < 0.5) sign = -1.0;

  while (1) {
    r = (int)(ran1() * effect_sizes.size());
    if (ran1()*max < effect_probabilities[r]) 
      return effect_sizes[r] * sign;
  }
}

/* mutate a new site */
void 
GenomeInfiniteSites::mutate_site(void) {
  double e = sample_effect_size();
  mutate_site(Population::create_site(e));
  return;
}

/* mutate a specific site. direction has default value 'up' */
void 
GenomeInfiniteSites::mutate_site(mutation_loc loc, double direction) {
  if (direction != up) throw SimError("infinite sites can only mutate up");
  switch (pop->sites[loc][individual]) {
    case homozygote_ancestral:
      mutant_sites.push_back(loc);
      pop->sites[loc].set_genotype(individual, heterozygote);
      break;
    case heterozygote:
      if (Site::ploidy_level == haploid)
        throw SimError("haploid populations can't mutate already mutated sites.\n");
      pop->sites[loc].set_genotype(individual, homozygote_derived);
      break;
    case homozygote_derived:
      throw SimError("can't mutate a homozygote-derived site");
    default:
      throw SimError("invalid genotype");
  }
  return;
}

/************************************************ 
 * Genome implementation for finite sites model *
 ************************************************/

GenomeFiniteSites::GenomeFiniteSites(Population *p, int indiv) : Genome(p, indiv) { 
}

/* mutate a random site */
void 
GenomeFiniteSites::mutate_site(void) {
  mutate_site( (mutation_loc)floor(ran1()*Population::num_loci) );
  return;
}

/* mutate a specific site. u can be passed to force it to mutate in a certain 
 * direction (the enum 'up' or 'down' can be used). This only matters in the 
 * case it's a heterozygote. When current genotype is homozygote only one 
 * direction is possible, and u is ignored. I only pass u when initially 
 * creating homozygote-derived genotypes, which I create by mutating a site 
 * twice 'up'. By default u has default value ran1() */
void 
GenomeFiniteSites::mutate_site(mutation_loc loc, double u) {
  /* Here I don't need to wory about the haploid case, because that's not yet 
   * supported for the finite sites model */
  switch (pop->sites[loc][individual]) {
    case homozygote_ancestral:
      mutant_sites.push_back(loc);
      pop->sites[loc].set_genotype(individual, heterozygote);
      break;
    case heterozygote:
      if (u < 0.5) {
        pop->sites[loc].set_genotype(individual, homozygote_ancestral);
        /* this find() is inefficient because it needs to search the whole 
         * vector, luckily the vectors are short (unless there are many 
         * intermediate-frequency sites) and mutations do not occur that 
         * often (compared to other operations) */
        vector<mutation_loc>::iterator x = find(mutant_sites.begin(), mutant_sites.end(), loc);
        if (x == mutant_sites.end()) {
          cerr << "looking for loc " << loc << " in individual " << individual << ", dumping list:";
          for (x = mutant_sites.begin(); x != mutant_sites.end(); x++) 
            cerr << " " << *x;
          cerr << endl;
          throw SimError("failed to find mutant site");
        }
        mutant_sites.erase(x);
      } else {
        pop->sites[loc].set_genotype(individual, homozygote_derived);
      }
      break;
    case homozygote_derived:
      pop->sites[loc].set_genotype(individual, heterozygote);
      break;
    default:
      throw SimError("invalid genotype");
  }
  return;
}

/* END */
