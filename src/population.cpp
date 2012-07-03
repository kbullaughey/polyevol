#include <vector>
#include <valarray>
#include <queue>
#include <iostream>

#include "population.h"
#include "genome.h"
#include "error_handling.h"
#include "sim_rand.h"
#include "statistic.h"

using std::valarray;
using std::cout;
using std::endl;
using std::ostream;
using std::vector;
using std::queue;
using std::map;

/* storage for class variables */
int Population::num_loci = 0;
queue<int> Population::lost;
bool Population::initialized = false;
int Population::popsize;
Model Population::sites_model;
vector<Population*> Population::pop_views;
int Population::generation = 0;

/* static storage used by population-level statistics */
map<double,int> Population::fixations;
vector<int> Population::visits;
RunningMean *Population::delta_p_first_moment;
RunningMean *Population::delta_p_second_moment;

/* Initialize the class variables of Population */
void Population::initialize(int N, Model m) {
  if (initialized) throw SimError("Population class already initialized");
  popsize = N;
  sites_model = m;
  initialized = true;
  if (Statistic::is_activated("visits")) {
    if (Site::ploidy_level == diploid) {
      visits = vector<int>(2*N-1, 0);
    } else {
      visits = vector<int>(N-1, 0);
    }
  }
  if (Statistic::is_activated("p_moments")) {
    int bins;
    if (Site::ploidy_level == diploid) {
      bins = 2*N;
    } else {
      bins = N;
    }
    /* Note, this memory is not freed until program exit */
    delta_p_first_moment = new RunningMean(bins);
    delta_p_second_moment = new RunningMean(bins);
  }
}

/* Create a population */
Population::Population(void) {
  /* populations can't be added after initialize has been called */
  if (!initialized) throw SimError("Population class must be initialized");

  /* no fitness can be lower than zero, this gets updated in by Genome class 
   * each time a fitness is updated */
  max_fitness = 0;

  /* allocate genomes of this populations */
  for (int i=0; i < popsize; i++) {
    if (sites_model == infinite_sites) {
      genomes.push_back(new GenomeInfiniteSites(this, i));
    } else {
      genomes.push_back(new GenomeFiniteSites(this, i));
    }
  }
  /* add this population to the class list */
  pop_views.push_back(this);
}

/* set up sites based on initial genotypes */
void Population::setup_initial_genotypes(valarray<int> &hets, valarray<int> &homs) {
  if (hets.size() != homs.size()) throw SimError("len(hets) != len(homs)");
  if (hets.max() > popsize) throw SimError("too many heterozygotes");
  if (homs.max() > popsize) throw SimError("too many homozygotes");
  if ((hets+homs).max() > popsize) throw SimError("negative homozygote-ancestral");

  /* used to randomize individuals */
  valarray<int> ranout(popsize);

  mutation_loc loc;
  int ind;
  for (int k=0; k < (int)hets.size(); k++) {
    ranint(popsize,ranout);

    if (sites_model == infinite_sites)
      loc = create_site(GenomeInfiniteSites::sample_effect_size());
    else
      loc = (mutation_loc)k;

    /* haploid can't have homozygote_derived */
    if (Site::ploidy_level == haploid && homs[k] > 0)
      throw SimError("haploid can't have homozygote derived site");

    /* mutate heterozygote sites once */
    for (ind=0; ind < hets[k]; ind++) 
      genomes[ranout[ind]]->mutate_site(loc, up);
    /* mutate homozygote sites twice */
    for (; ind < hets[k]+homs[k]; ind++) {
      genomes[ranout[ind]]->mutate_site(loc, up);
      genomes[ranout[ind]]->mutate_site(loc, up);
    }
  }

  /* compute fitnesses */
  for (ind = 0; ind < popsize; ind++) {
    genomes[ind]->update_phenotype();
    genomes[ind]->update_fitness();
  }
  return;
}

/* create the next generation (this object) from the parent generation */
void Population::populate_from(const Population &parpop) {
  int mom, dad;

  /* clear out all the children's genomes */
  for (int off = 0; off < popsize; off++) 
    genomes[off]->clear();

  /* loop over the offspring, creating each by mating two parents sampled 
   * according to their fitnesses */
  for (int off = 0; off < popsize; off++) {
    while (1) {
      mom = (int)(popsize*ran1());
      if (ran1()*parpop.max_fitness <= parpop.genomes[mom]->fitness)
        break;
    }
    while (1) {
      dad = (int)(popsize*ran1());
      if (ran1()*parpop.max_fitness <= parpop.genomes[dad]->fitness)
        break;
    }
    /* have some sex */
    genomes[off]->mate(parpop.genomes[mom], parpop.genomes[dad]);
  }
}

/* Create a new site. If there are lost sites, reuse one of these. Either way 
 * the result site gets a new mutation ID */
mutation_loc
Population::create_site(double e) {
  mutation_id id = Site::next_unique_id++;
  mutation_loc loc;
  if (lost.size() > 0) {
    loc = lost.front();
    lost.pop();
    /* go through each population view and renew the old site for this new mutation */
    for (vector<Population*>::iterator it = pop_views.begin(); it != pop_views.end(); it++)
      (*it)->sites[loc].renew(e, id, generation);   
  } else {
    loc = num_loci++;
    /* go through each population view and create a new site object */
    for (vector<Population*>::iterator it = pop_views.begin(); it != pop_views.end(); it++) 
      (*it)->sites.push_back(Site(popsize, e, id, generation));
  }

  /* dump the site from one of the pop views so we have a record of its creation */
  if (Statistic::is_activated("mutation"))
    cout << "gen: " << generation << " " << pop_views[0]->sites[loc] << endl;
  return loc;
}

/* cleanup unused sites */
void
Population::purge_lost(void) { 
  for (mutation_loc loc=0; loc < (mutation_loc)num_loci; loc++) {
    /* check the site in this population to see if it's empty but not 
     * already made reusable */
    if (sites[loc].derived_alleles_count == 0 && !sites[loc].reusable) {
      /* loop through the two population views and make the site reusable */
      for (vector<Population*>::iterator pit = pop_views.begin(); pit != pop_views.end(); pit++)
        (*pit)->sites[loc].reusable = true; /* make the site reusable */
      /* record this site as having been lost */
      lost.push(loc);
      if (Statistic::is_activated("sojourn")) {
        cout << "gen: " << generation << " absorption loss site: " << sites[loc].id 
          << " sojourn: " << generation-sites[loc].generation_created 
          << " effect: " << sites[loc].effect << endl;
      }
    } else if (sites[loc].derived_alleles_count == Site::ploidy_level*popsize && !sites[loc].reusable) {
      /* dealing with a fixed site is more complicated because we need to remove
       * it from all genomes and adjust the baseline to reflect this sites now 
       * perminant effect */
      for (vector<Population*>::iterator pit = pop_views.begin(); pit != pop_views.end(); pit++) {
        for (vector<Genome*>::iterator git = (*pit)->genomes.begin(); 
            git != (*pit)->genomes.end(); git++) {
          /* search the genome for this site and remove it */
          (*git)->purge_site(loc);
        }
        (*pit)->sites[loc].reset(); /* sets all genotypes back to homozygous ancestral */
        (*pit)->sites[loc].reusable = true; /* make the site reusable */
      }
      /* adjust the genomic baseline to reflect the fixation */
      Genome::baseline += Site::ploidy_level*sites[loc].effect;
      fixations[sites[loc].effect]++;
      if (Statistic::is_activated("sojourn")) {
        cout << "gen: " << generation << " absorption fixation site: " << sites[loc].id 
          << " sojourn: " << generation-sites[loc].generation_created 
          << " effect: " << sites[loc].effect << endl;
      }
    }
  }
}

/* Return the pointer to the other population view  (there are only ever two) */
Population* 
Population::other_view(void) {
  if (this == pop_views[0]) {
    return pop_views[1];
  } else {
    return pop_views[0];
  }
}

void
Population::stat_increment_visits(void) {
  if (!Statistic::is_activated("visits")) return;
  for (mutation_loc loc=0; loc < sites.size(); loc++) {
    if (!sites[loc].reusable) 
      visits[sites[loc].derived_alleles_count-1]++;
  }
  return;
}

void
Population::stat_print_visits(void) {
  if (!Statistic::is_activated("visits")) return;
  cout << "visits:";
  for (int i=0; i<(int)visits.size(); i++)
    cout << " " << visits[i];
  cout << endl;
  return;
}

void
Population::stat_fixations(void) {
  if (!Statistic::is_activated("fixations")) return;
  cout << "gen: " << generation << " fixations:";
  for (map<double,int>::iterator i=fixations.begin(); i!=fixations.end(); i++) {
    cout << " " << i->first << "," << i->second;
  }
  cout << endl;
  return;
}

/* print out the frequencies of all the sites that have mutated so far */
void
Population::stat_frequency_summary(void) {
  if (!Statistic::is_activated("frequencies")) return;
  cout << "gen: " << generation << " freqs:";
  for (mutation_loc loc=0; loc < sites.size(); loc++) {
    /* print only sites that haven't been recorded as lost */
    if (!sites[loc].reusable) {
      double f = sites[loc].frequency();
      if (f < 1.0) 
        cout << " " << sites[loc].id << ":" << f;
    }
  }
  cout << endl;
  return;
}

/* print out the number of segregating sites for each effect size */
void
Population::stat_segsites(void) {
  if (!Statistic::is_activated("segsites")) return;
  cout << "gen: " << generation << " segsites:";
  map<double,int> counts;
  for (mutation_loc loc=0; loc < sites.size(); loc++) {
    /* record only sites that haven't been lost */
    if (!sites[loc].reusable) counts[sites[loc].effect]++;
  }
  for (map<double,int>::iterator i=counts.begin(); i != counts.end(); i++) {
    cout << " " << i->first << "," << i->second;
  }
  cout << endl;
  return;
}

/* print out the phenotype mean and variance */
void
Population::stat_phenotype_summary(void) {
  if (!Statistic::is_activated("phenotype")) return;

  double sum, sumsq, p;
  sum = sumsq = 0.0;
  for (int ind=0; ind < popsize; ind++) {
    p = genomes[ind]->phenotype;
    sum += p;
    sumsq += p*p;
  }
  cout << "gen: " << generation << " pheno: " << sum/popsize 
    << " " << sumsq/popsize - (sum/popsize)*(sum/popsize) << endl;
  return;
}

/* Update the estimates for the first and second moment of the change in allele
 * frequency. Notice that the quantity we're keeping track of is not actually 
 * the change in allele frequency, it's the change in derived allele counts. 
 * It's left up to the user to divide by the population size */
void
Population::stat_update_p_moments(void) {
  if (!Statistic::is_activated("p_moments")) return;

  int delta, current_p, previous_p;
  for (mutation_loc loc=0; loc < sites.size(); loc++) {
    /* We only consider sites that are currently in use (sites can be waiting 
     * to be resused if they've been lost from the population) */
    if (!sites[loc].reusable) {
      /* Compute the change in allele frequency between this population view and 
       * the other view, which will be parent generation when this function is 
       * called. In the case when the site is new, the derived_alleles_count in
       * the parent generation (accessible via other_view) will be zero. */
      current_p = sites[loc].derived_alleles_count;
      previous_p = other_view()->sites[loc].derived_alleles_count;
      delta =  current_p - previous_p;
      delta_p_first_moment->post(current_p, (double)delta);
      delta_p_second_moment->post(current_p, pow((double)delta, 2.0));
    }
  }
}

/* Print out the first and second moments for the change in allele frequency */
void
Population::stat_print_p_moments(void) {
  if (!Statistic::is_activated("p_moments")) return;
  cout << "delta_p_first_moment: " << *delta_p_first_moment << endl;
  cout << "delta_p_second_moment: " << *delta_p_second_moment << endl;
}

/* print out the segregating sites of all individuals in the population */
ostream& 
operator<<(ostream &s, const Population &p) {
  Genome *g;
  s << "Dumping population (" << &p << "):" << endl;
  for (int i=0; i < p.popsize; i++) {
    g = p.genomes[i];
    s << " [" << i << "] " << "phenotype: " << g->phenotype << " fitness: " << g->fitness << " genotypes: " << *g << endl;
  }
  s << "there are " << p.pop_views.size() << " views" << endl;
  return s;
}

/* END */
