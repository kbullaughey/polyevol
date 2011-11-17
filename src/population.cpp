#include <vector>
#include <valarray>
#include <queue>
#include <iostream>

#include "population.h"
#include "genome.h"
#include "error_handling.h"
#include "sim_rand.h"

using std::valarray;
using std::cout;
using std::endl;
using std::ostream;
using std::vector;
using std::queue;

/* storage for class variables */
int Population::num_loci = 0;
queue<int> Population::lost;
bool Population::initialized = false;
int Population::popsize;
Model Population::sites_model;
vector<Population*> Population::pop_views;

/* Initialize the class variables of Population */
void Population::initialize(int N, Model m) {
  popsize = N;
  sites_model = m;
  initialized = true;
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
    /* go through each population view and renew use the old site for this new mutation */
    for (vector<Population*>::iterator it = pop_views.begin(); it != pop_views.end(); it++)
      (*it)->sites[loc].renew(e, id);   
  } else {
    loc = num_loci++;
    /* go through each population view and create a new site object */
    for (vector<Population*>::iterator it = pop_views.begin(); it != pop_views.end(); it++) 
      (*it)->sites.push_back(Site(popsize, e, id));
  }

  /* dump the site from one of the pop views so we have a record of its creation */
  cout << pop_views[0]->sites[loc] << endl;
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
    }
  }
}


/* print out the frequencies of all the sites that have mutated so far */
void
Population::print_frequency_summary(void) {
  cout << "freqs:";
  for (mutation_loc loc=0; loc < sites.size(); loc++) {
    /* print only sites that haven't been recorded as lost */
    if (!sites[loc].reusable)
      cout << " " << sites[loc].id << ":" << sites[loc].frequency();
  }
  cout << endl;
  return;
}

/* print out the phenotype mean and variance */
void
Population::print_phenotype_summary(void) {
  double sum, sumsq, p;

  sum = sumsq = 0.0;
  for (int ind=0; ind < popsize; ind++) {
    p = genomes[ind]->phenotype;
    sum += p;
    sumsq += p*p;
  }
  cout << "pheno: " << sum/popsize 
    << " " << sumsq/popsize - (sum/popsize)*(sum/popsize) << endl;
  return;
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
