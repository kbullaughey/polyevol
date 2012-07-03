#ifndef __POPULATION_H__
#define __POPULATION_H__

#include <valarray>
#include <ostream>
#include <vector>
#include <map>

#include "common.h"
#include "genome.h"
#include "running_mean.h"

class Population {
public:
  Population(void);
  ~Population() { }
  void setup_initial_genotypes(std::valarray<int> &hets, std::valarray<int> &homs);
  void stat_frequency_summary(void);
  void stat_phenotype_summary(void);
  void stat_increment_visits(void);
  void stat_fixations(void);
  void stat_segsites(void);
  static void stat_print_visits(void);
  void record_genotype(int indiv, mutation_loc loc, genotype g);
  void populate_from(const Population &parpop);
  void purge_lost(void);
  friend std::ostream& operator<<(std::ostream &s, const Population &p);

  /* I need a few class functions */
  static mutation_loc create_site(double e);
  static void initialize(int N, Model m);

  /* maximum fitness in the population, used for rejection sampling */
  double max_fitness;

  /* this is the length of the sites vectors, which are all the same length */
  static int num_loci;
  static int generation;

  /* I keep records in two ways: A list of genomes, each of which contains 
   * the loci that have derived alleles in that individual, and a list of
   * sites which contain the genotypes of all the individuals for that site.
   * I need two storage containers to make book-keeping efficient */
  std::vector<Site> sites;
private:
  std::vector<Genome*> genomes;

  /* stuff below here is common to all populations */

  /* A Population object, is actually a view onto a single population. I use
   * just two views, one for the current generation (parents) and one for the 
   * subsequent generation (offspring). However, these are stored as separate
   * population objects. But certain things, like the the vector of Site objects
   * need to be kept in sync, i.e., when a site is added, it should be added to
   * all Population objects, as these are really views into the same population
   * with the same segregating variation. For these reasons, I keep a class 
   * variable, pop_views, that allows me to modify all populations when 
   * necessary */
  static std::vector<Population*> pop_views;

  /* I do my own book keeping of sites that have been lost, so I can reuse 
   * them. This prevents me from allocating new memory every time a site drifts 
   * out of existence. It also allows me to keep existing sites in the same 
   * locations in the Site vector, as this is how I find them. As sites are 
   * lost, this opens up holes in the vector, so I reuse them, keeping track 
   * of which ones I can reuse in this queue container. */
  static std::queue<int> lost;

  static bool initialized;
  static int popsize;
  static Model sites_model;

  /* use by statistics */
  static std::vector<int> visits;
  static std::map<double,int> fixations;
  static RunningMean *delta_p_first_moment;
  static RunningMean *delta_p_second_moment;
};

#endif /* __POPULATION_H__ */
