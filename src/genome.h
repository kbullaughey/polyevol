#ifndef __GENOME_H__
#define __GENOME_H__

#include <iostream>
#include <valarray>
#include <vector>
#include <queue>

#include "site.h"
#include "sim_rand.h"

/************************* 
 * Genome abstract class *
 *************************/

typedef unsigned int mutation_loc;
class Population;

/* used to mutate a genotype:
 *  - up (add a derived allele) 
 *  - down (remove a derived allele) */
enum { down, up };

class Genome {
public:
  /* public member functions */
  Genome(Population *p, int indiv);
  virtual ~Genome() { }
  double genvalue(void);
  void mate(Genome *mother, Genome *father);
  double update_phenotype(void);
  double update_fitness(void);
  void clear(void);
  void mutate_genome();
  virtual void mutate_site(mutation_loc loc, double direction) = 0;
  virtual void mutate_site(void) = 0;

  /* public class function */
  static void new_optimum(double);
  static void initialize(double u, double sig, double opt, double env);

  /* operators */
  friend std::ostream& operator<<(std::ostream &s, Genome &g);

  /* genome-specific public data */
  double fitness;
  double phenotype;

  /* Since I only keep track of derived alleles, the baseline genvalue needs 
   * to be sum_i(-a_i) so that a genome that is entirely homozygotes ancestral,
   * has the correct genvalue */
  static double baseline;

protected:
  /* list of sites containing one or more derived allele in this individual */
  std::vector<mutation_loc> mutant_sites;

  /* index indicating which individual this is in the ordered population */
  int individual;

  /* each Genome is associated with a population, so it can use the population's
   * lookup table and other data */
  Population *pop;
  
  /* these static member variables store global class information */
  static double mu;
  static double environmental_noise;
  static double optimum;
  static double sig;
};

/************************************************** 
 * Genome implementation for infinite sites model *
 **************************************************/

class GenomeInfiniteSites : public Genome {
public:
  /* public member functions */
  GenomeInfiniteSites(Population *p, int indiv);
  ~GenomeInfiniteSites() { }
  void mutate_site(void);
  void mutate_site(mutation_loc loc, double direction = up);
  /* public class functions */
  static void setup_effect_probabilities(std::valarray<double> &ep, 
    std::valarray<double> &es);
  static double sample_effect_size(void);

private:
  /* private class variables */
  static std::valarray<double> effect_sizes;
  static std::valarray<double> effect_probabilities;
};

/************************************************ 
 * Genome implementation for finite sites model *
 ************************************************/

class GenomeFiniteSites : public Genome {
public:
  /* public member functions */
  GenomeFiniteSites(Population *p, int indiv);
  ~GenomeFiniteSites() { }
  void mutate_site(void);
  void mutate_site(mutation_loc loc, double direction = ran1());
};

#endif /* __GENOME_H__ */

