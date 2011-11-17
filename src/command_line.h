#ifndef __COMMAND_LINE_H
#define __COMMAND_LINE_H

#include <iostream>
#include <valarray>

#include "common.h"

class Args {
public:
   /* public member functions */
   Args(int argc, char *argv[]);
   ~Args() { }
   friend std::ostream& operator<<(std::ostream &s, const Args &a);
   
   /* command line parameters */
  int popsize;
  unsigned int rand_seed;
  int burnin;                               /* generations of burnin */
  double mu;                                /* mutation rate */
  double s;                                 /* Barton s parameter */
  double env;                               /* environmental phenotypic variance */
  Model sites_model;                        /* infinite or finite sites models */
  std::valarray<double> effect_sizes;       /* vector of loci effect sizes */
  std::valarray<double> opts;               /* vector of the optima for each epoch */
  std::valarray<int> times;                 /* times (in generations) when the epochs end */
  std::string cmd;

  /* for fixed number of loci model */
  std::valarray<int> loci_counts;           /* vector of loci count for each effect size */

  /* for infinite sites model */
  std::valarray<double> effect_probabilities;  /* vector of probabilities of effect sizes */
  int initial_num_loci;                     /* number of loci, and lines of input */

private:
   /* private functions */
   bool has_option(char *);
   void fix_negatives(char *str);
};

#endif /* __COMMAND_LINE_H */
