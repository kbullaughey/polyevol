#ifndef __SITE_H__
#define __SITE_H__

#include <valarray>

enum genotype { homozygote_ancestral, heterozygote, homozygote_derived };
typedef unsigned int mutation_id;

class Site {
public:
  Site(int N, double e, mutation_id sid);
  ~Site() { }
  void set_genotype(int i, genotype g);
  double frequency(void);
  int count(void);
  void renew(double e, mutation_id sid);

  /* operators */
  friend std::ostream& operator<<(std::ostream &o, Site &s);
  /* access a particular individual's genotype */
  genotype operator[](int i);

  /* effect size of a derived allele at this site */
  double effect;

  /* this should only be read publically, not altered */
  int derived_alleles_count;
  
  /* used to mark a site as free for reuse */
  bool reusable;

  /* used to keep track of this mutation's unique id, shouldn't be altered */
  mutation_id id;

  /* used to dish out unique IDs */
  static mutation_id next_unique_id;

private:
  /* genotypes for each individual in the population */
  std::valarray<char> genotypes;
};

#endif /* __SITE_H__ */
