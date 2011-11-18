#include <iostream>
#include <ostream>

#include "error_handling.h"
#include "site.h"

using std::ostream;

/* storage for (static) class variables */
mutation_id Site::next_unique_id = 0;

Site::Site(int N, double e, mutation_id sid) : genotypes(N) {
  effect = e;
  derived_alleles_count = 0;
  reusable = false;
  id = sid;
}

/* get a particular individual's genotype at this site */
genotype Site::operator[](int i) {
  return (genotype)genotypes[i];
}

/* used for assignment to a particular individual's genotype at this site */
void Site::set_genotype(int i, genotype g) {
  if (i > (int)genotypes.size()) 
    throw SimError(0, "invalid individual: %d", i);
  /* update the the allele count */
  derived_alleles_count += g - genotypes[i];
  genotypes[i] = g;
  return;
}

/* compute the frequency of derived alleles at this site */
double Site::frequency(void) {
  return derived_alleles_count / (2.0 * genotypes.size());
}

/* take a site that was lost and use it for a new mutation */
void Site::renew(double e, mutation_id sid) {
  if (derived_alleles_count != 0 || !reusable) 
    throw SimError("site cannot be reused");
  effect = e;
  reusable = false;
  id = sid;
}

/* print out information related to this site */
ostream& operator<<(ostream &o, Site &s) {
  o << "site: id: " << s.id << " effect: " << s.effect;
  return o;
}


/* END */

