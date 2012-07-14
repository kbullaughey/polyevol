#include "running_mean.h"
#include "error_handling.h"

using std::vector;
using std::ostream;

/* Constructor. Allocate storage and set the size */
RunningMean::RunningMean(int size) : means(size), counts(size) { 
  _size = size;
}

/* For the special case where we're only storing a single running mean in this 
 * object, add a sample to the mean at index 0 */
void RunningMean::post(double s) {
  post(0, s);
}

/* add a sample to the mean at index i */
void RunningMean::post(int i, double s) {
  if (i < 0) throw SimError("cannot tolerate a negative index");
  if (i >= _size) throw SimError("index exceeds size of RunningMean instance"); 
  means[i] = ((means[i] * counts[i]) + s) / (counts[i] + 1.0);
  counts[i] += 1.0;
}

/* Get the mean at index i */
double RunningMean::operator[](const int i) const {
  if (i < 0) throw SimError("cannot tolerate a negative index");
  if (i >= _size) throw SimError("index exceeds size of RunningMean instance"); 
  return means[i];
}

/* How many times the mean at index i was posted to */
int RunningMean::count(int i) const {
  if (i < 0) throw SimError("cannot tolerate a negative index");
  if (i >= _size) throw SimError("index exceeds size of RunningMean instance"); 
  return (int)counts[i];
}

/* Number of means we're keeping track of */
int RunningMean::size(void) const {
  return _size;
}

/* print out the segregating sites of all individuals in the population */
ostream& 
operator<<(ostream &s, const RunningMean &m) {
  /* save the original precision setting, so we can set it back at the end */
  std::streamsize old_precision = s.precision();
  s.precision(8);

  for (int i=0; i < m.size(); i++) 
    s << " " << m[i] << "," << m.count(i);

  s.precision(old_precision);
  return s;
}
/* END */
