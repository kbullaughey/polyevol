#include "running_mean.h"

using std::valarray;

/* Constructor. Allocate storage and set the size */
RunningMean::RunningMean(int size) : means(size), counts(size) { 
  this.size = size;
}

/* add a sample to the mean at index i */
void RunningMean::post(int i, double s) {
  if (i >= size) throw SimError("index exceeds size of RunningMean instance"); 
  means[i] = ((means[i] * counts[i]) + s) / (counts[i] + 1.0);
  counts[i] += 1.0;
}

/* END */
