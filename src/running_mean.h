#ifndef __RUNNING_MEAN_H__
#define __RUNNING_MEAN_H__

#include <vector>
#include <ostream>

/* The RunningMean class provides a generic way to keep track of a set of related running mean statistics. I use this for keeping track of the first and second moments of the change in allele frequency as a function of the present frequency. As implemented, the collection of means are a zero-indexed array. */
class RunningMean {
public:
  RunningMean(int size);
  ~RunningMean() { }

  void post(int i, double value);
  void post(double value);
  int size(void) const;
  int count(int i) const;

  double operator[] (const int i) const;
  friend std::ostream& operator<<(std::ostream &s, const RunningMean &m);
  
private:
  std::vector<double> means;    /* keep track of the mean */
  std::vector<double> counts;   /* keep track of how many samples comprise this estimate */
  int _size;                       /* number of means to keep track of */
};

#endif /* __RUNNING_MEAN_H__ */

