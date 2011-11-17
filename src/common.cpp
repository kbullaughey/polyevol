#include <iostream>
#include <valarray>
#include <string>
#include <sstream>

#include "common.h"

using std::cout;
using std::endl;
using std::valarray;
using std::string;
using std::stringstream;

void
print_double_vector(valarray<double> &x, const char *label) {
  cout << label;
  for (int i=0; i < (int)x.size(); i++) 
    cout << " " << x[i];
  cout << endl;
}
   
/* I print vectors that are designed to be read in directly to R using
 * some R code. That way the output can be evaluated as an R expression
 * to load the appropriate data. This allows for vectors of paremeters.
 * In order to account for both double and int, this is a template function */
template <typename T>
string&
print_r_vector_aux(const valarray<T> &x, const char *label, string &out) {
  stringstream s;
  s << label << "=c(";
  for (int i=0; i < (int)x.size(); i++) {
    if (i > 0) s << ",";
    s << x[i];
  }
  s << ")";
  out = s.str();
  return out;
}

/* the next two functions are wrappers to get around the restriction that 
 * template functions must be defined in every file in which they are used */
string& 
print_r_vector(const valarray<int> &x, const char *label, string &s) {
  return print_r_vector_aux<int>(x, label, s);
}
string& 
print_r_vector(const valarray<double> &x, const char *label, string &s) {
  return print_r_vector_aux<double>(x, label, s);
}



/* END */
