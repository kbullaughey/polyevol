#ifndef __COMMON_H__
#define __COMMON_H__

#include <valarray>
#include <string>

enum Model {unspecified, infinite_sites, finite_sites };

void print_double_vector(std::valarray<double> &x, const char *label);
std::string& print_r_vector(const std::valarray<int> &x, const char *label, std::string &s);
std::string& print_r_vector(const std::valarray<double> &x, const char *label, std::string &s);

#endif /* __COMMON_H__ */
