#ifndef __STATISTIC_H__
#define __STATISTIC_H__

#include "population.h"
#include <map>
#include <string>

/* The Statistic class provides an interface to turn on printing of statistics and a means to check if a statistic should be printed. It doesn't actually print anything. For a function to 'implement' the Statistic interface, it just needs to abide by checking whether printing has been activated and only do anyting when it's been activated */
class Statistic {
public:
  Statistic() { }
  ~Statistic() { }
  static void initialize_defaults(void);
  static void activate(const char *key);
  static void deactivate(const char *key);
  static bool is_activated(const char *key);
  static void deactivate_all(void);

  static std::map<std::string,bool> directory;
};

#endif /* __STATISTIC_H__ */
