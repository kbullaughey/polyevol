#include <map>
#include <string>

#include "statistic.h"
#include "error_handling.h"

using std::map;
using std::string;

/* populate the directory with the valid statistics, and their defaults */
map<string,bool> Statistic::directory;

/* given we only have class variables we use an initializer */
void Statistic::initialize_defaults(void) {
  directory[string("frequencies")] = true;
  directory[string("phenotype")] = true;
  directory[string("burnin")] = true;
  directory[string("sojourn")] = true;
  directory[string("mutation")] = true;
  directory[string("visits")] = false;
}

/* turn off all the statistics */
void Statistic::deactivate_all(void) {
  for (map<string,bool>::iterator it=directory.begin(); it != directory.end(); it++)
    (*it).second = false;
}

/* activate statistic by key */
void Statistic::activate(const char *key) {
  string k(key);
  if (directory.count(key) == 0)
    throw SimError(0, "invalid statistic: %s", key);
  directory[key] = true;
}

/* deactivate a statistic by key */
void Statistic::deactivate(const char *key) {
  string k(key);
  if (directory.count(key) == 0)
    throw SimError(0, "invalid statistic: %s", key);
  directory[key] = false;
}

/* used to check if a statistic should be printed. All statistics should call 
 * this first to see if they should be printed */
bool Statistic::is_activated(const char *key) {
  string k(key);
  if (directory.count(key) == 0)
    throw SimError(0, "invalid statistic: %s", key);
  return directory[key];
}

/* END */
