#include "command_line.h"
#include "error_handling.h"
#include "common.h"

#include <getopt.h>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#include <map>
#include <valarray>
#include <stdlib.h>

/* identifiers for arguments, below 255 is for ascii characters */
#define RAND_SEED     300
#define BURNIN        301
#define LOCI          302
#define ENV           303
#define EFFECTS       304
#define OPTS          305
#define TIMES         306
#define EFFECT_PROBS  307
#define FREQS         308

using std::cerr;
using std::cin;
using std::endl;
using std::ostream;
using std::string;
using std::vector;
using std::valarray;
using std::map;

void strsplit(const string &s, vector<string> &res, char sep);
void valdouble_from_string(const string &s, valarray<double> &x);
void valdouble_from_string(const char *s, valarray<double> &x);
void valint_from_string(const string &s, valarray<int> &x);
void valint_from_string(const char *s, valarray<int> &x);

/* a map to make converting from model string to constant easier */
static map<string,Model> model_lookup;
string model_reverse_lookup[] = { string("unspecified"), string("infinite"), string("finite") };

static map<string,freq_input> freq_lookup;
string freq_reverse_lookup[2] = { string("file"), string("even") };


/* set default options */
Args::Args(int argc, char *argv[]) {

  model_lookup[string("infinite")] = infinite_sites;
  model_lookup[string("finite")] = finite_sites;
  freq_lookup[string("file")] = freqfile;
  freq_lookup[string("even")] = freqeven;

  /* defaults */
  popsize = 5000;
  burnin = popsize;
  mu = 0.0001;
  s = 0.01;
  env = 0.0;
  rand_seed = 0;
  sites_model = unspecified;
  freqin = freqfile;

  /* process all the arguments from argv[] */
  int c;

  /* re-create the original command line as best we can */
  string sep("");
  for (int i=0; i<argc; i++) {
    cmd += sep+argv[i];
    sep = " ";
  }

  /* the following is a kluge to get negative arguments to work. Basically 
   * if -X occurs at the beginning of the string, I replace -X with NX 
   * where X is in [0-9]. Then later, before sending the optarg string on 
   * to a function that expects negative numbers I change it back to -X. */
  int pos, pos2;
  for (int i=0; i<argc; i++) {
    string s(argv[i]);
    pos = s.find_first_of("0123456789.");
    if (pos == 1 && *(argv[i]) == '-') {
      *(argv[i]) = 'N';
    } else {
      /* perhaps this arg is of the form --option=-X */
      pos = s.find("=");
      if (pos != (int)string::npos) {
        pos2 = s.find_first_of("0123456789.", pos);
        if (pos2-pos == 2 && (argv[i])[pos2-1] == '-') {
          (argv[i])[pos2-1] = 'N';
        }  
      }
    }
  }

  /* we don't want option errors, as we detect our own */
  opterr = 1;

  while (1) {
    static struct option long_options[] = {
      {"seed",  required_argument, 0, RAND_SEED},
      {"popsize",  required_argument, 0, 'N'},
      {"burnin", required_argument, 0, BURNIN},
      {"loci", required_argument, 0, LOCI},
      {"mu", required_argument, 0, 'u'},
      {"s", required_argument, 0, 's'},
      {"env", required_argument, 0, ENV},
      {"effects", required_argument, 0, EFFECTS},
      {"opts", required_argument, 0, OPTS},
      {"times", required_argument, 0, TIMES},
      {"eprobs", required_argument, 0, EFFECT_PROBS},
      {"freqs", required_argument, 0, FREQS},
      {"model", required_argument, 0, 'm'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;
  
    /* For the abreviated options, the ':' indicates a parameter is required following the option */
    c = getopt_long(argc, argv, "N:u:s:m:", long_options, &option_index);

    char *end;
    /* Detect the end of the options. */
    if (c == -1) break;
    switch (c) {
      case '?': {
        throw SimUsageError("extraneous option");
        break;
      }
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0) break;

        /* every non-flag option should not have c==0 */
        throw SimUsageError("unrecognized option");
        break;

      case 'm':
        if (!has_option(optarg))
          throw SimUsageError("must specify model");
        if (model_lookup.count(string(optarg)) == 0)
          throw SimUsageError("invalid model");
        sites_model = model_lookup[string(optarg)];
        break;

      case FREQS:
        if (!has_option(optarg))
          throw SimUsageError("must specify starting frequencies");
        if (freq_lookup.count(string(optarg)) == 0)
          throw SimUsageError("invalid frequency strategy");
        freqin = freq_lookup[string(optarg)];
        break;

      case 'u':
        if (!has_option(optarg))
          throw SimUsageError("must specify mutation rate");
        mu = strtod(optarg, &end);
        if (optarg == end) 
          throw SimUsageError("non-numeric mutation rate");
        break;

      case 's':
        if (!has_option(optarg))
          throw SimUsageError("must specify selection parameter s");
        s = strtod(optarg, &end);
        if (optarg == end) 
          throw SimUsageError("non-numeric s");
        break;

      case 'N':
        if (!has_option(optarg))
          throw SimUsageError("must specify population size");
        popsize = strtoul(optarg, &end, 10);
        if (optarg == end) 
          throw SimUsageError("non-numeric population size");
        break;

      case RAND_SEED:
        if (!has_option(optarg))
          throw SimUsageError("must specify random number seed");
        rand_seed = strtoul(optarg, &end, 10);
        if (optarg == end) 
          throw SimUsageError("non-numeric random number seed");
        break;

      case BURNIN:
        if (!has_option(optarg))
          throw SimUsageError("must specify number of burnin steps");
        burnin = strtoul(optarg, &end, 10);
        if (optarg == end) 
          throw SimUsageError("non-numeric number of burnin steps");
        break;

      case ENV:
        if (!has_option(optarg))
          throw SimUsageError("must specify environmental variance");
        env = strtod(optarg, &end);
        if (optarg == end) 
          throw SimUsageError("non-numeric environmental variance");
        break;

      case LOCI:
        if (!has_option(optarg))
          throw SimUsageError("must specify number of loci of each effect size");
        fix_negatives(optarg);
        valint_from_string(optarg, loci_counts);
        break;

      case EFFECTS:
        if (!has_option(optarg))
          throw SimUsageError("must specify effects");
        fix_negatives(optarg);
        valdouble_from_string(optarg, effect_sizes);
        break;

      case OPTS:
        if (!has_option(optarg))
          throw SimUsageError("must specify optima");
        fix_negatives(optarg);
        valdouble_from_string(optarg, opts);
        break;

      case TIMES:
        if (!has_option(optarg))
          throw SimUsageError("must specify times of optima shifts");
        fix_negatives(optarg);
        valint_from_string(optarg, times);
        break;

      case EFFECT_PROBS:
        if (!has_option(optarg))
          throw SimUsageError("must specify effect probabilities");
        fix_negatives(optarg);
        valdouble_from_string(optarg, effect_probabilities);
        break;

      default:
        char o[10];
        snprintf(o, 10, "%d", c);
        string s(o);
        throw SimUsageError("invalid command option" + s);
    }
  }

  /* throw exception if there are any extra options */
  if (optind < argc) {
    throw SimUsageError("unrecognized additional command options");
  }

  /* first check to make sure we give a model */
  if (sites_model == unspecified) 
    throw SimUsageError("must explicitly give a sites model (infinite/finite)");

  /* require certain parameters to be present, as there is no default */
  if (effect_sizes.size() == 0) throw SimUsageError("must specify effects");
  if (opts.size() == 0) throw SimUsageError("must specify opts");
  if (times.size() == 0) throw SimUsageError("must specify times");

 /* some minimal parameter checking for consistency */
 if (opts.size() != times.size()) throw SimUsageError("opts and times must be same length");

  /* check model-specific configuration */
  switch (sites_model) {
    case infinite_sites:
      if (effect_probabilities.size() == 0) {
        /* if there are no probabilities given and there's only one effect size, 
         * we know 100% are this effect size */
        if (effect_sizes.size() == 1) {
          effect_probabilities.resize(1);
          effect_probabilities[0] = 1.0;
        } else {
          throw SimUsageError("must specify effect size probabilities");
        }
      } 
      if (effect_probabilities.size() != effect_sizes.size()) {
        throw SimUsageError("effect sizes and effect probabilities must be same length");
      }
      if (loci_counts.size() != 1)
        throw SimUsageError("infinite sites model takes exactly 1 loci count");
      if (loci_counts.min() < 0)
        throw SimUsageError("negative number of loci"); 
      break;
    case finite_sites:
      if (loci_counts.size() == 0) throw SimUsageError("must specify loci count(s)");
      if (loci_counts.size() != effect_sizes.size()) 
        throw SimUsageError("loci and effects must be same length");
      if (loci_counts.min() <= 0)
        throw SimUsageError("non-positive number of loci"); 
      break;
    default:
        throw SimUsageError("invalid sites model");
  }
  nloci = (int)loci_counts.sum();

  /* initialize the random number generator */
  srand48(rand_seed);
  return;
}

/* negative numbers can look like parameters that begin with '-' signs. So I 
 * convert all these to the character 'N' and then here, if the string begins 
 * with N, change this back to '-' */
void Args::fix_negatives(char *str) {
  if (str[0] == 'N') str[0] = '-';
  int pos;
  string s(str);
  pos = s.find("=");
  if (pos != (int)string::npos) 
    if ((int)s.length() > pos+1 && s[pos+1] == 'N') str[pos+1] = '-';
}

/* make sure optarg doesn't point to the next option */
bool Args::has_option(char *optarg) {
  return !((optarg == 0) || (strlen(optarg) > 1 && optarg[0] == '-'));
}

/* print the arguments to the screen */
ostream& operator<<(ostream &s, const Args &a) {
  s << "cmd: " << a.cmd << endl;
  s << "params:"
    << " seed=" << a.rand_seed
    << " popsize=" << a.popsize 
    << " mu=" << a.mu
    << " s=" << a.s
    << " env=" << a.env
    << " model=\"" << model_reverse_lookup[a.sites_model] << "\""
    << " freqs=\"" << freq_reverse_lookup[a.freqin] << "\""
    << " burnin=" << a.burnin;

  string tmp;
  s << " " << print_r_vector(a.effect_sizes, "effects", tmp);
  s << " " << print_r_vector(a.opts, "opts", tmp);
  s << " " << print_r_vector(a.times, "times", tmp);
  if (a.sites_model == finite_sites)
    s << " " << print_r_vector(a.loci_counts, "loci", tmp);
  if (a.sites_model == infinite_sites)
    s << " " << print_r_vector(a.effect_probabilities, "eprobs", tmp);
  return s;
}

/* I use this function to fetch the next allele frequency. In addition to 
 * passing these as standard input, I'm implementing a few stock ways that 
 * these can be started. That way execution is not dependent on some file 
 * that may or may not be subsequently deleted or forgotten which was was 
 * used.
 * This function returns negative when it's out of frequencies. */
double Args::get_initial_frequency(void) {
  static double delivered = 0;
  double f;
  switch (freqin) {
    case freqfile:
      cin >> f;
      if (cin.eof() != 0) return -1;
      break;
    case freqeven:
      if (delivered == nloci) return -1;
      f = (delivered+1) / (nloci+1);
      break;
    default:
      throw SimError("invalid frequency initiation strategy");
  }
  delivered++;
  return f;
}


/* split a string by sep filling the vector */
void strsplit(const string &s, vector<string> &res, char sep) {
  int newpos = 0, oldpos = 0;
  while (1) {
    newpos = s.find(sep, oldpos);
    if (newpos == (int)string::npos) break;
    if (newpos-oldpos <= 0) throw SimError("nothing between separators");
    res.push_back(s.substr(oldpos, newpos-oldpos));
    oldpos = newpos+1;
  }
  res.push_back(s.substr(oldpos, s.length()-oldpos));
}

/* allocate a valarray for doubles, filling it from the comma-separated str */
void valdouble_from_string(const string &s, valarray<double> &x) {
  vector<string> pieces;
  strsplit(s, pieces, ',');
  if (pieces.size() == 0) throw SimError("error parsing comma-separated list");
  if (x.size() != pieces.size()) x.resize(pieces.size());
  unsigned int k = 0;
  for (vector<string>::iterator i = pieces.begin(); i != pieces.end(); i++) {
    char *end;
    x[k++] = strtod(i->c_str(), &end);
    if (i->c_str() == end) throw SimError("error in comma-separated list");
  }
}

/* same as above but for char* */
void valdouble_from_string(const char *s, valarray<double> &x) {
  string tmp(s);
  valdouble_from_string(tmp, x);
}

/* allocate a valarray for base-10 integers, filling it from the comma-separated str */
void valint_from_string(const string &s, valarray<int> &x) {
  vector<string> pieces;
  strsplit(s, pieces, ',');
  if (pieces.size() == 0) throw SimError("error parsing comma-separated list");
  if (x.size() != pieces.size()) x.resize(pieces.size());
  unsigned int k = 0;
  for (vector<string>::iterator i = pieces.begin(); i != pieces.end(); i++) {
    char *end;
    x[k++] = strtol(i->c_str(), &end, 10);
    if (i->c_str() == end) throw SimError("error in comma-separated list");
  }
}

/* same as above but for char* */
void valint_from_string(const char *s, valarray<int> &x) {
  string tmp(s);
  valint_from_string(tmp, x);
}

/* END */



