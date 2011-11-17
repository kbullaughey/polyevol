#include "error_handling.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include <sstream>

using std::string;
using std::stringstream;

SimError::SimError(const char *s) {
   detail = s;
}

SimError::SimError(const string &s) {
   detail = s;
}

SimError::SimError(const stringstream &s) {
  detail = s.str();
}

#define SIM_ERROR_BUF 1024
SimError::SimError(int unused, const char *fmt, ...) {
  char *p;
  va_list ap;
  if ((p = (char*)malloc(SIM_ERROR_BUF)) == NULL)
    throw SimError("malloc failed");
  va_start(ap, fmt);
  vsnprintf(p, SIM_ERROR_BUF, fmt, ap);
  va_end(ap);
  detail = string(p);
}

SimUsageError::SimUsageError(const char *s) 
      : SimError(s) {
}

SimUsageError::SimUsageError(const string &s) 
      : SimError(s) {
}

SimUsageError::SimUsageError(const stringstream &s) 
      : SimError(s.str()) {
}


/* END */
