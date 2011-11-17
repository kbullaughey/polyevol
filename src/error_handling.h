#ifndef __ERROR_HANDLING_H__
#define __ERROR_HANDLING_H__

#include <string>
#include <sstream>

class SimError {
public:
  SimError(const char *s = "unknown error");
  SimError(const std::string &s);
  SimError(const std::stringstream &s);
  SimError(int unused, const char *fmt, ...);
  std::string detail;
};

class SimUsageError : public SimError {
public:
  SimUsageError(const char *s);
  SimUsageError(const std::string &s);
  SimUsageError(const std::stringstream &s);
};

#endif /* __ERROR_HANDLING_H__ */
