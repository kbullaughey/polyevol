CC = g++
HEADERS = command_line.h error_handling.h sim_rand.h common.h genome.h population.h site.h statistic.h running_mean.h
OBJS = quant.o command_line.o error_handling.o sim_rand.o common.o genome.o population.o site.o statistic.o running_mean.o
CFLAGS = -Wall
LIBS = -lm
PLATFORM := $(shell uname -s)
ROOT := $(shell pwd)
TEST_SUPPORT = test/support
TESTS = $(shell find test -name '*_test.cpp')
TEST_OBJECTS = $(TESTS:.cpp=.o)

# I use Google Test as my testing framework, which is in the repository at src/vendor/gtest
GTEST_DIR = vendor/gtest

ifeq ($(profile),yes)
  CFLAGS += -pg -O3
else
  ifeq ($(debug),yes)
    CFLAGS += -g
  else
    CFLAGS += -O3  
  endif
endif

ifeq ($(PLATFORM),Linux)
  GSLLIBS = -lgsl -lgslcblas
else ifeq ($(PLATFORM),Darwin)
  GSLLIBS = -lgsl
else
  $(error Error: Unsupported platform)
endif

all: quant qapprox $(TEST_SUPPORT)/libgtest.a test/runner

quant: $(OBJS) $(HEADERS)
	$(CC) $(CFLAGS) $(INC) $(LIBS) -o $@ $(OBJS)

qapprox: qapprox.c
	gcc -o qapprox qapprox.c -lm $(GSLLIBS)

test/%.o: test/%.cpp
	$(CC) -I. -Ivendor/gtest/include $(CFLAGS) -c $< -o $@

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@ 

$(TEST_SUPPORT)/libquant.a: $(OBJS) $(HEADERS)
	ar -rvs $(TEST_SUPPORT)/libquant.a $(OBJS)

# If the gtest Makefile doesn't exist, assume we need to run configure
$(GTEST_REL)/Makefile: 
	(cd vendor/gtest && ./configure)

# Make the gtest-all.o library
$(TEST_SUPPORT)/gtest-all.o: $(GTEST_DIR)/Makefile
	$(MAKE) -C vendor/gtest
	$(CC) -I$(GTEST_DIR)/include -I$(GTEST_DIR) -c $(GTEST_DIR)/src/gtest-all.cc -o $(TEST_SUPPORT)/gtest-all.o

# Make the gtest static library
$(TEST_SUPPORT)/libgtest.a: $(TEST_SUPPORT)/gtest-all.o
	ar -rvs $(TEST_SUPPORT)/libgtest.a $(TEST_SUPPORT)/gtest-all.o

# Make the tests
test/runner: test/runner.o $(TEST_OBJECTS) $(TEST_SUPPORT)/libquant.a
	$(CC) -o test/runner $(TEST_OBJECTS) test/runner.o -Ltest/support -lgtest -lquant

blah:
	echo $(TESTS)
	echo $(TEST_OBJECTS)

clean:
	-rm *.o quant
	-rm libgtest.a gtest-all.o
	$(MAKE) -C vendor/gtest clean

# END
