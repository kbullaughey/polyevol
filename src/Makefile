CC = g++
HEADERS = command_line.h error_handling.h sim_rand.h common.h genome.h population.h site.h statistic.h running_mean.h
OBJS = quant.o command_line.o error_handling.o sim_rand.o common.o genome.o population.o site.o statistic.o running_mean.o
CFLAGS = -Wall
LIBS = -lm
PLATFORM := $(shell uname -s)

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

all: quant qapprox gtest

quant: $(OBJS) $(HEADERS)
	$(CC) $(CFLAGS) $(INC) $(LIBS) -o $@ $(OBJS)

qapprox: qapprox.c
	gcc -o qapprox qapprox.c  -lm $(GSLLIBS)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@ 

gtest:
	$(MAKE) -C vendor/gtest

clean:
	-rm *.o quant

# END
