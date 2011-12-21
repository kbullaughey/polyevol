CC = g++
HEADERS = command_line.h error_handling.h sim_rand.h common.h genome.h population.h site.h
OBJS = quant.o command_line.o error_handling.o sim_rand.o common.o genome.o population.o site.o
CFLAGS = -Wall
LIBS = -lm

ifeq ($(profile),yes)
CFLAGS += -pg -O3
else
ifeq ($(debug),yes)
CFLAGS += -g
else
CFLAGS += -O3  
endif
endif

all: quant qapprox

quant: $(OBJS) $(HEADERS)
	$(CC) $(CFLAGS) $(INC) $(LIBS) -o $@ $(OBJS)

qapprox: qapprox.c
	gcc -o qapprox qapprox.c  -lm -lgsl

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@ 

clean:
	-rm *.o quant

# END
