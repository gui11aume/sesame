P= main

OBJECTS= mem_seed_prob.o random.o

CC= gcc
CFLAGS= -std=gnu99 -Wall
LDLIBS= -lm

all: CFLAGS += -O3
all: $(P)

debug: CFLAGS += -DDEBUG -g -O0
debug: $(P)

profile: CFLAGS += -g -pg -O3 -fprofile-arcs -ftest-coverage
profile: $(P)


$(P): main.c $(OBJECTS)
	gcc $(CFLAGS) $^ -lm -o $@

clean:
	rm -f $(P) $(OBJECTS)
