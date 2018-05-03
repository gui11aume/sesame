vpath %.c ..
vpath %.c lib
vpath %.h ..
vpath %.h lib

P= runtests

OBJECTS= libunittest.so unittests_compute_mem_prob.o
SOURCES= compute_mem_prob.c

CC= gcc
INCLUDES= -I.. -Ilib
COVERAGE= -fprofile-arcs -ftest-coverage
PROFILE= -pg
CFLAGS= -std=gnu99 -g -Wall -O0 $(INCLUDES) $(COVERAGE) $(PROFILE)
LDLIBS= -L. -lunittest -lz -lm -lpthread
# Use different flags on Linux and MacOS.
ifeq ($(shell uname -s),Darwin)
	libflag= -dynamiclib
else
	libflag= -shared
endif

$(P): $(OBJECTS) $(SOURCES) runtests.c
	$(CC) $(CFLAGS) runtests.c $(OBJECTS) $(LDLIBS) -o $@

clean:
	rm -f $(P) $(OBJECTS) *.gcda *.gcno *.gcov gmon.out .inspect.gdb

libunittest.so: unittest.c
	$(CC) -fPIC $(libflag) $(CFLAGS) -o libunittest.so lib/unittest.c

test: $(P)
	./$(P)

inspect: $(P)
	gdb --command=.inspect.gdb --args $(P)

debug: $(P)
	gdb --command=.debug.gdb --args $(P)

valgrind: $(P)
	#valgrind --vgdb=yes --vgdb-error=0 ./$(P)
	valgrind --leak-check=full --show-leak-kinds=all ./$(P)

unittests_compute_mem_prob.gcda gmon.out: $(P)
	-./$(P)

cov: unittests_compute_mem_prob.gcda
	gcov unittests_compute_mem_prob.c
	cat compute_mem_prob.c.gcov

prof: gmon.out
	gprof ./$(P) --flat-profile --brief
