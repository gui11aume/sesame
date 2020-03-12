vpath %.c ..
vpath %.c lib
vpath %.h ..
vpath %.h lib

RUN=run_tests_for_coverage
DBG=run_tests_for_debug
VAL=run_tests_for_valgrind

OBJECTS= libunittest.so unittests_sesame.o
SOURCES= sesame.c

ifdef TRAVIS_COMPILER
CC= $(TRAVIS_COMPILER)
else
CC= gcc
endif
INCLUDES= -I.. -Ilib
CFLAGS= -std=gnu99 -g -Wall -Wextra $(INCLUDES)
COVERAGE= -fprofile-arcs -ftest-coverage
PROFILE= -pg
LDLIBS= -L. -Wl,-rpath,. -lunittest -lz -lm -lpthread
# Use different flags on Linux and MacOS.
ifeq ($(shell uname -s),Darwin)
	libflag= -dynamiclib
else
	libflag= -shared
endif

$(RUN): CFLAGS += -O2 $(COVERAGE) $(PROFILE)
$(RUN): $(OBJECTS) $(SOURCES) $(OBJECTS)
	$(CC) $(CFLAGS) runtests.c $(OBJECTS) $(LDLIBS) -o $@

$(DBG): CFLAGS += -O0 -DDEBUG
$(DBG): $(OBJECTS) $(SOURCES) $(OBJECTS)
	$(CC) $(CFLAGS) runtests.c $(OBJECTS) $(LDLIBS) -o $@

$(VAL): CFLAGS += -O2 -DVALGRIND
$(VAL): $(OBJECTS) $(SOURCES) $(OBJECTS)
	$(CC) $(CFLAGS) runtests.c $(OBJECTS) $(LDLIBS) -o $@

clean:
	rm -rf $(RUN) $(DBG) $(VAL) $(OBJECTS) \
		*.gcda *.gcno *.gcov gmon.out \
		.inspect.gdb libunittest.so.dSYM runtests.dSYM

libunittest.so: unittest.c
	$(CC) -fPIC $(libflag) $(CFLAGS) -o libunittest.so lib/unittest.c

test: $(RUN)
	./$(RUN)

inspect: $(DBG)
	gdb --command=.inspect.gdb --args $(DBG)

debug: $(DBG)
	gdb --command=.debug.gdb --args $(DBG)

valgrind: $(VAL)
	valgrind --leak-check=full --show-leak-kinds=all ./$(VAL)

vgdb: $(VAL)
	valgrind --vgdb=yes --vgdb-error=0 ./$(VAL)

unittests_sesamerob.gcda gmon.out: $(RUN)
	-./$(RUN)

cov: unittests_sesame.gcda
	gcov unittests_sesame.c
	cat sesame.c.gcov

prof: gmon.out
	gprof ./$(RUN) --flat-profile --brief
