P= run_example

OBJECTS= sesame.o

ifdef TRAVIS_COMPILER
CC= $(TRAVIS_COMPILER)
else
CC= gcc
endif
CFLAGS= -std=gnu99 -Wall -Wextra
LDLIBS= -lm

all: CFLAGS += -O3
all: $(P)

debug: CFLAGS += -DDEBUG -g3 -O0
debug: $(P)

$(P): example_sesame_dump.c $(OBJECTS)
	gcc $(CFLAGS) $^ -lm -o $@

analyze: CC= clang --analyze
analyze: CFLAGS += -DDEBUG -g -O0
analyze: $(OBJECTS)

tidy:
	clang-tidy sesame.c --

warnings: CFLAGS += -g -O0 -Wunused-parameter -Wredundant-decls \
        -Wreturn-type -Wswitch-default -Wunused-value -Wimplicit \
        -Wimplicit-function-declaration -Wimplicit-int -Wimport \
        -Wunused  -Wunused-function -Wunused-label -Wbad-function-cast \
        -Wno-int-to-pointer-cast -Wpointer-sign -Wnested-externs \
        -Wold-style-definition -Wstrict-prototypes -Wredundant-decls \
        -Wunused -Wunused-function -Wunused-parameter -Wunused-value \
        -Wformat -Wunused-variable -Wformat-nonliteral -Wparentheses \
        -Wundef -Wsequence-point -Wuninitialized -Wbad-function-cast \
	-Wall -Wextra
warnings: $(OBJECTS)

clean:
	rm -f $(P) $(OBJECTS)
