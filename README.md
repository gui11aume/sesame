## Sesame: Computing seeding probabilities ##
[![Build Status](https://travis-ci.org/gui11aume/sesame.svg?branch=master)](https://travis-ci.org/gui11aume/sesame)
---
## Contents: ##
    1. What is sesame?
    2. Compilation
    3. Using sesame
    4. License
    5. Citation

---

## I. What is sesame?   ##

Sesame is an open-source C library to compute seeding probabilities in
short DNA reads.

The most common method to identify a sequencing read is to align it to a
reference sequence, typically a genome. The alignment aims to find the
subsequence of the genome that is the closest to the read. There exist
algorithms to find the best subsequence of the genome, but they are too
slow when the reference sequence is the size of a genome.

There are ways to speed up the alignment using heuristics. They provide a
huge speed boost at the cost of a small chance of missing the best
subsequence. One of the most popular heuristics is to use _seeds_, i.e.
short matches between the read and the genome, to create a shortlist of
candidate subsequences. Checking only the shortlisted subsequences speeds
up the process, but if the best subsequence was not shortlisted then it is
missed entirely.

Sesame computes the probabilities of mapping errors when using different
seeding strategies. These probabilities depend on the type of seed
(sesame supports exact seeds, skip seeds and MEM seeds), and whether the
read was taken from a repeated sequence.

Knowing the probability that the seeding heuristic fails allows the user to:

  1. design a rational seeding strategy,
  1. optimize the speed/sensitivity trade-off,
  1. give accurate confidence score for the reads.


II. Compilation
---------------

To use sesame, clone this git repository.

    git clone https://github.com/gui11aume/sesame

The files should be downloaded in a folder named `sesame`. The library
consists of a header file `sesame.h` and a source file `sesame.c`.
Sesame does not depend on any external library. It works on Linux and
on Mac.

The repository contains an example program using sesame. Check that
you can compile sesame on your machine by running `make` (Mac users
require 'xcode', available at the Mac Appstore).

    cd sesame
    make

You can also run the test cases. If any of them fails, it means that
sesame is not running as expected in your environment (in this case
you can open an issue).

    cd unittests
    make test


III. Using sesame
-----------------

Sesame runs on Linux and Mac. It has not been tested on Windows.
It is a C library that you can include in your own C or C++ projects.
Copy `sesame.h` and `sesame.c` in the appropriate directories of
your project. Import the sesame functions as usual in all the source
files that must have access to them.

```C
#include "sesame.h"
```

### Initialization and cleanup

Before calling any function of sesame, you must first initialize
the _static parameters_. Those are the parameters that are constant
for a given sequencing run. They are the minimum seed size, the
maximum size of the reads, and the error rate of the sequencer.
For instance, to initialize the library for seeds of size 17 or
greater in reads of size 100 from a sequencer with 1% error, you
would use the following:

```C
int success = sesame_set_static_params(17, 100, .01);
```

What to do in case the library cannot be initialized is up to you
(for instance this can happen if you run out of memory). This means
that you will not be able to use sesame at run time, so most of the
time you will just bail out.

```C
if (!success) {
   exit(EXIT_FAILURE);
}
```

When you no longer need sesame, it is recommended that you free the
memory used internally.

```C
sesame_clean();
```

Below are the prototypes of the two functions:

```C
void sesame_clean(void);
int sesame_set_static_params(size_t, size_t, double);
```

### Automatic functions

Sesame has six user-friendly functions that can be used to compute
seeding probabilities without having to worry about the detail.
The functions return the probability of interest as a single number
of type `double`. But internally, sesame computes the probabilities
for multiple read sizes and stores the results for rapid recall.
It also maps the input parameters to a coarse-grain grid to speed
up the recall, so it is important to understand that the results
are not exactly the ones that are requested (for this, see the
functions described in the next section). Those functions are thus
for users who want a quick and approximate result with relatively
little overhead in their own code.

Sesame computes the probabilties of errors for exact seeds, skip
seeds and MEM seeds. Briefly, exact seeds are perfect matches of
a given size, skip seeds are the same except that they are disallowed
to start on a certain number of consecutive positions, and MEM seeds
(Maximal Exact Match) seeds are exact matches that cannot be extended
left or right (see reference article for explanations). Sesame also
distinguishes the probability that the seeding step returns no
candidate (null seeding), or that it returns only wrong candidates
(off-terget seeding).

The arguments that are required by all the functions are the size
of the sequence where the seed is sought, the divergence rate
between the duplicates of the sequence, and the number of said
duplicates.

#### Off-target probability

For instance, to compute the probability that **MEM seeding** will be
off target for a read of `size 100` from a sequence with `50 extra
duplicates` in the genome that typically `differ at 6%` of the positions
is computed by

```C
double p = auto_mem_seed_offp (100, .06, 50);
```
    
#### Null seeding probability

To compute the probability that **skip-8 seeding** will be null (no true
candidate is found) for a read of `size 100` from a sequence with `50 extra`
duplicates in the genome that typically `differ at 6%` of the positions is 
computed by

```C
double p = auto_skip_seed_nullp (100, 8, .06, 50);
```
    
Note that when using skip seeds, one extra argument is the amount of 
skipping (the number of consecutive positions where a new seed is 
disallowed).

#### Function prototypes

Below are the protypes of the six automatic functions.
```C
// Exact seeds
double auto_exact_seed_nullp (int size, double perror, int duplicates);
double auto_exact_seed_offp (int size, double perror, int duplicates);
// Skip seeds
double auto_skip_seed_nullp (int size, int skip, double perror, int duplicates);
double auto_skip_seed_offp (int size, int skip, double perror, int duplicates);
// MEM seeds
double auto_mem_seed_nullp (int size, double perror, int duplicates);
double auto_mem_seed_offp (int size, double perror, int duplicates);
```


### Non-automatic functions

For the users who want to have complete control over the results,
sesame has six functions that are similar to automatic functions,
except that the parameters are not coarse grainded and that the
results are not stored.

Importantly, those functions do not return a single number but an
array of probabilities for all read sizes up to the secified
maximum (see `sesame_set_static_params`). This depends on internal
calls to `malloc` and it is the responsability of the caller to
free the requeted memory. Sesame offers an option to store the
results for recall and to dump them to a text file (see next
sections).

Non-automatic functions have one parater fewer than automatic
functions because the read size is not required (the probabilities
for all the sizes up to the maximum are computed).

Below are the protypes of the six non-automatic functions.

```C
double * exact_seed_null (double perror, int duplicates);
double * exact_seed_offp (double perror, int duplicates);
double * skip_seed_null (int skip, double perror, int duplicates);
double * skip_seed_offp (int skip, double perror, int duplicates);
double * mem_seed_null (double perror, int duplicates);
double * mem_seed_offp (double perror, int duplicates);
double * mem_seed_offp_mcmc (double perror, int duplicates);
```

### Storing results in memory

To facilitate the use of non-automatic functions, sesame allows the
user to store the results of the computations in an internal hash for
fast recall. The keys of the hash have three parameters. A suggestion
is to use the amount of skipping, the divergence rate and the
number of duplicates, but they can be set to any value specified by the
user.

The internal hash does not make a copy of the results, it only stores
a pointer to them. Storing results with the same key will overwrite the
associated pointer, which can cause memory leaks. Freeing the memory
is the responsability of the user if a hash entry is overwritten.

Below is an example of commands to compute the probability of **null
seeding** for **skip-8 seeds** and sequences with `50 extra duplicates` that
typically `differ from the target by 6%`, and then to store the results
for later recall:

```C
double * p = skip_seed_nullp(8, .06, 50);
int success = store_prob(8, .06, 50, p);
```

The probabilities can then be recalled with

```C
double * q = fetch_prob(8, .06, 50);
```

If no record has been previously stored in the requested key, `fetch_prob` 
returns a `NULL` pointer. Below are the prototypes of the functions to store 
and retrieve arrays of probabilities:

```C
double * fetch_prob (int, double, int);
int store_prob (int, double, int, double *);
```

To erase the entries of the internal hash and to free the computed
probabilities, use the following function:

```C
clean_prob_store();
```

This operation is automatically performed upon cleanup by
`sesame_clean()`.


### Saving and loading from disk

Sesame allows the users to dump the internal hash to a text file.
The arrays of probabilities stored in the internal hash are then
converted to a simple tab-separated text format. The header of the
file contains the associated values of the static parameters.

To give an example, storing the computed results in a file called
`seed_probs.txt` can be done with the following commands:

```C
FILE * f = fopen("seed_probs.txt", "w");
dump_prob_to_file(f);
```

After this, the file `seed_probs.txt` can be loaded in future
runs with the commands:

```C
FILE * f = fopen("seed_probs.txt", "r");
load_prob_from_file(f);
```

The internal hash is erased before loading a file from disk, which
will also free any stored array of probabilities. Loading a file
from disk will also reset the static parameters and erase the
precomputations of the automatic functions. For these reasons,
loading a file from disk should be done before calling any other
sesame functions, even `sesame_set_static_params` (which is no
longer required after loading the file).

Below are the prototypes for the input/ouput functions.

```C
void dump_prob_to_file (FILE *);
void load_prob_from_file (FILE *);
```

IV. License
-----------

Sesame is open (of course it is). It is licensed under the GNU General
Public License, version 3 (GPLv3), for more information read the LICENSE
file or refer to:

  http://www.gnu.org/licenses/


V. Citation
--------------

Sesame is in the process of peer review. For now you can cite it as a bioRxiv preprint.

Filion GJ, Cortini R, Zorita E 2019.
[Calibrating seed-based heuristics to map short DNA reads](https://www.biorxiv.org/content/10.1101/619155v1.full)
bioRxiv DOI:10.1101/619155.
