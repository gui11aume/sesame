#ifndef _MEMSEEDP_DECLARED_
#define _MEMSEEDP_DECLARED_
// Initialization and clean up.
int memseedp_set_static_params (size_t, size_t, double);
void memseedp_clean (void);

// Options.
void memseedp_set_method (int);
void memseedp_set_max_prcsn (void);
void memseedp_unset_max_prcsn (void);

// MEM seeding probabilities.
double memseedp_auto_false_positive(int, double, int);
#endif
