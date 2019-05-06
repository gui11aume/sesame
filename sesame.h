#ifndef _SESAME_DECLARED_ 
#define _SESAME_DECLARED_ 
// Initialization and clean up.
int sesame_set_static_params(size_t, size_t, double);
void sesame_clean(void);

// Options.
void sesame_set_method(int);
void sesame_set_max_prcsn(void);
void sesame_unset_max_prcsn(void);

// MEM seeding probabilities.
double sesame_auto_false_positive(int, double, int);
#endif
