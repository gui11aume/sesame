#ifndef _MEM_SEED_PROB_DECLARED_
#define _MEM_SEED_PROB_DECLARED_

// Visible functions.

void   clean_mem_prob (void);
double mem_seed_prob(size_t, size_t);
int    get_mem_prob_error_code (void);
void   reset_mem_prob_error (void);
int    set_params_mem_prob (size_t, size_t, double, double);

#endif
