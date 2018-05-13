#ifndef __COMPUTE_MEM_PROB_HEADER
#define __COMPUTE_MEM_PROB_HEADER

// Visible functions.

void   clean_mem_prob (void);
double compute_mem_prob(size_t, size_t);
int    get_mem_prob_error_code (void);
void   reset_mem_prob_error (void);
int    set_params_mem_prob (size_t, size_t, double, double);
void   set_mem_prob_max_precision_on (void);
void   set_mem_prob_max_precision_off (void);

#endif
