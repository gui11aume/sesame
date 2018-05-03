#include <stdio.h>
#include "compute_mem_prob.h"

int main(void) {

   // First fix parameters to compute seeding probabilities.
   // In this example, the parameters are
   // 1. minimum MEM seed size: 17
   // 2. maximum read size: 50
   // 3. sequencing error rate: 0.01
   // 4. divergence rate: 0.05
   int success = set_params_mem_prob(17, 50, 0.01, 0.05);

   // Success of initialization can be checked
   // from the return value of the function.
   if (!success) return -1;

   double x = compute_mem_prob(2, 18);

   // If something went wrong, the function returns 'nan'.
   // This can be checked with the condition 'x != x'.
   if (x != x) return -1;

   // Alternatively, one can directly access the last error.
   if (get_mem_prob_error_code()) return -1;

   // After using the library, one needs to call 'clean_mem_prob()'
   // to free internal variables and reset the internal state.
   clean_mem_prob(); 

   fprintf(stdout, "%.5f\n", x);
   return 0;

}
