#include <stdio.h>
#include "sesame.h"

int main(void) {

   // Compute off-target probabilities for MEM seeds, with
   // seeds of minimum size gamma = 17, on reads of size up
   // to 50 and with sequencing error rate 0.01.
   int success = sesame_set_static_params(17, 50, .01);
   if (!success) return -1;

   double u[] = {.01, .02, .03, .04, .05, .10, .20};
   size_t N[] = {0, 1, 2, 3, 4, 5, 10, 20, 100, 500, 10000};

   // Compute the probabilities.
   for (int i = 0 ; i < 7 ; i++) {
      // Use 'mem_seed_offp' for N up to 20.
      for (int j = 0 ; j < 8 ; j++) {
         double *prob = mem_seed_offp(u[i], N[j]);
         // Store the results in 'H1'.
         if (!store_prob(0, u[i], N[j], prob)) return -1;
      }
      // Above 20, switch to 'mem_seed_offp_mcmc'.
      for (int j = 8 ; j < 11 ; j++) {
         double *prob = mem_seed_offp_mcmc(u[i], N[j]);
         // Store the results in 'H1'.
         if (!store_prob(0, u[i], N[j], prob)) return -1;
      }
   }

   // Save probabilities to file.
   FILE *f = fopen("mem_offp_17_50_001.txt", "w");
   if (f == NULL) return -1;

   dump_prob_to_file(f);
   fclose(f);

   // Retrieve probablities from files.
   FILE *g = fopen("mem_offp_17_50_001.txt", "r");
   if (g == NULL) return -1;

   // This calls 'sesame_set_static_params' again.
   load_prob_from_file(g);
   fclose(g);

   double *prob = fetch_prob(0, .01, 5);
   fprintf(stderr, "u=0.01, N=5, prob:%.8f\n", prob[20]);

   sesame_clean();

   return 0;

}
