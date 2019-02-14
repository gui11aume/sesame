#include <stdio.h>
#include "memseedp.h"

int main(void) {

   // Worst mu: 0.016
   int success = memseedp_set_static_params(17, 100, .01);
   if (!success) return -1;

   // N curves.
   size_t N[] = {0,1,5,20,100,500,10000};

   for (int i = 2 ; i <= 2 ; i++) {
      for (int k = 17 ; k <= 100 ; k++) {
         fprintf(stdout, "\t%.12f",
               memseedp_auto_false_positive(k, .01, N[i]));
      }
      fprintf(stdout, "\n");
   }

#if 0
   int success = set_params_mem_prob(20, 50, .08, .05);
   if (!success) return -1;

   double x = mem_seed_prob(50,0);

   fprintf(stdout, "%.6f %.12f %.12f %.12f %.12f\n",
         0.0, x, x, x, x);

   for (double mu = 0.00125 ; mu <= .10 ; mu += .00125) {

      int success = set_params_mem_prob(20, 50, .08, mu);
      if (!success) return -1;

      double p1 = mem_seed_prob(50, 1);
      double p2 = mem_seed_prob(50, 5);
      double p3 = mem_seed_prob(50, 20);
      double p4 = mem_seed_prob(50, 100);

      fprintf(stdout, "%.4f %.12f %.12f %.12f %.12f\n",
            mu, p1, p2, p3, p4);

      memseedprob_clean(); 
   }
#endif

   return 0;

}
