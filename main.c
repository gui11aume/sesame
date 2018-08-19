#include <stdio.h>
#include "mem_seed_prob.h"

int main(void) {

   int success = set_params_mem_prob(17, 100, .01, .016);
   if (!success) return -1;

   size_t N[] = {0,1,5,20,100,500,10000};

   for (int i = 0 ; i <= 6 ; i++) {
      fprintf(stdout, "%.12f", mem_seed_prob(17, N[i]));
      for (int k = 18 ; k <= 100 ; k++) {
         fprintf(stdout, "\t%.12f", mem_seed_prob(k, N[i]));
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

      clean_mem_prob(); 
   }
#endif

   return 0;

}
