#include "unittest.h"
#include "mem_seed_prob.c"

void
test_set_params_mem_prob
(void)
{

   int success;

   success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   test_assert(G == 17);
   test_assert(K == 50);
   test_assert(P == 0.01);
   test_assert(U == 0.05);

   test_assert(KSZ == sizeof(trunc_pol_t) + 51*sizeof(double));

   for (int i = 0 ; i < MAXN ; i++) {
      test_assert(ARRAY[i] == NULL);
   }

   clean_mem_prob();

   return;

}


void
test_error_set_params_mem_prob
(void)
{

   int success;

   // Case 1.
   redirect_stderr();
   // The error is that 'p' is 0.0.
   success = set_params_mem_prob(17, 50, 0.00, 0.05);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   // Case 2.
   redirect_stderr();
   // The error is that 'mu' is 0.0.
   success = set_params_mem_prob(17, 50, 0.01, 0.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   // Case 3.
   redirect_stderr();
   // The error is that 'p' is 1.0.
   success = set_params_mem_prob(17, 50, 0.01, 0.00);
   success = set_params_mem_prob(17, 50, 1.00, 0.05);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   // Case 4.
   redirect_stderr();
   // The error is that 'mu' is 1.0.
   success = set_params_mem_prob(17, 50, 0.01, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   // Case 5.
   redirect_stderr();
   // The error is that 'G' (gamma) is 0.
   success = set_params_mem_prob(0, 50, 0.01, 0.05);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   // Case 6.
   redirect_stderr();
   // The error is that 'K' is 0.
   success = set_params_mem_prob(17, 0, 0.01, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   // Casae 7. Test memory error.
   set_alloc_failure_rate_to(1.0);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   success = set_params_mem_prob(17, 50, 0.01, 0.05);
   unredirect_stderr();
   reset_alloc();
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   return;

}


void
test_uninitialized_error
(void)
{

   // Do not call 'set_params_mem_prob()'.
   redirect_stderr();
   double x = mem_seed_prob(5, 20);
   unredirect_stderr();
   test_assert_stderr("[mem_seed_prob] error in function `mem_");
   test_assert(x != x);

   return;

}


void
test_omega
(void)
{

   // Only used to set global 'U'.
   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Simple tests for small values.
   test_assert(fabs(omega(0,0)-1.0) < 1e-9);
   test_assert(fabs(omega(1,0)-0.0) < 1e-9);
   test_assert(fabs(omega(2,0)-0.0) < 1e-9);
   test_assert(fabs(omega(0,1)-(1-.05/3)) < 1e-9);
   test_assert(fabs(omega(1,1)-.05/3) < 1e-9);
   test_assert(fabs(omega(2,1)-0.0) < 1e-9);
   test_assert(fabs(omega(0,2)-(1-.05/3)*(1-.05/3)) < 1e-9);
   test_assert(fabs(omega(1,2)-2*(1-.05/3)*.05/3) < 1e-9);
   test_assert(fabs(omega(2,2)-.05/3*.05/3) < 1e-9);

   // Tests for high values checked with R.
   test_assert(fabs(omega(0,9)-0.8596206731) < 1e-9);
   test_assert(fabs(omega(1,9)-0.1311285773) < 1e-9);
   test_assert(fabs(omega(5,9)-1.515016e-07) < 1e-9);
   test_assert(fabs(omega(0,23)-0.6793874326) < 1e-9);
   test_assert(fabs(omega(1,23)-0.2648459483) < 1e-9);
   test_assert(fabs(omega(5,23)-3.197640e-05) < 1e-9);
   test_assert(fabs(omega(0,58)-0.3772629470) < 1e-9);
   test_assert(fabs(omega(1,58)-0.3708686598) < 1e-9);
   test_assert(fabs(omega(5,58)-0.0024179659) < 1e-9);

   // Tests for small values checked with R.
   // Here we use relative error because the values are small.
   test_assert(fabs(omega(9,9)/9.922903012752172135e-17 - 1) < 1e-6);
   test_assert(fabs(omega(18,23)/3.046165264390310362e-28 - 1) < 1e-6);
   test_assert(fabs(omega(23,23)/1.266255198051510610e-41 - 1) < 1e-6);
   test_assert(fabs(omega(18,58)/2.262004384436928424e-18 - 1) < 1e-6);
   test_assert(fabs(omega(34,58)/2.992106656854436556e-45 - 1) < 1e-6);
   test_assert(fabs(omega(44,58)/4.627213397203421023e-66 - 1) < 1e-6);
   test_assert(fabs(omega(58,58)/7.365928141161108488e-104 - 1) < 1e-6);

   clean_mem_prob();

}


void
test_psi
(void)
{

   // Only used to set global 'U'.
   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   test_assert(fabs(psi(1,0,0,0,1)-1.0) < 1e-9);
   test_assert(fabs(psi(17,0,1,1,1)-1.0) < 1e-9);
   test_assert(fabs(psi(15,1,3,7,12)-0.0) < 1e-9);
   test_assert(fabs(psi(15,7,1,3,12)-0.0) < 1e-9);
   test_assert(fabs(psi(30,14,19,15,23)-0.0) < 1e-9);

   // Tests for values checked with R.
   test_assert(fabs(psi(17,1,3,2,5)-2.1637593296) < 1e-9);
   test_assert(fabs(psi(15,1,7,3,12)-14.466081668) < 1e-9);

}


void
test_zeta
(void)
{

   // Only used to set global 'U'.
   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Simple test for small values.
   test_assert(fabs(zeta(20,0,0,1)-0.0) < 1e-9);
   test_assert(fabs(zeta(17,0,1,1)-pow(.95,17)) < 1e-9);

   // Tests for values checked with R.
   test_assert(fabs(zeta(20,0,1,1)-0.35848592241) < 1e-9);
   test_assert(fabs(zeta(20,4,5,8)-0.65255894306) < 1e-9);
   test_assert(fabs(zeta(20,6,3,8)-0.25509528217) < 1e-9);
   test_assert(fabs(zeta(30,14,19,23)-0.83001525276) < 1e-9);

   clean_mem_prob();

}


void
test_new_zero_trunc_pol
(void)
{

   size_t ksz = 50;

   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *a = new_zero_trunc_pol();
   test_assert_critical(a);

   test_assert(a->monodeg == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   free(a);
   clean_mem_prob();

   return;

}


void
test_trunc_pol_mult
(void)
{
   
   // Test multiplications between zero polynomials.

   size_t ksz = 50;

   set_params_mem_prob(17, ksz, 0.01, 0.05);
   trunc_pol_t *a = new_zero_trunc_pol();

   test_assert_critical(a != NULL);

   // Note: 'trunc_pol_mult' returns NULL but 'a' is not set to NULL.
   test_assert(trunc_pol_mult(a, NULL, NULL) == NULL);

   test_assert_critical(a != NULL);
   test_assert(a->monodeg == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   trunc_pol_t *b = new_zero_trunc_pol();
   test_assert_critical(b != NULL);

   test_assert(trunc_pol_mult(a, b, NULL) == NULL);

   // Same remark as above.
   test_assert_critical(a != NULL);
   test_assert(a->monodeg == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   test_assert(trunc_pol_mult(a, NULL, b) == NULL);

   // Same remark as above.
   test_assert_critical(a != NULL);
   test_assert(a->monodeg == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   trunc_pol_t *c = new_zero_trunc_pol();
   test_assert_critical(c != NULL);

   // Here 'b' and 'c' are still zero polynomials,
   // so same remark as above.
   test_assert(trunc_pol_mult(a, b, c) == NULL);

   // Same remark as above.
   test_assert_critical(a != NULL);
   test_assert(a->monodeg == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }


   // Test multiplications between monomials (b = 5z and c = z^2).
   b->monodeg = 1; b->coeff[1] = 5;
   c->monodeg = 2; c->coeff[2] = 1;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->monodeg == 3);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = i == 3 ? 5 : 0;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->monodeg == 3);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = i == 3 ? 5 : 0;
      test_assert(a->coeff[i] == target);
   }

   // Test overflow (multiplying two monomials of sufficient
   // degree yields a zero polynomial).
   bzero(b, KSZ);
   bzero(c, KSZ);
   b->monodeg = ksz-1; b->coeff[1] = 0; b->coeff[ksz-1] = 5;
   c->monodeg = ksz-1; c->coeff[2] = 0; c->coeff[ksz-1] = 1;

   test_assert(trunc_pol_mult(a, b, c) == NULL);
   test_assert(a != NULL);


   // Test multiplications between a monomial and a
   // polynomial (b = 5z and c = z^2 + 2z^3).
   bzero(b, KSZ);
   bzero(c, KSZ);
   b->monodeg = 1; b->coeff[1] = 5;
   c->monodeg = ksz+1; c->coeff[2] = 1; c->coeff[3] = 2;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->monodeg > ksz);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 10;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->monodeg > ksz);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 10;
      test_assert(a->coeff[i] == target);
   }

   // Test multiplications between two polynomials
   // (b = 5z + 3z^2 and c = z^2 + 2z^3).
   bzero(b, KSZ);
   bzero(c, KSZ);
   b->monodeg = ksz+1; b->coeff[1] = 5; b->coeff[2] = 3;
   c->monodeg = ksz+1; c->coeff[2] = 1; c->coeff[3] = 2;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->monodeg > ksz);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 13;
      if (i == 5) target = 6;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->monodeg > ksz);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 13;
      if (i == 5) target = 6;
      test_assert(a->coeff[i] == target);
   }

   // Test that higher terms overflow.
   // (b = 5z + 3z^2 + z^49 and c = z^2 + 2z^3 + z^49).
   b->coeff[ksz-1] = 1;
   c->coeff[ksz-1] = 1;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->monodeg > ksz);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3)  target = 5;
      if (i == 4)  target = 13;
      if (i == 5)  target = 6;
      if (i == 50) target = 5;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->monodeg > ksz);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3)  target = 5;
      if (i == 4)  target = 13;
      if (i == 5)  target = 6;
      if (i == 50) target = 5;
      test_assert(a->coeff[i] == target);
   }

   free(a);
   free(b);
   free(c);
   clean_mem_prob();

   return;

}


void
test_matrix_mult
(void)
{

   size_t ksz = 50;

   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert_critical(success);

   matrix_t *mat1 = new_zero_matrix(2);
   matrix_t *mat2 = new_zero_matrix(2);
   test_assert_critical(mat1 != NULL);
   test_assert_critical(mat2 != NULL);

   
   //  | 1   z |    |  0   z  |
   //  |z^2 z^3|    |2z^2 3z^3|

   // Fill dummy matrices mat1 and mat2.
   for (int i = 0; i < 4 ; i++) {
      trunc_pol_t *w1 = mat1->term[i];
      trunc_pol_t *w2 = mat2->term[i];
      w1->monodeg = i;
      w1->coeff[w1->monodeg] = 1;
      w2->monodeg = i;
      w2->coeff[w2->monodeg] = i;
   }

   matrix_t *tmp1 = new_zero_matrix(2);
   test_assert_critical(tmp1 != NULL);

   matrix_mult(tmp1, mat1, mat2);

   const double tmp1_array1[51] = {0,0,0,2};
   const double tmp1_array2[51] = {0,1,0,0,3};
   const double tmp1_array3[51] = {0,0,0,0,0,2};
   const double tmp1_array4[51] = {0,0,0,1,0,0,3};

   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tmp1->term[i]->monodeg > ksz);
   }
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(tmp1->term[0]->coeff[i] == tmp1_array1[i]);
      test_assert(tmp1->term[1]->coeff[i] == tmp1_array2[i]);
      test_assert(tmp1->term[2]->coeff[i] == tmp1_array3[i]);
      test_assert(tmp1->term[3]->coeff[i] == tmp1_array4[i]);
   }

   matrix_t *tmp2 = new_zero_matrix(2);
   test_assert_critical(tmp2 != NULL);

   matrix_mult(tmp2, mat2, mat1);

   const double tmp2_array1[51] = {0,0,0,1};
   const double tmp2_array2[51] = {0,0,0,0,1};
   const double tmp2_array3[51] = {0,0,2,0,0,3};
   const double tmp2_array4[51] = {0,0,0,2,0,0,3};

   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tmp2->term[i]->monodeg > ksz);
   }
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(tmp2->term[0]->coeff[i] == tmp2_array1[i]);
      test_assert(tmp2->term[1]->coeff[i] == tmp2_array2[i]);
      test_assert(tmp2->term[2]->coeff[i] == tmp2_array3[i]);
      test_assert(tmp2->term[3]->coeff[i] == tmp2_array4[i]);
   }

   destroy_mat(mat1);
   destroy_mat(mat2);
   destroy_mat(tmp1);
   destroy_mat(tmp2);
   clean_mem_prob();

}


void
test_error_matrix_mult
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   matrix_t *mat1 = new_zero_matrix(2);
   matrix_t *mat2 = new_zero_matrix(3);
   matrix_t *tmp = new_zero_matrix(3);
   test_assert_critical(mat1 != NULL);
   test_assert_critical(mat2 != NULL);
   test_assert_critical(tmp != NULL);

   redirect_stderr();
   // The error is that matrices are not congruent.
   test_assert(matrix_mult(tmp, mat1, mat2) == NULL);
   unredirect_stderr();

   test_assert_stderr("[mem_seed_prob] error in function `matrix_");

   destroy_mat(mat1);
   destroy_mat(mat2);
   destroy_mat(tmp);
   clean_mem_prob();

}


void
test_error_new_zero_trunc_pol
(void)
{
   trunc_pol_t *w;

   redirect_stderr();
   w = new_zero_trunc_pol();
   unredirect_stderr();

   test_assert(w == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   w = new_zero_trunc_pol();
   unredirect_stderr();
   reset_alloc();

   test_assert(w == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

}


void
test_new_trunc_pol_A
(void)
{

   size_t ksz = 50;
   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *A;
   double factor;

   A = new_trunc_pol_A(0,0,1);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg > ksz);
   factor = .01 * (1-.05/3);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= ksz ; i++) {
      test_assert(fabs(A->coeff[i]) < 1e-9);
   }

   free(A);
   A = NULL;

   A = new_trunc_pol_A(1,0,1);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg > ksz);
   factor = .01 * (1-.05/3);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   // The terms below are exactly the same as above. They
   // are separated only to highlight the two terms of the sum.
   for (int i = 18 ; i <= ksz ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   free(A);
   A = NULL;

   A = new_trunc_pol_A(0,1,1);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg > ksz);
   factor = .01 * .05/3;
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= ksz ; i++) {
      double target = pow(.95 * .99, i-1) * .01 * .05/3;
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   free(A);
   A = NULL;

   A = new_trunc_pol_A(1,1,1);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg > ksz);
   factor = .01 * .05/3;
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   // The terms below are exactly the same as above. They
   // are separated only to highlight the two terms of the sum.
   for (int i = 18 ; i <= ksz ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   free(A);
   A = NULL;

   // This one was checked with R.
   A = new_trunc_pol_A(3,5,9);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg > ksz);
   // This is "p times omega sub n".
   factor = .01 * 1.5150164144e-07;
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double xi_term = 1 - pow(1-pow(.95, i-1), 9);
      double target = factor * xi_term * pow(.99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   // The zeta terms computed with R.
   double zeta_terms[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0.821756652886299, 0.801506619929721, 0.780761747499825,
      0.759627103233651, 0.73820271041281, 0.716582893873948,
      0.694855827355588, 0.673103252282386, 0.651400340568036,
      0.629815676711063, 0.608411337149625, 0.587243047451526,
      0.566360400389608, 0.545807120255763, 0.525621360879725,
      0.50583602673333, 0.486479108216981, 0.467574023748514,
      0.449139962615167, 0.431192223719037, 0.41374254635902,
      0.396799430061806, 0.380368441215413, 0.36445250488482,
      0.349052180713631, 0.3341659222508, 0.319790319398603,
      0.305920323967661, 0.292549458556336, 0.279670009153575,
      0.267273202003828, 0.255349365376582, 0.243888076957123,
      0.232878297624443 };
   for (int i = 18 ; i <= ksz ; i++) {
      double xi_term = pow(1-pow(.95, i-1), 3);
      double target = (1 - xi_term * (1 - zeta_terms[i-1])) *
            pow(.99, i-1) * factor;
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   free(A);
   A = NULL;

   // Test special case N = 0.
   A = new_trunc_pol_A(0,0,0);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg == 1);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 1) target = .01;
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   clean_mem_prob();

}


void
test_error_new_trunc_pol_A
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *A;


   // Test error for 'new_trunc_pol_A()'.

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   A = new_trunc_pol_A(0,0,1);
   unredirect_stderr();
   reset_alloc();

   test_assert(A == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_B
(void)
{

   size_t ksz = 50;

   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *B = NULL;

   // Test i = 2, N = 1.
   B = new_trunc_pol_B(2,1);
   test_assert_critical(B != NULL);
   test_assert(B->monodeg == 2);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0.0;
      if (i == 2) target = 0.05*(1-0.05) * .99*.99;
      test_assert(fabs(B->coeff[i]-target) < 1e-9);
   }

   free(B);
   B = NULL;

   // Test i = 3, N = 1.
   B = new_trunc_pol_B(3,1);
   test_assert_critical(B != NULL);
   test_assert(B->monodeg == 3);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0.0;
      if (i == 3) target = 0.05*(1-0.05)*(1-0.05) * .99*.99*.99;
      test_assert(fabs(B->coeff[i]-target) < 1e-9);
   }

   free(B);
   B = NULL;

   // Test the special case N = 0.
   
   // Only the case i = 1 is non-zero.
   B = new_trunc_pol_B(2, 0);
   test_assert_critical(B != NULL);
   test_assert(B->monodeg == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(B->coeff[i] == 0);
   }

   free(B);
   B = NULL;

   B = new_trunc_pol_B(1, 0);
   test_assert_critical(B != NULL);
   test_assert(B->monodeg == 1);
   for (int i = 0 ; i <= 50 ; i++) {
      double target = 0;
      if (i == 1) target = 0.99;
      test_assert(fabs(B->coeff[i]-target) < 1e-9);
   }

   free(B);

   clean_mem_prob();

}


void
test_error_new_trunc_pol_B
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *B;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   B = new_trunc_pol_B(1,2);
   unredirect_stderr();
   reset_alloc();

   test_assert(B == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   redirect_stderr();
   // The error is that 'i = 0'.
   B = new_trunc_pol_B(0,2);
   unredirect_stderr();

   test_assert(B == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_trunc");

   clean_mem_prob();

}


void
test_new_trunc_pol_C
(void)
{

   // NOTE: The case m > N is not tested. Such polynomials are properly
   // defined, but they are not used in the present theory. They are
   // not forbidden, but also not used (and therefore not tested).
   
   size_t ksz = 50;

   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *C;
   
   C = new_trunc_pol_C(0, 1);
   test_assert_critical(C != NULL);
   test_assert(C->monodeg > ksz);
   for (int i = 0 ; i <= 16 ; i++) {
      double target = pow((1-.05) * .99, i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= ksz ; i++) {
      test_assert(C->coeff[i] == 0);
   }

   free(C);
   C = NULL;
   
   C = new_trunc_pol_C(1, 1);
   test_assert_critical(C != NULL);
   test_assert(C->monodeg > ksz);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = pow((1-.05) * .99, i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }

   free(C);
   C = NULL;
   
   C = new_trunc_pol_C(5, 10);
   test_assert_critical(C != NULL);
   test_assert(C->monodeg > ksz);
   for (int i = 0 ; i <= 16 ; i++) {
      double target = (1-pow(1-pow(.95,i),10)) * pow(.99,i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= ksz ; i++) {
      double target = (1-pow(1-pow(.95,i),5)) * pow(.99,i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }

   free(C);
   C = NULL;

   // Test the special case that 'N' is 0. Do not test
   // any other value than 'm' = 0 (see comment above).
   
   C = new_trunc_pol_C(0, 0);
   test_assert_critical(C != NULL);
   test_assert(C->monodeg == 0);
   test_assert(C->coeff[0] == 1.0);
   for (int i = 1 ; i <= ksz ; i++) {
      test_assert(C->coeff[i] == 0);
   }

   free(C);

   clean_mem_prob();

}


void
test_error_new_trunc_pol_C
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *C;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   C = new_trunc_pol_C(1, 3);
   unredirect_stderr();
   reset_alloc();

   test_assert(C == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_trunc_pol_update_add
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *w1 = new_zero_trunc_pol();
   trunc_pol_t *w2 = new_zero_trunc_pol();
   test_assert_critical(w1 != NULL);
   test_assert_critical(w2 != NULL);

   w1->monodeg = 1;
   w1->coeff[1] = 1;

   // Add null 'trunc_poly_t'.
   trunc_pol_update_add(w1, NULL);

   double array1[51] = {0,1};

   // Check that nothing has changed.
   test_assert(w1->monodeg == 1);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(w1->coeff[i] == array1[i]);
   }

   // Add another definition of null polynomial.
   trunc_pol_update_add(w1, w2);

   // Check that nothing has changed.
   test_assert(w1->monodeg == 1);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(w1->coeff[i] == array1[i]);
   }

   // Add polynomial of same degree.
   w2->monodeg = 1;
   w2->coeff[1] = 2;
   trunc_pol_update_add(w1, w2);

   double array2[51] = {0,3};
   test_assert(w1->monodeg == 1);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(w1->coeff[i] == array2[i]);
   }

   // Add non null polynomial of different degree.
   w2->monodeg = 2;
   w2->coeff[2] = 2;
   trunc_pol_update_add(w1, w2);

   // 'w1' is now {0,3,0,...} and 'w2' is {0,2,2,0,...}.
   double array3[51] = {0,5,2};
   test_assert(w1->monodeg == 51);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(w1->coeff[i] == array3[i]);
   }

   free(w1);
   free(w2);
   clean_mem_prob();

}



void
test_new_trunc_pol_D
(void)
{

   size_t ksz = 50;
   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *D;

   D = new_trunc_pol_D(1,0,1);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg > ksz);
   test_assert(D->coeff[0] == 0);
   for (int i = 1 ; i <= 16 ; i++) {
      double target = .01 * (1-.05/3) * pow(.99, i-1);
      test_assert(fabs(D->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= ksz ; i++) {
      test_assert(D->coeff[i] == 0);
   }
   
   free(D);
   D = NULL;

   D = new_trunc_pol_D(16,0,1);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg == 1);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 1) target = .01 * (1-.05/3);
      test_assert(fabs(D->coeff[i]-target) < 1e-9);
   }
   
   free(D);
   D = NULL;

   D = new_trunc_pol_D(13,18,20);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg > ksz);
   test_assert(D->coeff[0] == 0);
   // This is "p times omega sub m" computed with R.
   const double factor = 1.80897521494885666448e-30;
   for (int i = 1 ; i <= 4 ; i++) {
      double target = factor * .01 * pow(.99, i-1);
      // Use the relative error because the numbers are tiny.
      test_assert(fabs(D->coeff[i]/target - 1) < 1e-6);
   }
   for (int i = 5 ; i <= ksz ; i++) {
      test_assert(D->coeff[i] == 0);
   }
   
   free(D);
   D = NULL;

   // Test the case that 'm' is greater than 'N'.

   D = new_trunc_pol_D(1,10,5);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(D->coeff[i] == 0);
   }

   // Test the special cases that 'N' is 0.

   D = new_trunc_pol_D(1,0,0);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg > ksz);
   test_assert(D->coeff[0] == 0);
   for (int i = 1 ; i <= 16 ; i++) {
      double target = .01 * pow(.99, i-1);
      test_assert(fabs(D->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= ksz ; i++) {
      test_assert(D->coeff[i] == 0);
   }
   
   free(D);
   D = NULL;

   D = new_trunc_pol_D(16,0,0);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(D->coeff[i] == 0);
   }
   
   free(D);
   D = NULL;

   D = new_trunc_pol_D(1,1,0);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(D->coeff[i] == 0);
   }
   
   free(D);

   clean_mem_prob();

}


void
test_error_new_trunc_pol_D
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *D;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   D = new_trunc_pol_D(5,1,3);
   unredirect_stderr();
   reset_alloc();

   test_assert(D == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   redirect_stderr();
   // The error is that 'j' is greater than G-1
   D = new_trunc_pol_D(17,1,3);
   unredirect_stderr();
   reset_alloc();

   test_assert(D == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_tr");

   clean_mem_prob();

}


void
test_new_trunc_pol_E
(void)
{

   size_t ksz = 50;

   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *E;

   E = new_trunc_pol_E(16);
   test_assert_critical(E != NULL);

   // E polynomials with j = G-1 are monomials.
   test_assert(E->monodeg == 0);
   test_assert(E->coeff[0] == 1);
   for (int i = 1 ; i <= ksz ; i++) {
      test_assert(E->coeff[i] == 0);
   }

   free(E);
   E = NULL;

   E = new_trunc_pol_E(1);
   test_assert_critical(E != NULL);

   test_assert(E->monodeg > ksz);
   for (int i = 0 ; i <= 15 ; i++) {
      double target = pow(1-P,i);
      test_assert(fabs(E->coeff[i]-target) < 1e-9);
   }
   for (int i = 16 ; i <= ksz ; i++) {
      test_assert(E->coeff[i] == 0);
   }

   free(E);

   clean_mem_prob();

}


void
test_error_new_trunc_pol_E
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *E;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   E = new_trunc_pol_E(2);
   unredirect_stderr();
   reset_alloc();

   test_assert(E == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   redirect_stderr();
   // The error is that 'j' is greater than G-1
   E = new_trunc_pol_E(17);
   unredirect_stderr();

   test_assert(E == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_tr");

   clean_mem_prob();

}


void
test_new_null_matrix
(void)
{

   // Test a matrix of dimension 50.
   matrix_t *matrix = new_null_matrix(50);

   test_assert_critical(matrix != NULL);
   for (int i = 0 ; i < 50*50 ; i++) {
      test_assert(matrix->term[i] == NULL);
   }

   destroy_mat(matrix);

}


void
test_error_new_null_matrix
(void)
{

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   matrix_t *matrix = new_null_matrix(50);
   unredirect_stderr();
   reset_alloc();

   test_assert(matrix == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_n");

}


void
test_new_zero_matrix
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a matrix of dimension 10.
   matrix_t *matrix = new_zero_matrix(10);

   test_assert_critical(matrix != NULL);
   for (int i = 0 ; i < 10*10 ; i++) {
      trunc_pol_t *w = matrix->term[i];
      test_assert_critical(w != NULL);
      test_assert(w->monodeg == 0);
      for (int j = 0 ; j <= 10 ; j++) {
         test_assert(w->coeff[j] == 0);
      }
   }

   destroy_mat(matrix);

}


void
test_error_new_zero_matrix
(void)
{

   matrix_t *matrix;

   set_alloc_failure_countdown_to(0);
   redirect_stderr();
   matrix = new_zero_matrix(10);
   unredirect_stderr();
   reset_alloc();

   test_assert(matrix == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_n");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   matrix = new_zero_matrix(10);
   unredirect_stderr();
   reset_alloc();

   test_assert(matrix == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

}


void
test_new_matrix_M
(void)
{

   matrix_t *M;

   trunc_pol_t *A;
   trunc_pol_t *B;
   trunc_pol_t *C;
   trunc_pol_t *D;
   trunc_pol_t *E;

   size_t ksz = 50;

   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert_critical(success);

   // Test martrix M with 0 duplicate because it is a special case.
   M = new_matrix_M(0);
   test_assert_critical(M != NULL);

   const int dim0 = 18; // 17+0+1
   test_assert(M->dim == dim0);
   

   // -- First row -- //

   // First term (A polynomial).
   A = M->term[0];
   test_assert_critical(A != NULL);
   test_assert(A->monodeg == 1);
   test_assert(A->coeff[1] == .01);

   // Next 16 terms are B polynomials.
   for (int j = 1 ; j <= 16 ; j++) {
      B = M->term[j];
      test_assert_critical(B != NULL);
      if (j == 1) {
         test_assert(B->monodeg == 1);
         test_assert(B->coeff[1] == .99);
      }
      else {
         test_assert(iszero(B));
      }
   }
   
   // Final term of the row (C polynomial).
   C = M->term[17];
   test_assert_critical(C != NULL);
   test_assert(C->monodeg == 0);
   test_assert(C->coeff[0] == 1.0);
   for (int i = 1 ; i <= ksz ; i++) {
      test_assert(C->coeff[i] == 0);
   }

   // -- Next 'G-1' rows -- //
   
   for (int i = 1 ; i <= G-1 ; i++) {

      // First term (D polynomial).
      D = M->term[i*dim0];
      test_assert_critical(D != NULL);
      if (i == 1) {
         test_assert(D->monodeg > ksz);
         test_assert(D->coeff[0] == 0);
         for (int j = 1 ; j <= G-1 ; j++) {
            double target = pow(.99, j-1) * .01;
            test_assert(fabs(D->coeff[j]-target) < 1e-9);
         }
         for (int j = G ; j <= ksz ; j++) {
            test_assert(D->coeff[j] == 0);
         }
      }
      else {
         test_assert(iszero(D));
      }

      // Next 16 terms are null.
      for (int j = 1 ; j <= G-1 ; j++) {
         test_assert(M->term[(i+1)*dim0+j] == NULL);
      }

      // Final term of the row (E polynomial).
      E = M->term[i*dim0+17];
      test_assert_critical(E != NULL);
      if (i == 1) {
         test_assert(E->monodeg > ksz);
         for (int j = 0 ; j <= G-2 ; j++) {
            double target = pow(.99, j);
            test_assert(fabs(E->coeff[j]-target) < 1e-9);
         }
         for (int j = G-1 ; j <= ksz ; j++) {
            test_assert(E->coeff[j] == 0);
         }
      }
   }


   // -- Last row -- //
   
   for (int j = 0 ; j < dim0 ; j++) {
      test_assert(M->term[17*dim0+j] == NULL);
   }

   destroy_mat(M);
   M = NULL;

   // Test martrix M with one duplicate because the
   // polynomials are particularly simple in this case.
   M = new_matrix_M(1);
   test_assert_critical(M != NULL);

   const int dim1 = 19; // 17+1+1
   test_assert(M->dim == dim1);
   

   // -- First row -- //

   // First two terms (A polynomials).
   A = M->term[0];
   test_assert_critical(A != NULL);
   test_assert(A->monodeg > ksz);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = pow(.95 * .99, i-1) * .01 * (1-.05/3);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= ksz ; i++) {
      test_assert(A->coeff[i] == 0);
   }

   A = M->term[1];
   test_assert_critical(A != NULL);
   test_assert(A->monodeg > ksz);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= ksz ; i++) {
      double target = pow(.95 * .99, i-1) * .01 * .05/3;
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   // Next 16 terms are B polynomials.
   for (int j = 1 ; j <= 16 ; j++) {
      B = M->term[j+1];
      test_assert_critical(B != NULL);
      test_assert(B->monodeg == j);
      double target = (pow(.95,j-1)-pow(.95,j)) * pow(.99,j);
      test_assert(fabs(B->coeff[j]-target) < 1e-9);
   }
   
   // Final term of the row (C polynomial).
   C = M->term[18];
   test_assert_critical(C != NULL);
   test_assert(C->monodeg > ksz);
   for (int i = 0 ; i <= 16 ; i++) {
      double target = pow(.95 * .99, i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= ksz ; i++) {
      test_assert(C->coeff[i] == 0);
   }

   
   // -- Second row -- //

   // First two terms (A polynomials).
   A = M->term[1*dim1+0];
   test_assert_critical(A != NULL);
   test_assert(A->monodeg > ksz);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= ksz ; i++) {
      double target = pow(.95 * .99, i-1) * .01 * (1-.05/3);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   A = M->term[1*dim1+1];
   test_assert_critical(A != NULL);
   test_assert(A->monodeg > ksz);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= ksz ; i++) {
      double target = pow(.95 * .99, i-1) * .01 * .05/3;
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   // Next 16 terms are B polynomials.
   for (int j = 1 ; j <= 16 ; j++) {
      B = M->term[1*dim1+j+1];
      test_assert_critical(B != NULL);
      test_assert(B->monodeg == j);
      double target = (pow(.95,j-1)-pow(.95,j)) * pow(.99,j);
      test_assert(fabs(B->coeff[j]-target) < 1e-9);
   }
   
   // Final term of the row (C polynomial).
   C = M->term[1*dim1+18];
   test_assert_critical(C != NULL);
   test_assert(C->monodeg > ksz);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = pow(.95 * .99, i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }

   // -- Next 'G-1' rows -- //
   
   for (int i = 1 ; i <= G-1 ; i++) {
      // First two terms (D polynomials).
      D = M->term[(i+1)*dim1];
      test_assert_critical(D != NULL);
      if (i == 16)
         test_assert(D->monodeg == 1);
      else
         test_assert(D->monodeg > ksz);
      test_assert(D->coeff[0] == 0);
      for (int j = 1 ; j <= G-i ; j++) {
         double target = pow(.99, j-1) * .01 * (1-.05/3);
         test_assert(fabs(D->coeff[j]-target) < 1e-9);
      }
      for (int j = G-i+1 ; j <= ksz ; j++) {
         test_assert(D->coeff[j] == 0);
      }

      D = M->term[(i+1)*dim1+1];
      test_assert_critical(D != NULL);
      if (i == 16)
         test_assert(D->monodeg == 1);
      else
         test_assert(D->monodeg > ksz);
      test_assert(D->coeff[0] == 0);
      for (int j = 1 ; j <= G-i ; j++) {
         double target = pow(.99, j-1) * .01 * .05/3;
         test_assert(fabs(D->coeff[j]-target) < 1e-9);
      }
      for (int j = G-i+1 ; j <= ksz ; j++) {
         test_assert(D->coeff[j] == 0);
      }
      
      // Next 16 terms are null.
      for (int j = 1 ; j <= G-1 ; j++) {
         test_assert(M->term[(i+1)*dim1+j+1] == NULL);
      }

      // Final term of the row (E polynomial).
      E = M->term[(i+1)*dim1+18];
      test_assert_critical(E != NULL);
      if (i == 16)
         test_assert(E->monodeg == 0);
      else
         test_assert(E->monodeg > ksz);
      for (int j = 0 ; j <= G-i-1 ; j++) {
         double target = pow(.99, j);
         test_assert(fabs(E->coeff[j]-target) < 1e-9);
      }
      for (int j = G-i+1 ; j <= ksz ; j++) {
         test_assert(E->coeff[j] == 0);
      }
   }


   // -- Last row -- //
   
   for (int j = 0 ; j < dim1 ; j++) {
      test_assert(M->term[18*dim1+j] == NULL);
   }

   destroy_mat(M);
   clean_mem_prob();

}


void
test_error_new_matrix_M
(void)
{

   matrix_t *M;

   set_alloc_failure_countdown_to(0);
   redirect_stderr();
   M = new_matrix_M(1);
   unredirect_stderr();
   reset_alloc();

   test_assert(M == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_n");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   M = new_zero_matrix(1);
   unredirect_stderr();
   reset_alloc();

   test_assert(M == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

}


void
test_mem_seed_prob
(void)
{

   int success = set_params_mem_prob(17, 19, 0.01, 0.05);
   test_assert_critical(success);

   // The first terms can be computed directly.
   for (int i = 0 ; i < 17 ; i++) {
      test_assert(fabs(mem_seed_prob(0,i)-1) < 1e-9);
   }

   for (int i = 0 ; i < 17 ; i++) {
      test_assert(fabs(mem_seed_prob(1,i)-1) < 1e-9);
   }

   for (int i = 0 ; i < 17 ; i++) {
      test_assert(fabs(mem_seed_prob(2,i)-1) < 1e-9);
   }

   // N = 2, k = 17.
   double target_17 = 1-pow(.99,17);
   test_assert(fabs(mem_seed_prob(2,17)-target_17) < 1e-9);


   // N = 2, k = 18.
   double target_18;
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,2);
   test_assert(fabs(mem_seed_prob(2,18)-target_18) < 1e-9);

   // N = 3, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,3);
   // FIXME this test fails.
   test_assert(fabs(mem_seed_prob(3,18)-target_18) < 1e-9);

   // N = 4, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,4);
   // FIXME this test fails.
   test_assert(fabs(mem_seed_prob(4,18)-target_18) < 1e-9);

   // N = 20, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,20);
   // FIXME this test fails.
   test_assert(fabs(mem_seed_prob(20,18)-target_18) < 1e-9);


   // N = 1, k = 19.
   double target_19;
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * (1-pow(.95,18)*.05/3) - \
      2*.01*pow(.99,18) * (1-pow(.95,17)*.05/3) - \
      2*.01*.01*pow(.99,17) * (1-pow(.95,17)*.05/3) - \
      .01*.01*pow(.99,17) * (1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3);
   test_assert(fabs(mem_seed_prob(1,19)-target_19) < 1e-9);

   // N = 2, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,2) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,2) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,2) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,2);
   test_assert(fabs(mem_seed_prob(2,19)-target_19) < 1e-9);

   // N = 3, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,3) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,3) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,3) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,3);
   // FIXME this test fails.
   test_assert(fabs(mem_seed_prob(3,19)-target_19) < 1e-9);

   // N = 4, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,4) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,4) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,4) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,4);
   // FIXME this test fails.
   test_assert(fabs(mem_seed_prob(4,19)-target_19) < 1e-9);

   // N = 20, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,20) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,20) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,20) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,20);
   // FIXME this test fails.
   test_assert(fabs(mem_seed_prob(20,19)-target_19) < 1e-9);


   // Other test case with longer seeds.

   success = set_params_mem_prob(20, 50, 0.02, 0.05);
   test_assert_critical(success);

   for (int i = 0 ; i < 20 ; i++) {
      test_assert(fabs(mem_seed_prob(0,i)-1) < 1e-9);
      test_assert(fabs(mem_seed_prob(1,i)-1) < 1e-9);
      test_assert(fabs(mem_seed_prob(2,i)-1) < 1e-9);
   }

   const double target_20 = 1-pow(.98,20);
   test_assert(fabs(mem_seed_prob(0,20)-target_20) < 1e-9);
   test_assert(fabs(mem_seed_prob(1,20)-target_20) < 1e-9);
   test_assert(fabs(mem_seed_prob(2,20)-target_20) < 1e-9);
   test_assert(fabs(mem_seed_prob(3,20)-target_20) < 1e-9);

   double target_21;

   // Special case N = 0.
   target_21 = 1-pow(.98,21) - 2*.02*pow(.98,20);
   test_assert(fabs(mem_seed_prob(0,21)-target_21) < 1e-9);

   // Special case N = 1.
   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * (1-pow(.95,20)*.05/3);
   test_assert(fabs(mem_seed_prob(1,21)-target_21) < 1e-9);

   // Cases N > 1.
   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,2);
   test_assert(fabs(mem_seed_prob(2,21)-target_21) < 1e-9);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,3);
   // FIXME this test fails.
   test_assert(fabs(mem_seed_prob(3,21)-target_21) < 1e-9);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,4);
   // FIXME this test fails.
   test_assert(fabs(mem_seed_prob(4,21)-target_21) < 1e-9);

   clean_mem_prob();

}


void
test_error_mem_seed_prob
(void)
{

   double x;

   redirect_stderr();
   // The error is that the parameters are not initialized.
   x = mem_seed_prob(2,2);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `mem_");

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   redirect_stderr();
   // The error is that 'N' is greater than 'MAXN'.
   x = mem_seed_prob(1025,2);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `mem_");

   redirect_stderr();
   // The error is that 'k' is greater than specified value above.
   x = mem_seed_prob(2,51);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `mem_");

   set_alloc_failure_countdown_to(0);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   x = mem_seed_prob(2,10);
   unredirect_stderr();
   reset_alloc();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail (somewhere else).
   x = mem_seed_prob(2,10);
   unredirect_stderr();
   reset_alloc();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `new_n");

   clean_mem_prob();

}


void
test_misc_correctness
(void)
// The aim of this test is to compute the probability that a large
// read without MEM seed has exactly one mismatch. The principle is
// to count the paths that go through the states down/m exactly once
// and that do not go through the state down/0. Such paths can go
// through the states up/i before and / or after going through the
// state down/m, so the path consists of two to four segments.
//
// For a single error at position [j] less than or equal to 'G'
// from the left end of the read, the probability that the main
// thread is covered on the right is
//
//  I.    1 - (1 - (1-mu)^k-j * mu/3)^N = 1 - eta(k-j)^N.
//
// For a single error at position [j] more than 'G' from eiher
// end, the probability that the main thread is covered is 
//
//  II.   1 - ( eta(j-1)^N + eta(k-j)^N -
//                   (1-mu/3 + mu/3 xi(j-1) xi (k-j))^N ).
//
// The probability of the read is computed as the sum of those
// terms (the terms of the first kind are computed two times for
// symmetry). We also need to multiply by the probability that
// there is exactly one error at position [j].
//
// For term I. the state of the read before the error is irrelevant
// because there cannot be a seed. There is a seed if each duplicate
// has a mismatch between the error and the end of the read, which
// has probability (1 - (1-mu)^k-j * mu/3)^N = eta(k-j)^N. Term II.
// is more complex to derive. When the error is in the middle of the
// read, there can be a seed on the left or on the right of it. Thus
// there is a seed if all the duplicates leave the left side uncovered
// or all the duplicates leave the right side uncovered. The first
// event has probability (1 - (1-mu)^j-1 * mu/3)^N = eta(j-1)^N, the
// second has probability (1 - (1-mu)^k-j * mu/3)^N = eta(k-j)^N, but
// we have counted two times the event that all duplicates leave both
// ends uncovered. A duplicate covers neither sides if it matches the
// error -- probability 1-mu/3 -- or if it mismatches the error and
// has a mismatch on the left and on the right of it -- probability
// mu/3 (1 - (1-mu)^j-1) (1 - (1-mu)^k-j).
//
{
   
   size_t ksz = 35;

   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert_critical(success);

   matrix_t *M  = NULL;

   trunc_pol_t *tmp = tmp = new_zero_trunc_pol();
   test_assert_critical(tmp != NULL);

   for (int N = 1 ; N < 36 ; N++) {

      int dim = 17+N+1;

      // Transfer matrix.
      M = new_matrix_M(N);
      test_assert_critical(M != NULL);
      test_assert(M->dim == dim);

      // State vectors.
      trunc_pol_t **u = malloc(dim * sizeof(trunc_pol_t *));
      trunc_pol_t **v = malloc(dim * sizeof(trunc_pol_t *));
      test_assert_critical(u != NULL);
      test_assert_critical(v != NULL);
      for (int i = 0 ; i < dim ; i++) {
         u[i] = new_zero_trunc_pol();
         v[i] = new_zero_trunc_pol();
         test_assert_critical(u[i] != NULL);
         test_assert_critical(v[i] != NULL);
      }

      // One segment from head to all down/m states (except 'm' = 0).
      for (int j = 1 ; j <= N ; j++) {
         // Store the end state in 'S[j]'.
         trunc_pol_update_add(u[j], M->term[j]);
      }
      // One segment from head to all up/i and
      // one segment to all down/m states.
      for (int i = 1 ; i <= 16 ; i++) {
      for (int j = 1 ; j <= N ; j++) {
         // Store the end state in 'S[j]'.
         trunc_pol_update_add(u[j], trunc_pol_mult(tmp,
            M->term[N+i], M->term[(N+i)*dim+j]));
      }
      }
      // One segment from all down/m states to tail.
      for (int j = 1 ; j <= N ; j++) {
         trunc_pol_update_add(v[j], M->term[j*dim+(dim-1)]);
      }
      // One segment from all down/m to all up/i
      // and one segment to tail.
      for (int i = 1 ; i <= 16 ; i++) {
      for (int j = 1 ; j <= N ; j++) {
         trunc_pol_update_add(v[j], trunc_pol_mult(tmp,
            M->term[j*dim+(N+i)], M->term[(N+i)*dim+(dim-1)]));
      }
      }

      // Combine.
      trunc_pol_t *w = new_zero_trunc_pol();
      for (int j = 1 ; j <= N ; j++) {
         trunc_pol_update_add(w, trunc_pol_mult(tmp, u[j], v[j]));
      }

      double target = 0.0;
      // The position of the error is [j] (1-based).
      for (int j = 1 ; j <= 17 ; j++) {
         target += 2 * (1.0 - pow(eta(ksz-j),N));
      }
      for (int j = 18 ; j <= ksz-17 ; j++) {
         target +=  1 - ( pow(eta(j-1),N) + pow(eta(ksz-j),N) -
             pow(1-.05/3 + .05/3 * xi(j-1) * xi(ksz-j),N) );
      }
      // Multiply by the probability that there is one error.
      target *= .01 * pow(.99, ksz-1);
      test_assert(fabs(target - w->coeff[ksz]) < 1e-12);

      destroy_mat(M);
      for (int i = 0 ; i < dim ; i++) free(u[i]);
      for (int i = 0 ; i < dim ; i++) free(v[i]);
      free(u);
      free(v);
      free(w);

   }

   free(tmp);

   clean_mem_prob();

}


// Test cases for export.
const test_case_t test_cases_mem_seed_prob[] = {
   {"omega",                       test_omega},
   {"psi",                         test_psi},
   {"zeta",                        test_zeta},
   {"set_params_mem_prob",         test_set_params_mem_prob},
   {"error_set_params_mem_prob",   test_error_set_params_mem_prob},
   {"uninitialized_error",         test_uninitialized_error},
   {"new_zero_trunc_pol",          test_new_zero_trunc_pol},
   {"error_new_zero_trunc_pol",    test_error_new_zero_trunc_pol},
   {"trunc_pol_updated_add",       test_trunc_pol_update_add},
   {"trunc_pol_mult",              test_trunc_pol_mult},
   {"matrix_mult",                 test_matrix_mult},
   {"error_matrix_mult",           test_error_matrix_mult},
   {"new_trunc_pol_A",             test_new_trunc_pol_A},
   {"error_new_trunc_pol_A",       test_error_new_trunc_pol_A},
   {"new_trunc_pol_B",             test_new_trunc_pol_B},
   {"error_new_trunc_pol_B",       test_error_new_trunc_pol_B},
   {"new_trunc_pol_C",             test_new_trunc_pol_C},
   {"error_new_trunc_pol_C",       test_error_new_trunc_pol_C},
   {"new_trunc_pol_D",             test_new_trunc_pol_D},
   {"error_new_trunc_pol_D",       test_error_new_trunc_pol_D},
   {"new_trunc_pol_E",             test_new_trunc_pol_E},
   {"error_new_trunc_pol_E",       test_error_new_trunc_pol_E},
   {"new_null_matrix",             test_new_null_matrix},
   {"error_new_null_matrix",       test_error_new_null_matrix},
   {"new_zero_matrix",             test_new_zero_matrix},
   {"error_new_zero_matrix",       test_error_new_zero_matrix},
   {"new_matrix_M",                test_new_matrix_M},
   {"error_new_matrix_M",          test_error_new_matrix_M},
   {"mem_seed_prob",               test_mem_seed_prob},
   {"error_mem_seed_prob",         test_error_mem_seed_prob},
   {"misc_correctness",            test_misc_correctness},
   {NULL, NULL},
};
