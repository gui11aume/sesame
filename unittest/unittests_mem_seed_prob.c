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

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);


   // Test 'new_trunc_pol_A()'.

   // Test special case N = 0.
   
   clean_mem_prob();
   test_assert(0);

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
   // FIXME: give better paramters.
   A = new_trunc_pol_A(2,2,2);
   unredirect_stderr();
   reset_alloc();

   test_assert(A == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");


   // Test error for 'new_trunc_pol_A()'.

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // FIXME: give better parameters.
   A = new_trunc_pol_A(2,2,2);
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
   trunc_pol_t *B = NULL;

   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert_critical(success);

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

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a D polynomial of degree 10 with N = 2.
#if 0
   trunc_pol_t *C_ddown = new_trunc_pol_C_ddown(10, 2);
   test_assert_critical(C_ddown != NULL);
   test_assert(C_ddown->mono.deg == 0);
   test_assert(C_ddown->mono.coeff == 0);
   test_assert(C_ddown->coeff[0] == 0);
   const double omega = .01 * pow(1-.05/3,2);
   for (int i = 1 ; i <= 10 ; i++) {
      double target = omega * pow(.99, i-1);
      test_assert(fabs(C_ddown->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(C_ddown->coeff[i] == 0);
   }

   trunc_pol_t *C_down = new_trunc_pol_C_down(10, 2);
   test_assert_critical(C_down != NULL);
   test_assert(C_down->mono.deg == 0);
   test_assert(C_down->mono.coeff == 0);
   test_assert(C_down->coeff[0] == 0);
   const double _omega = .01 * (1 - pow(1-.05/3,2));
   for (int i = 1 ; i <= 10 ; i++) {
      double target = _omega * pow(.99, i-1);
      test_assert(fabs(C_down->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(C_down->coeff[i] == 0);
   }

   // Test the special cases N = 0.
   trunc_pol_t *C_ddown0 = new_trunc_pol_C_ddown(50, 0);
   trunc_pol_t *C_down0 = new_trunc_pol_C_down(50, 0);
   test_assert_critical(C_ddown0 != NULL);
   test_assert_critical(C_down0 != NULL);
   test_assert(C_ddown0->mono.deg == 0);
   test_assert(C_ddown0->mono.coeff == 0);
   test_assert(C_down0->mono.deg == 0);
   test_assert(C_down0->mono.coeff == 0);

   test_assert(C_ddown0->coeff[0] == 0);
   test_assert(C_down0->coeff[0] == 0);
   for (int i = 1 ; i <= 50 ; i++) {
      double target = pow(1-.01,i-1) * 0.01;
      test_assert(fabs(C_ddown0->coeff[i]-target) < 1e-9);
      test_assert(C_down0->coeff[i] == 0);
   }

   free(C_ddown);
   free(C_down);
   free(C_ddown0);
   free(C_down0);

#endif

   clean_mem_prob();
   test_assert(0);

}


void
test_error_new_trunc_pol_C
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *C;

   // Test errors in 'new_trunc_pol_C_ddown()'.
#if 0
   redirect_stderr();
   C = new_trunc_pol_C_ddown(0, 0);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_t");

   redirect_stderr();
   C = new_trunc_pol_C_ddown(0, 2);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_C_ddown(51, 2);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   C = new_trunc_pol_C_ddown(18, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(C == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");


   // Test errors in 'new_trunc_pol_C_down()'.

   redirect_stderr();
   C = new_trunc_pol_C_down(0, 0);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_t");

   redirect_stderr();
   C = new_trunc_pol_C_down(0, 2);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_C_down(51, 2);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   C = new_trunc_pol_C_down(18, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(C == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");
#endif

   clean_mem_prob();
   test_assert(0);

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
   double factor = 1.80897521494885666448e-30;
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

   trunc_pol_t *nspct;
   matrix_t *M;

   int success = set_params_mem_prob(17, 100, 0.01, 0.05);
   test_assert_critical(success);

#if 0
   const int dim = 17+2;

   // Test martrix M with one duplicate because the
   // polynomials are particularly simple in this case.
   M = new_matrix_M(1);
   test_assert_critical(M != NULL);
   test_assert(M->dim == 17+2);
   
   // Test first row.
   for (int j = 0 ; j < dim ; j++) {
      test_assert(M->term[j] == NULL);
   }

   // Second row, first term (T polynomial).
   test_assert_critical(M->term[dim] != NULL);
   nspct = M->term[dim];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   for (int i = 0 ; i <= 16 ; i++) {
      double target = pow(.99*.95,i);
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= 100 ; i++) {
      test_assert(nspct->coeff[i] == 0);
   }

   // Second row, second term.
   test_assert_critical(M->term[dim+1] != NULL);
   nspct = M->term[dim+1];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   test_assert(nspct->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = pow(.99*.95,i-1) * .01*(1-.05/3);
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 100 ; i++) {
      test_assert(nspct->coeff[i] == 0);
   }

   // Second row, third term.
   test_assert_critical(M->term[dim+2] != NULL);
   nspct = M->term[dim+2];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   test_assert(nspct->coeff[0] == 0);
   for (int i = 1 ; i <= 100 ; i++) {
      double target = pow(.99*.95,i-1) * .01*.05/3;
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }

   // Second row, u terms.
   for (int j = 1 ; j <= G-1 ; j++) {
      test_assert_critical(M->term[dim+2+j] != NULL);
      nspct = M->term[dim+2+j];
      double target = pow(.99*.95,j-1) * .99*.05;
      test_assert(nspct->mono.deg == j);
      test_assert(fabs(nspct->mono.coeff-target) < 1e-9);
      for (int i = 0 ; i <= 100 ; i++) {
         if (i == j)
            test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
         else
            test_assert(nspct->coeff[i] == 0);
      }
   }

   // Third row, first term (T polynomial).
   test_assert_critical(M->term[2*dim] != NULL);
   nspct = M->term[2*dim];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   for (int i = 0 ; i <= 100 ; i++) {
      double target = pow(.99*.95,i);
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }

   // Third row, second term.
   test_assert_critical(M->term[2*dim+1] != NULL);
   nspct = M->term[2*dim+1];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   test_assert(nspct->coeff[0] == 0);
   for (int i = 1 ; i <= 100 ; i++) {
      double target = pow(.99*.95,i-1) * .01*(1-.05/3);
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }

   // Third row, third term.
   test_assert_critical(M->term[2*dim+2] != NULL);
   nspct = M->term[2*dim+2];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   test_assert(nspct->coeff[0] == 0);
   for (int i = 1 ; i <= 100 ; i++) {
      double target = pow(.99*.95,i-1) * .01*.05/3;
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }

   // Third row, u terms.
   for (int j = 1 ; j <= G-1 ; j++) {
      test_assert_critical(M->term[2*dim+2+j] != NULL);
      nspct = M->term[2*dim+2+j];
      double target = pow(.99*.95,j-1) * .99*.05;
      test_assert(nspct->mono.deg == j);
      test_assert(fabs(nspct->mono.coeff-target) < 1e-9);
      for (int i = 0 ; i <= 100 ; i++) {
         if (i == j)
            test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
         else
            test_assert(nspct->coeff[i] == 0);
      }
   }

   // Middle series of rows.
   for (int j = 1 ; j <= G-1 ; j++) {
      // T polynomials.
      test_assert_critical(M->term[(j+2)*dim] != NULL);
      nspct = M->term[(j+2)*dim];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      for (int i = 0 ; i <= G-j-1 ; i++) {
         double target = pow(.99,i);
         test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
      }
      for (int i = G-j ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // C polynomials.
      test_assert_critical(M->term[(j+2)*dim+1] != NULL);
      nspct = M->term[(j+2)*dim+1];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      test_assert(nspct->coeff[0] == 0);
      for (int i = 1 ; i <= G-j ; i++) {
         double target = pow(.99,i-1) * .01*(1-.05/3);
         test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
      }
      for (int i = G-j+1 ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // Tilde C polynomials.
      test_assert_critical(M->term[(j+2)*dim+2] != NULL);
      nspct = M->term[(j+2)*dim+2];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      test_assert(nspct->coeff[0] == 0);
      for (int i = 1 ; i <= G-j ; i++) {
         double target = pow(.99,i-1) * .01*.05/3;
         test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
      }
      for (int i = G-j+1 ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // Rest of the rows.
      for (int i = 1 ; i <= G-1 ; i++) {
         test_assert(M->term[(j+2)*dim+2+i] == NULL);
      }

   }

   destroy_mat(M);



   // Test martrix M with two duplicates because the
   // polynomials are also simple in this case.
   M = new_matrix_M(2);
   test_assert_critical(M != NULL);
   
   // Test first row.
   for (int j = 0 ; j < dim ; j++) {
      test_assert(M->term[j] == NULL);
   }

   // Second row, first term (T polynomial).
   test_assert_critical(M->term[dim] != NULL);
   nspct = M->term[dim];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   for (int i = 0 ; i <= 16 ; i++) {
      double target = (1-pow(1-pow(.95,i),2))*pow(.99,i);
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= 100 ; i++) {
      test_assert(nspct->coeff[i] == 0);
   }

   // Second row, second term.
   test_assert_critical(M->term[dim+1] != NULL);
   nspct = M->term[dim+1];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   test_assert(nspct->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = (1-pow(1-pow(.95,i-1),2))*pow(.99,i-1) * \
         .01 * pow(1-.05/3,2);
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 100 ; i++) {
      test_assert(nspct->coeff[i] == 0);
   }

   // Second row, third term.
   test_assert_critical(M->term[dim+2] != NULL);
   nspct = M->term[dim+2];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   test_assert(nspct->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = (1-pow(1-pow(.95,i-1),2)) * pow(.99,i-1) * \
         0.01*(1-pow(1-.05/3,2));
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 100 ; i++) {
      double target = (1-pow(1-pow(.95,i-1)*.05/3,2)) * pow(.99,i-1)*.01;
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }

   // Second row, u terms.
   for (int j = 1 ; j <= G-1 ; j++) {
      test_assert_critical(M->term[dim+2+j] != NULL);
      nspct = M->term[dim+2+j];
      double target = (2*.05*pow(.95,j-1)*(1-pow(.95,j-1)) + \
         pow(.05*pow(.95,j-1),2)) * pow(.99,j);
      test_assert(nspct->mono.deg == j);
      test_assert(fabs(nspct->mono.coeff-target) < 1e-9);
      for (int i = 0 ; i <= 100 ; i++) {
         if (i == j)
            test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
         else
            test_assert(nspct->coeff[i] == 0);
      }
   }

   // Third row, first term (T polynomial).
   test_assert_critical(M->term[2*dim] != NULL);
   nspct = M->term[2*dim];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   test_assert(nspct->coeff[0] == 1.0);
   for (int i = 1 ; i <= 16 ; i++) {
      double xi_term = 2*pow(.95,i) - pow(.95,2*i);
      double target = pow(.99,i) * xi_term;
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= 100 ; i++) {
      // Declare constant in for loop to keep scope local.
      const double denom = 1-pow(1-.05/3,2);
      double alpha_i_sq = pow(1-pow(.95,i) * .05/3,2);
      double target = (1-alpha_i_sq) * pow(.99,i)/ denom;
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }

   // Third row, second term.
   test_assert_critical(M->term[2*dim+1] != NULL);
   nspct = M->term[2*dim+1];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   test_assert(nspct->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      // Declare constant in for loop to keep scope local.
      const double omega = .01 * pow(1-.05/3,2);
      double xi_term = 2*pow(.95,i-1) - pow(.95,2*i-2);
      double target = omega * xi_term * pow(.99,i-1);
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 100 ; i++) {
      // Declare constants in for loop to keep scope local.
      const double denom = 1-pow(1-.05/3,2);
      const double omega = .01 * pow(1-.05/3,2);
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double target = omega * (1-alpha_i_sq) * pow(.99,i-1) / denom;
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }

   // Third row, third term.
   test_assert_critical(M->term[2*dim+2] != NULL);
   nspct = M->term[2*dim+2];
   test_assert(nspct->mono.deg == 0);
   test_assert(nspct->mono.coeff == 0);
   test_assert(nspct->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      // Declare constant in for loop to keep scope local.
      const double _omega = .01 * (1-pow(1-.05/3,2));
      double xi_term = 2*pow(.95,i-1) - pow(.95,2*i-2);
      double target = _omega * xi_term * pow(.99,i-1);
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 100 ; i++) {
      // Declare constant in for loop to keep scope local.
      const double denom = 1-pow(1-.05/3,2);
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double alpha_0_sq = pow(1-.05/3,2);
      double zeta_ii_sq = pow(1-pow(.95,i-1) * .05/3 - 
         pow(.95,i-1) * .05/3 * (1-.05/3),2);
      double zeta_0i_sq = pow(1-.05/3 - 
         pow(.95,i-1) * .05/3 * (1-.05/3),2);
      double target = .01 * pow(.99, i-1) * (1-alpha_i_sq + \
         (alpha_i_sq - alpha_0_sq -zeta_ii_sq + zeta_0i_sq) / denom);
      test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
   }

   // Third row, u terms.
   for (int j = 1 ; j <= G-1 ; j++) {
      test_assert_critical(M->term[2*dim+2+j] != NULL);
      nspct = M->term[2*dim+2+j];
      double target = (2*.05*pow(.95,j-1)*(1-pow(.95,j-1)) + \
         pow(.05*pow(.95,j-1),2)) * pow(.99,j);
      test_assert(nspct->mono.deg == j);
      test_assert(fabs(nspct->mono.coeff-target) < 1e-9);
      for (int i = 0 ; i <= 100 ; i++) {
         if (i == j)
            test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
         else
            test_assert(nspct->coeff[i] == 0);
      }
   }

   // Middle series of rows.
   for (int j = 1 ; j <= G-1 ; j++) {
      // T polynomials.
      test_assert_critical(M->term[(j+2)*dim] != NULL);
      nspct = M->term[(j+2)*dim];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      for (int i = 0 ; i <= G-j-1 ; i++) {
         double target = pow(.99,i);
         test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
      }
      for (int i = G-j ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // C polynomials.
      test_assert_critical(M->term[(j+2)*dim+1] != NULL);
      nspct = M->term[(j+2)*dim+1];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      test_assert(nspct->coeff[0] == 0);
      for (int i = 1 ; i <= G-j ; i++) {
         double target = pow(.99,i-1) * .01 * pow(1-.05/3,2);
         test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
      }
      for (int i = G-j+1 ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // Tilde C polynomials.
      test_assert_critical(M->term[(j+2)*dim+2] != NULL);
      nspct = M->term[(j+2)*dim+2];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      test_assert(nspct->coeff[0] == 0);
      for (int i = 1 ; i <= G-j ; i++) {
         double target = pow(.99,i-1) * .01 * (1-pow(1-.05/3,2));
         test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
      }
      for (int i = G-j+1 ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // Rest of the rows.
      for (int i = 1 ; i <= G-1 ; i++) {
         test_assert(M->term[(j+2)*dim+2+i] == NULL);
      }

   }

   destroy_mat(M);
#endif

   clean_mem_prob();
   test_assert(0);

}


void
test_error_new_matrix_M
(void)
{

   matrix_t *M;

#if 0
   set_alloc_failure_countdown_to(0);
   redirect_stderr();
   M = new_matrix_M(10);
   unredirect_stderr();
   reset_alloc();

   test_assert(M == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_n");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   M = new_zero_matrix(10);
   unredirect_stderr();
   reset_alloc();

   test_assert(M == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");
#endif

}


void
test_mem_seed_prob
(void)
{

   int success = set_params_mem_prob(17, 19, 0.01, 0.05);
   test_assert_critical(success);

#if 0
   // Set full precision for these tests.
   set_mem_prob_max_precision_on();

   // The first terms can be computed directly.
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
   test_assert(fabs(mem_seed_prob(3,18)-target_18) < 1e-9);

   // N = 4, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,4);
   test_assert(fabs(mem_seed_prob(4,18)-target_18) < 1e-9);

   // N = 150, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,150);
   test_assert(fabs(mem_seed_prob(150,18)-target_18) < 1e-9);


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
   test_assert(fabs(mem_seed_prob(3,19)-target_19) < 1e-9);

   // N = 4, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,4) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,4) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,4) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,4);
   test_assert(fabs(mem_seed_prob(4,19)-target_19) < 1e-9);

   // N = 150, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,150) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,150) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,150) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,150);
   test_assert(fabs(mem_seed_prob(150,19)-target_19) < 1e-9);


   // Other test case with long seeds.

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
   test_assert(fabs(mem_seed_prob(3,21)-target_21) < 1e-9);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,4);
   test_assert(fabs(mem_seed_prob(4,21)-target_21) < 1e-9);
#endif

   clean_mem_prob();
   test_assert(0);

}


void
test_error_mem_seed_prob
(void)
{

   double x;

#if 0
   redirect_stderr();
   x = mem_seed_prob(2,2);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `mem_");

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   redirect_stderr();
   x = mem_seed_prob(1025,2);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `mem_");

   redirect_stderr();
   x = mem_seed_prob(2,51);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `mem_");

   set_alloc_failure_countdown_to(0);
   redirect_stderr();
   x = mem_seed_prob(2,10);
   unredirect_stderr();
   reset_alloc();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   x = mem_seed_prob(2,10);
   unredirect_stderr();
   reset_alloc();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `new_n");
#endif

   clean_mem_prob();
   test_assert(0);

}


void
test_misc_correctness
(void)
// The aim of this test is to compute the probability that a large
// read without MEM seed has exactly one mismatch. The principle is
// to count the paths that go through the state down exactly once
// and that do not go through the state ddown. Such paths can go
// through the states up before and / or after going through the
// state down, so the path consists of two to four segments.
//
// For a single error at position [j] less than or equal to [gamma],
// from the left end of the read, the probability that the main
// thread is covered is on the right is
//
//       1 - (1 - (1-mu)^k-j * mu/3)^N = 1 - alpha(k-j)^N.
//
// For a single error at position [j] more than [gamma] from eiher
// end, the probability that the main thread is covered is 
//
//     (1 - alpha(j-1)^N) (1 - alpha(k-j)^N / (1-(1-mu/3)^N). 
//
//     WARNING: The number above it only an approximation.
//
// The probability of the read is computed as the sum of those
// terms (the terms of the first kind are computed two times for
// symmetry). We also need to multiply by the probability that
// there is exactly one error at position [j].
{
   
   const int dim = 17+2; // Dimension of the matrix 'M' throughout.

   matrix_t *M  = NULL;
   trunc_pol_t *w1 = NULL;
   trunc_pol_t *w2 = NULL;
   trunc_pol_t *tmp = NULL;

   int success = set_params_mem_prob(17, 150, 0.01, 0.05);
   test_assert_critical(success);

   tmp = new_zero_trunc_pol();
   test_assert_critical(tmp != NULL);

#if 0
   for (int N = 0 ; N < 150 ; N++) {

      w1 = new_zero_trunc_pol();
      w2 = new_zero_trunc_pol();
      test_assert_critical(w1 != NULL);
      test_assert_critical(w2 != NULL);

      M = new_matrix_M(N);
      test_assert_critical(M != NULL);

      // One segment from head to down.
      trunc_pol_update_add(w1, M->term[1*dim+2]);
      for (int i = 1 ; i <= 16 ; i++) {
         // One segment from head to up(i) and one segment to down.
         trunc_pol_update_add(w1, trunc_pol_mult(tmp,
            M->term[1*dim+2+i], M->term[(2+i)*dim+2]));
      }
      // One segment from down to tail.
      trunc_pol_update_add(w2, M->term[2*dim+0]);
      for (int i = 1 ; i <= 16 ; i++) {
         // One segment from down to up(i) and one segment to tail.
         trunc_pol_update_add(w2, trunc_pol_mult(tmp,
            M->term[2*dim+2+i], M->term[(2+i)*dim+0]));
      }
      // Combine 'w1' and 'w2'.
      trunc_pol_mult(tmp, w1, w2);

      double target = 0.0;
      if (N > 0) {
         // The position of the error is [j] (1-based).
         for (int j = 1 ; j <= 17 ; j++) {
            target += 2 * (1.0 - aN(150-j));
         }
         for (int j = 18 ; j <= 133 ; j++) {
            target += (1-aN(j-1)) * (1-aN(150-j)) / (1-pow(1-.05/3,N));
         }
         target *= .01 * pow(.99,149);
      }

      test_assert(fabs(target - tmp->coeff[150]) < 1e-12);

      destroy_mat(M);
      free(w1);
      free(w2);

   }

   free(tmp);
#endif

   clean_mem_prob();
   test_assert(0);

}


// Test cases for export.
const test_case_t test_cases_mem_seed_prob[] = {
   {"set_params_mem_prob",         test_set_params_mem_prob},
   {"error_set_params_mem_prob",   test_error_set_params_mem_prob},
   {"uninitialized_error",         test_uninitialized_error},
   {"omega",                       test_omega},
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
