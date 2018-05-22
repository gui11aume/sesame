#include "unittest.h"
#include "compute_mem_prob.c"

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
   success = set_params_mem_prob(17, 50, 0.00, 0.05);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 2.
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.01, 0.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 3.
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 1.00, 0.05);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 4.
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.01, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 5.
   redirect_stderr();
   success = set_params_mem_prob(0, 50, 0.01, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 6.
   redirect_stderr();
   success = set_params_mem_prob(17, 0, 0.01, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 7.
   redirect_stderr();
   success = set_params_mem_prob(50, 20, 0.99, 0.95);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Test memory error.
   set_alloc_failure_rate_to(1.0);
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.01, 0.05);
   unredirect_stderr();
   reset_alloc();
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   return;

}


void
test_uninitialized_error
(void)
{

   // Do not call 'set_params_mem_prob()'.
   redirect_stderr();
   double x = compute_mem_prob(5, 20);
   unredirect_stderr();
   test_assert_stderr("[compute_mem_prob] error in function `comp");
   test_assert(x != x);

   return;

}


void
test_new_zero_trunc_pol
(void)
{

   size_t ksz = 50;
   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert(success);

   trunc_pol_t *a = new_zero_trunc_pol();
   test_assert_critical(a);

   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
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

   test_assert(trunc_pol_mult(a, NULL, NULL) == NULL);

   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   trunc_pol_t *b = new_zero_trunc_pol();
   test_assert_critical(b != NULL);

   test_assert(trunc_pol_mult(a, b, NULL) == NULL);

   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   test_assert(trunc_pol_mult(a, NULL, b) == NULL);

   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   trunc_pol_t *c = new_zero_trunc_pol();
   test_assert_critical(c != NULL);

   // Returns 'a' if arguments are not NULL.
   test_assert(trunc_pol_mult(a, b, c) == a);

   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }


   // Test multiplications between monomials (b = 5z and c = z^2).
   b->mono.deg = 1; b->mono.coeff = 5; b->coeff[1] = 5;
   c->mono.deg = 2; c->mono.coeff = 1; c->coeff[2] = 1;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->mono.deg == 3);
   test_assert(a->mono.coeff == 5);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = i == 3 ? 5 : 0;
      test_assert(a->coeff[i] == target);
   }


   // Test multiplications between a monomial and a
   // polynomial (b = 5z and c = z^2 + 2z^3).
   b->mono.deg = 1; b->mono.coeff = 5; b->coeff[1] = 5;
   c->mono.deg = 0; c->mono.coeff = 0; c->coeff[2] = 1; c->coeff[3] = 2;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 10;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 10;
      test_assert(a->coeff[i] == target);
   }

   // Test multiplications between two polynomials
   // (b = 5z + 3z^2 and c = z^2 + 2z^3).
   b->mono.deg = 0; b->mono.coeff = 0; b->coeff[1] = 5; b->coeff[2] = 3;
   c->mono.deg = 0; c->mono.coeff = 0; c->coeff[2] = 1; c->coeff[3] = 2;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 13;
      if (i == 5) target = 6;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 13;
      if (i == 5) target = 6;
      test_assert(a->coeff[i] == target);
   }

   free(a);
   free(b);
   free(c);
   clean_mem_prob();

   return;

}


void
test_special_matrix_mult
(void)
{


   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
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
      w1->mono.deg = i;
      w1->mono.coeff = 1;
      w1->coeff[w1->mono.deg] = w1->mono.coeff;
      w2->mono.deg = i;
      w2->mono.coeff = i;
      w2->coeff[w2->mono.deg] = w2->mono.coeff;
   }

   matrix_t *tmp1 = new_zero_matrix(2);
   test_assert_critical(tmp1 != NULL);

   special_matrix_mult(tmp1, mat1, mat2);

   const double tmp1_array1[51] = {0,0,0,2};
   const double tmp1_array2[51] = {0,1,0,0,3};
   const double tmp1_array3[51] = {0,0,0,0,0,2};
   const double tmp1_array4[51] = {0,0,0,1,0,0,3};

   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tmp1->term[i]->mono.deg == 0);
      test_assert(tmp1->term[i]->mono.coeff == 0);
   }
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(tmp1->term[0]->coeff[i] == tmp1_array1[i]);
      test_assert(tmp1->term[1]->coeff[i] == tmp1_array2[i]);
      test_assert(tmp1->term[2]->coeff[i] == tmp1_array3[i]);
      test_assert(tmp1->term[3]->coeff[i] == tmp1_array4[i]);
   }

   matrix_t *tmp2 = new_zero_matrix(2);
   test_assert_critical(tmp2 != NULL);

   special_matrix_mult(tmp2, mat2, mat1);

   const double tmp2_array1[51] = {0,0,0,1};
   const double tmp2_array2[51] = {0,0,0,0,1};
   const double tmp2_array3[51] = {0,0,2,0,0,3};
   const double tmp2_array4[51] = {0,0,0,2,0,0,3};

   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tmp2->term[i]->mono.deg == 0);
      test_assert(tmp2->term[i]->mono.coeff == 0);
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
test_error_special_matrix_mult
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
   test_assert(special_matrix_mult(tmp, mat1, mat2) == NULL);
   unredirect_stderr();

   test_assert_stderr("[compute_mem_prob] error in function `spe");

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
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   w = new_zero_trunc_pol();
   unredirect_stderr();
   reset_alloc();

   test_assert(w == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

}


void
test_new_trunc_pol_A
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert(success);


   // Test 'new_trunc_pol_A_ddown()'.

   trunc_pol_t *A_ddown = new_trunc_pol_A_ddown(2);
   test_assert_critical(A_ddown != NULL);
   test_assert(A_ddown->mono.deg == 0);
   test_assert(A_ddown->mono.coeff == 0);
   test_assert(A_ddown->coeff[0] == 0);
   const double omega = .01 * pow(1-.05/3,2);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = omega * pow(.99, i-1) * (1-pow(1-pow(.95,i-1),2)); 
      test_assert(fabs(A_ddown->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 50 ; i++) {
      test_assert(A_ddown->coeff[i] == 0);
   }

   trunc_pol_t *A_down = new_trunc_pol_A_down(2);
   test_assert_critical(A_down != NULL);
   test_assert(A_down->mono.deg == 0);
   test_assert(A_down->mono.coeff == 0);
   test_assert(A_down->coeff[0] == 0);
   const double _omega = .01 * (1-pow(1-.05/3,2));
   for (int i = 1 ; i <= 17 ; i++) {
      double target = _omega * pow(.99, i-1) * (1-pow(1-pow(.95,i-1),2)); 
      test_assert(fabs(A_down->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double target = 0.01 * pow(.99, i-1) * (1-alpha_i_sq);
      test_assert(fabs(A_down->coeff[i]-target) < 1e-9);
   }

   // Test special case N = 0.
   trunc_pol_t *A_ddown0 = new_trunc_pol_A_ddown(0);
   test_assert_critical(A_ddown0 != NULL);
   test_assert(A_ddown0->mono.deg == 1);
   test_assert(A_ddown0->mono.coeff == .01);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(A_ddown0->coeff[i] == (i == 1 ? 0.01 : 0));
   }

   trunc_pol_t *A_down0 = new_trunc_pol_A_down(0);
   test_assert_critical(A_down0 != NULL);
   test_assert(A_down0->mono.deg == 0);
   test_assert(A_down0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(A_down0->coeff[i] == 0);
   }

   free(A_ddown);
   free(A_down);
   free(A_ddown0);
   free(A_down0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_A
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *A;


   // Test error for 'new_trunc_pol_A_ddown()'.

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   A = new_trunc_pol_A_ddown(2);
   unredirect_stderr();
   reset_alloc();

   test_assert(A == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");


   // Test error for 'new_trunc_pol_A_ddown()'.

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   A = new_trunc_pol_A_down(2);
   unredirect_stderr();
   reset_alloc();

   test_assert(A == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_B
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert(success);

   // Test 'new_trunc_pol_B_ddown()'.
   trunc_pol_t *B_ddown = new_trunc_pol_B_ddown(2);
   test_assert_critical(B_ddown != NULL);
   test_assert(B_ddown->mono.deg == 0);
   test_assert(B_ddown->mono.coeff == 0);
   test_assert(B_ddown->coeff[0] == 0);
   const double omega = .01 * pow(1-.05/3,2);
   const double denom = 1-pow(1-.05/3,2);
   for (int i = 1 ; i <= 17 ; i++) {
      double xi_term = 2*pow(.95,i-1) - pow(.95,2*i-2);
      double target = omega * pow(.99, i-1) * xi_term;
      test_assert(fabs(B_ddown->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double target = omega * pow(.99, i-1) * (1-alpha_i_sq) / denom;
      test_assert(fabs(B_ddown->coeff[i]-target) < 1e-9);
   }

   // Test 'new_trunc_pol_B_down()'.
   trunc_pol_t *B_down = new_trunc_pol_B_down(2);
   test_assert_critical(B_down != NULL);
   test_assert(B_down->mono.deg == 0);
   test_assert(B_down->mono.coeff == 0);
   test_assert(B_down->coeff[0] == 0);
   const double _omega = .01 * (1-pow(1-.05/3,2));
   for (int i = 1 ; i <= 17 ; i++) {
      double xi_term = 2*pow(.95,i-1) - pow(.95,2*i-2);
      double target = _omega * pow(.99, i-1) * xi_term;
      test_assert(fabs(B_down->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double alpha_0_sq = pow(1-.05/3,2);
      double zeta_ii_sq = pow(1-pow(.95,i-1) * .05/3 - 
         pow(.95,i-1) * .05/3 * (1-.05/3),2);
      double zeta_0i_sq = pow(1-.05/3 - 
         pow(.95,i-1) * .05/3 * (1-.05/3),2);
      double target = .01 * pow(.99, i-1) * (1-alpha_i_sq + \
         (alpha_i_sq - alpha_0_sq -zeta_ii_sq + zeta_0i_sq) / denom);
      test_assert(fabs(B_down->coeff[i]-target) < 1e-9);
   }

   // Test the special case N = 0.
   trunc_pol_t *B_ddown0 = new_trunc_pol_B_ddown(0);
   test_assert_critical(B_ddown0 != NULL);
   test_assert(B_ddown0->mono.deg == 0);
   test_assert(B_ddown0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(B_ddown0->coeff[i] == 0);
   }

   trunc_pol_t *B_down0 = new_trunc_pol_B_down(0);
   test_assert_critical(B_down0 != NULL);
   test_assert(B_down0->mono.deg == 0);
   test_assert(B_down0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(B_down0->coeff[i] == 0);
   }

   free(B_ddown);
   free(B_down);
   free(B_ddown0);
   free(B_down0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_B
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *B;


   // Test error in 'new_trunc_pol_B_ddown()'.
   
   set_alloc_failure_rate_to(1);
   redirect_stderr();
   B = new_trunc_pol_B_ddown(2);
   unredirect_stderr();
   reset_alloc();

   test_assert(B == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");


   // Test error in 'new_trunc_pol_B_down()'.
   
   set_alloc_failure_rate_to(1);
   redirect_stderr();
   B = new_trunc_pol_B_down(2);
   unredirect_stderr();
   reset_alloc();

   test_assert(B == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_C
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a D polynomial of degree 10 with N = 2.
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
   clean_mem_prob();

}


void
test_error_new_trunc_pol_C
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *C;


   // Test errors in 'new_trunc_pol_C_ddown()'.

   redirect_stderr();
   C = new_trunc_pol_C_ddown(0, 0);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   C = new_trunc_pol_C_ddown(0, 2);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_C_ddown(51, 2);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   C = new_trunc_pol_C_ddown(18, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");


   // Test errors in 'new_trunc_pol_C_down()'.

   redirect_stderr();
   C = new_trunc_pol_C_down(0, 0);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   C = new_trunc_pol_C_down(0, 2);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_C_down(51, 2);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   C = new_trunc_pol_C_down(18, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

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

   w1->mono.deg = 1;
   w1->mono.coeff = 1;
   w1->coeff[1] = 1;

   w2->mono.deg = 2;
   w2->mono.coeff = 2;
   w2->coeff[2] = 2;

   trunc_pol_update_add(w1, w2);

   double array[51] = {0,1,2};
   test_assert(w1->mono.deg == 0);
   test_assert(w1->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(w1->coeff[i] == array[i]);
   }

   trunc_pol_update_add(w1, NULL);

   // Check that nothing has changed.
   test_assert(w1->mono.deg == 0);
   test_assert(w1->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(w1->coeff[i] == array[i]);
   }

   free(w1);
   free(w2);
   clean_mem_prob();

}


void
test_new_trunc_pol_u
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a u polynomial of degree 10 with N = 2.
   trunc_pol_t *u = new_trunc_pol_u(10, 2);
   test_assert_critical(u != NULL);

   const double target = pow(1-.01,10) * \
      (pow(1-pow(1-.05,10),2) - pow(1-pow(1-.05,9),2));
   test_assert(u->mono.deg == 10);
   test_assert(fabs(u->mono.coeff-target) < 1e-9);
   for (int i = 0 ; i <= 50 ; i++) {
      if (i == 10)
         test_assert(fabs(u->coeff[i]-target) < 1e-9);
      else
         test_assert(u->coeff[i] == 0);
   }

   // Test the special case N = 0.
   trunc_pol_t *u0 = new_trunc_pol_u(1, 0);
   test_assert_critical(u0 != NULL);
   test_assert(u0->mono.deg == 1);
   test_assert(fabs(u0->mono.coeff-.99) < 1e-9);
   for (int i = 0 ; i <= 50 ; i++) {
      if (i == 1)
         test_assert(fabs(u0->coeff[i]-.99) < 1e-9);
      else
         test_assert(u0->coeff[i] == 0);
   }
   
   free(u);
   free(u0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_u
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *u;

   redirect_stderr();
   u = new_trunc_pol_u(0, 0);
   unredirect_stderr();

   test_assert(u == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   u = new_trunc_pol_u(0, 2);
   unredirect_stderr();

   test_assert(u == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_u(18, 2);
   unredirect_stderr();

   test_assert(u == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   u = new_trunc_pol_u(10, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(u == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_T_down
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a T_down polynomial with N = 2.
   trunc_pol_t *T_down = new_trunc_pol_T_down(2);
   test_assert_critical(T_down != NULL);

   const double denom = 1-pow(1-.05/3,2);
   test_assert(T_down->mono.deg == 0);
   test_assert(T_down->mono.coeff == 0);
   test_assert(T_down->coeff[0] == 1);
   for (int i = 1 ; i <= 16 ; i++) {
      double xi_term = 2*pow(.95,i) - pow(.95,2*i);
      double target = pow(.99,i) * xi_term;
      test_assert(fabs(T_down->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i) * .05/3,2);
      double target = pow(.99,i) * (1-alpha_i_sq) / denom;
      test_assert(fabs(T_down->coeff[i]-target) < 1e-9);
   }

   // Test the special cases N = 0.
   trunc_pol_t *T_down0 = new_trunc_pol_T_down(0);
   test_assert_critical(T_down0 != NULL);
   test_assert(T_down0->mono.deg == 0);
   test_assert(T_down0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(T_down0->coeff[i] == 0);
   }
   
   free(T_down);
   free(T_down0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_T_down
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *T_down;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   T_down = new_trunc_pol_T_down(2);
   unredirect_stderr();
   reset_alloc();

   test_assert(T_down == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_T_ddown
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a T_ddown polynomial with N = 2.
   trunc_pol_t *T_ddown = new_trunc_pol_T_ddown(2);
   test_assert_critical(T_ddown != NULL);

   test_assert(T_ddown->mono.deg == 0);
   test_assert(T_ddown->mono.coeff == 0);
   test_assert(T_ddown->coeff[0] == 1);
   for (int i = 1 ; i <= 16 ; i++) {
      double target = (1-pow(1-pow(1-.05,i),2)) * pow(1-.01,i);
      test_assert(fabs(T_ddown->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= 50 ; i++) {
      test_assert(T_ddown->coeff[i] == 0);
   }

   // Test the special cases N = 0.
   trunc_pol_t *T_ddown0 = new_trunc_pol_T_ddown(0);
   test_assert_critical(T_ddown0 != NULL);
   test_assert(T_ddown0->mono.deg == 0);
   test_assert(T_ddown0->mono.coeff == 1);
   test_assert(T_ddown0->coeff[0] == 1);
   for (int i = 1 ; i <= 50 ; i++) {
      test_assert(T_ddown0->coeff[i] == 0);
   }
   
   free(T_ddown);
   free(T_ddown0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_T_ddown
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *T_ddown;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   T_ddown = new_trunc_pol_T_ddown(2);
   unredirect_stderr();
   reset_alloc();

   test_assert(T_ddown == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_T_up
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a T_up polynomial of degree 10 with N = 2.
   trunc_pol_t *T_up = new_trunc_pol_T_up(10, 2);
   test_assert_critical(T_up != NULL);

   test_assert(T_up->mono.deg == 0);
   test_assert(T_up->mono.coeff == 0);
   test_assert(T_up->coeff[0] == 1);
   for (int i = 1 ; i <= 10 ; i++) {
      double target = pow(1-.01,i);
      test_assert(fabs(T_up->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(T_up->coeff[i] == 0);
   }

   // Test the special cases N = 0.
   trunc_pol_t *T_up0 = new_trunc_pol_T_up(10, 0);
   test_assert_critical(T_up0 != NULL);
   test_assert(T_up0->mono.deg == 0);
   test_assert(T_up0->mono.coeff == 0);
   for (int i = 0 ; i <= 10 ; i++) {
      double target = pow(1-.01,i);
      test_assert(fabs(T_up0->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(T_up0->coeff[i] == 0);
   }
   
   free(T_up);
   free(T_up0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_T_up
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *T_up;

   redirect_stderr();
   T_up = new_trunc_pol_T_up(17, 2);
   unredirect_stderr();

   test_assert(T_up == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   T_up = new_trunc_pol_T_up(10, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(T_up == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

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
   test_assert_stderr("[compute_mem_prob] error in function `new_n");

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
      test_assert(w->mono.deg == 0);
      test_assert(w->mono.coeff == 0);
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
   test_assert_stderr("[compute_mem_prob] error in function `new_n");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   matrix = new_zero_matrix(10);
   unredirect_stderr();
   reset_alloc();

   test_assert(matrix == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

}


void
test_new_matrix_M
(void)
{

   trunc_pol_t *nspct;
   matrix_t *M;

   int success = set_params_mem_prob(17, 100, 0.01, 0.05);
   test_assert_critical(success);

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

   // Write the other tests.
   int test_case_incomplete = 0;
   test_assert(test_case_incomplete);

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
   M = new_matrix_M(10);
   unredirect_stderr();
   reset_alloc();

   test_assert(M == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_n");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   M = new_zero_matrix(10);
   unredirect_stderr();
   reset_alloc();

   test_assert(M == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

}


void
test_compute_mem_prob
(void)
{

   int success = set_params_mem_prob(17, 19, 0.01, 0.05);
   test_assert_critical(success);

   // Set full precision for these tests.
   set_mem_prob_max_precision_on();

   // The first terms can be computed directly.
   for (int i = 0 ; i < 17 ; i++) {
      test_assert(fabs(compute_mem_prob(2,i)-1) < 1e-9);
   }

   // N = 2, k = 17.
   double target_17 = 1-pow(.99,17);
   test_assert(fabs(compute_mem_prob(2,17)-target_17) < 1e-9);


   // N = 2, k = 18.
   double target_18;
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,2);
   test_assert(fabs(compute_mem_prob(2,18)-target_18) < 1e-9);

   // N = 3, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,3);
   test_assert(fabs(compute_mem_prob(3,18)-target_18) < 1e-9);

   // N = 4, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,4);
   test_assert(fabs(compute_mem_prob(4,18)-target_18) < 1e-9);

   // N = 150, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,150);
   test_assert(fabs(compute_mem_prob(150,18)-target_18) < 1e-9);


   // N = 1, k = 19.
   double target_19;
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * (1-pow(.95,18)*.05/3) - \
      2*.01*pow(.99,18) * (1-pow(.95,17)*.05/3) - \
      2*.01*.01*pow(.99,17) * (1-pow(.95,17)*.05/3) - \
      .01*.01*pow(.99,17) * (1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3);
   test_assert(fabs(compute_mem_prob(1,19)-target_19) < 1e-9);

   // N = 2, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,2) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,2) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,2) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,2);
   test_assert(fabs(compute_mem_prob(2,19)-target_19) < 1e-9);

   // N = 3, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,3) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,3) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,3) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,3);
   test_assert(fabs(compute_mem_prob(3,19)-target_19) < 1e-9);

   // N = 4, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,4) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,4) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,4) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,4);
   test_assert(fabs(compute_mem_prob(4,19)-target_19) < 1e-9);

   // N = 150, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,150) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,150) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,150) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,150);
   test_assert(fabs(compute_mem_prob(150,19)-target_19) < 1e-9);


   // Other test case with long seeds.

   success = set_params_mem_prob(20, 50, 0.02, 0.05);
   test_assert_critical(success);

   for (int i = 0 ; i < 20 ; i++) {
      test_assert(fabs(compute_mem_prob(0,i)-1) < 1e-9);
      test_assert(fabs(compute_mem_prob(1,i)-1) < 1e-9);
      test_assert(fabs(compute_mem_prob(2,i)-1) < 1e-9);
   }

   const double target_20 = 1-pow(.98,20);
   test_assert(fabs(compute_mem_prob(0,20)-target_20) < 1e-9);
   test_assert(fabs(compute_mem_prob(1,20)-target_20) < 1e-9);
   test_assert(fabs(compute_mem_prob(2,20)-target_20) < 1e-9);
   test_assert(fabs(compute_mem_prob(3,20)-target_20) < 1e-9);

   double target_21;

   // Special case N = 0.
   target_21 = 1-pow(.98,21) - 2*.02*pow(.98,20);
   test_assert(fabs(compute_mem_prob(0,21)-target_21) < 1e-9);

   // Special case N = 1.
   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * (1-pow(.95,20)*.05/3);
   test_assert(fabs(compute_mem_prob(1,21)-target_21) < 1e-9);

   // Cases N > 1.
   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,2);
   test_assert(fabs(compute_mem_prob(2,21)-target_21) < 1e-9);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,3);
   test_assert(fabs(compute_mem_prob(3,21)-target_21) < 1e-9);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,4);
   test_assert(fabs(compute_mem_prob(4,21)-target_21) < 1e-9);

   clean_mem_prob();

}


void
test_error_compute_mem_prob
(void)
{

   double x;

   redirect_stderr();
   x = compute_mem_prob(2,2);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[compute_mem_prob] error in function `comp");

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   redirect_stderr();
   x = compute_mem_prob(1025,2);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[compute_mem_prob] error in function `comp");

   redirect_stderr();
   x = compute_mem_prob(2,51);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[compute_mem_prob] error in function `comp");

   set_alloc_failure_countdown_to(0);
   redirect_stderr();
   x = compute_mem_prob(2,10);
   unredirect_stderr();
   reset_alloc();

   test_assert(x != x);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   x = compute_mem_prob(2,10);
   unredirect_stderr();
   reset_alloc();

   test_assert(x != x);
   test_assert_stderr("[compute_mem_prob] error in function `new_n");

   clean_mem_prob();

}


void
test_misc_exactness
(void)
{

   matrix_t *M1  = NULL;
   matrix_t *M2 = NULL;

   int success = set_params_mem_prob(17, 150, 0.01, 0.05);
   test_assert_critical(success);

   // Test with N = 0.
   M1 = new_matrix_M(0);
   M2 = new_zero_matrix(17+2);
   test_assert_critical(M1 != NULL);
   test_assert_critical(M2 != NULL);

   special_matrix_mult(M2, M1, M1);

   test_assert(M1->term[17+2]->coeff[150] == 0);
   test_assert(M2->term[17+2]->coeff[150] == 0);

   int test_case_incomplete = 0;
   test_assert(test_case_incomplete);

   destroy_mat(M1);
   destroy_mat(M2);
   clean_mem_prob();

}


// Test cases for export.
const test_case_t test_cases_compute_mem_prob[] = {
   {"set_params_mem_prob",         test_set_params_mem_prob},
   {"error_set_params_mem_prob",   test_error_set_params_mem_prob},
   {"uninitialized_error",         test_uninitialized_error},
   {"new_zero_trunc_pol",          test_new_zero_trunc_pol},
   {"error_new_zero_trunc_pol",    test_error_new_zero_trunc_pol},
   {"trunc_pol_updated_add",       test_trunc_pol_update_add},
   {"trunc_pol_mult",              test_trunc_pol_mult},
   {"special_matrix_mult",         test_special_matrix_mult},
   {"error_special_matrix_mult",   test_error_special_matrix_mult},
   {"new_trunc_pol_A",             test_new_trunc_pol_A},
   {"error_new_trunc_pol_A",       test_error_new_trunc_pol_A},
   {"new_trunc_pol_B",             test_new_trunc_pol_B},
   {"error_new_trunc_pol_B",       test_error_new_trunc_pol_B},
   {"new_trunc_pol_C",             test_new_trunc_pol_C},
   {"error_new_trunc_pol_C",       test_error_new_trunc_pol_C},
   {"new_trunc_pol_u",             test_new_trunc_pol_u},
   {"error_new_trunc_pol_u",       test_error_new_trunc_pol_u},
   {"new_trunc_pol_T_down",        test_new_trunc_pol_T_down},
   {"error_new_trunc_pol_T_down",  test_error_new_trunc_pol_T_down},
   {"new_trunc_pol_T_ddown",       test_new_trunc_pol_T_ddown},
   {"error_new_trunc_pol_T_ddown", test_error_new_trunc_pol_T_ddown},
   {"new_trunc_pol_T_up",          test_new_trunc_pol_T_up},
   {"error_new_trunc_pol_T_up",    test_error_new_trunc_pol_T_up},
   {"new_null_matrix",             test_new_null_matrix},
   {"error_new_null_matrix",       test_error_new_null_matrix},
   {"new_zero_matrix",             test_new_zero_matrix},
   {"error_new_zero_matrix",       test_error_new_zero_matrix},
   {"new_matrix_M",                test_new_matrix_M},
   {"error_new_matrix_M",          test_error_new_matrix_M},
   {"compute_mem_prob",            test_compute_mem_prob},
   {"error_compute_mem_prob",      test_error_compute_mem_prob},
   {"misc_exactness",              test_misc_exactness},
   {NULL, NULL},
};
