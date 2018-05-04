#include "unittest.h"
#include "compute_mem_prob.c"

void
test_set_params_mem_prob
(void)
{

   int success;

   // Case 1.
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

   // Case 2.
   set_params_mem_prob(50, 20, 0.99, 0.95);
   test_assert_critical(success);

   test_assert(G == 50);
   test_assert(K == 20);
   test_assert(P == 0.99);
   test_assert(U == 0.95);

   test_assert(KSZ == sizeof(trunc_pol_t) + 21*sizeof(double));

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
test_new_trunc_pol_A
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert(success);

   trunc_pol_t *A = new_trunc_pol_A(17, 2, NO);
   test_assert_critical(A != NULL);
   test_assert(A->mono.deg == 0);
   test_assert(A->mono.coeff == 0);
   test_assert(A->coeff[0] == 0);
   const double omega = .01 * pow(1-.05/3,2);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = omega * pow(.99, i-1) * (1-pow(1-pow(.95,i-1),2)); 
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 50 ; i++) {
      test_assert(A->coeff[i] == 0);
   }

   trunc_pol_t *_A = new_trunc_pol_A(50, 2, YES);
   test_assert_critical(_A != NULL);
   test_assert(_A->mono.deg == 0);
   test_assert(_A->mono.coeff == 0);
   test_assert(_A->coeff[0] == 0);
   const double _omega = .01 * (1-pow(1-.05/3,2));
   for (int i = 1 ; i <= 17 ; i++) {
      double target = _omega * pow(.99, i-1) * (1-pow(1-pow(.95,i-1),2)); 
      test_assert(fabs(_A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double target = 0.01 * pow(.99, i-1) * (1-alpha_i_sq);
      test_assert(fabs(_A->coeff[i]-target) < 1e-9);
   }

   // Test special case N = 0.
   trunc_pol_t *A0 = new_trunc_pol_A(50, 0, NO);
   test_assert_critical(A0 != NULL);
   test_assert(A0->mono.deg == 1);
   test_assert(A0->mono.coeff == .01);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(A0->coeff[i] == (i == 1 ? 0.01 : 0));
   }

   trunc_pol_t *_A0 = new_trunc_pol_A(50, 0, YES);
   test_assert_critical(_A0 != NULL);
   test_assert(_A0->mono.deg == 0);
   test_assert(_A0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(_A0->coeff[i] == 0);
   }

   free(A);
   free(_A);
   free(A0);
   free(_A0);
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
test_error_new_trunc_pol_A
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *A;

   redirect_stderr();
   A = new_trunc_pol_A(0, 0, NO);
   unredirect_stderr();

   test_assert(A == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   A = new_trunc_pol_A(0, 2, NO);
   unredirect_stderr();

   test_assert(A == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_A(51, 2, NO);
   unredirect_stderr();

   test_assert(A == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   A = new_trunc_pol_A(18, 2, NO);
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

   trunc_pol_t *B = new_trunc_pol_B(50, 2, NO);
   test_assert_critical(B != NULL);
   test_assert(B->mono.deg == 0);
   test_assert(B->mono.coeff == 0);
   test_assert(B->coeff[0] == 0);
   const double omega = .01 * pow(1-.05/3,2);
   const double denom = 1-pow(1-.05/3,2);
   for (int i = 1 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double target = omega * pow(.99, i-1) * (1-alpha_i_sq) / denom;
      test_assert(fabs(B->coeff[i]-target) < 1e-9);
   }

   trunc_pol_t *_B = new_trunc_pol_B(50, 2, YES);
   test_assert_critical(_B != NULL);
   test_assert(_B->mono.deg == 0);
   test_assert(_B->mono.coeff == 0);
   test_assert(_B->coeff[0] == 0);
   const double _omega = .01 * (1-pow(1-.05/3,2));
   for (int i = 1 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double target = _omega * pow(.99, i-1) * (1-alpha_i_sq) / denom; 
      test_assert(fabs(_B->coeff[i]-target) < 1e-9);
   }

   // Test the special case N = 0.
   trunc_pol_t *B0 = new_trunc_pol_B(50, 0, NO);
   test_assert_critical(B0 != NULL);
   test_assert(B0->mono.deg == 0);
   test_assert(B0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(B0->coeff[i] == 0);
   }

   trunc_pol_t *_B0 = new_trunc_pol_B(50, 0, YES);
   test_assert_critical(_B0 != NULL);
   test_assert(_B0->mono.deg == 0);
   test_assert(_B0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(_B0->coeff[i] == 0);
   }

   free(B);
   free(_B);
   free(B0);
   free(_B0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_B
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *B;

   redirect_stderr();
   B = new_trunc_pol_B(0, 0, NO);
   unredirect_stderr();

   test_assert(B == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   B = new_trunc_pol_B(0, 2, NO);
   unredirect_stderr();

   test_assert(B == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_B(51, 2, NO);
   unredirect_stderr();

   test_assert(B == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   B = new_trunc_pol_B(18, 2, NO);
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

   // Test a C polynomial of degree 10 with N = 2.
   trunc_pol_t *C = new_trunc_pol_C(10, 2, NO);
   test_assert_critical(C != NULL);
   test_assert(C->mono.deg == 0);
   test_assert(C->mono.coeff == 0);
   test_assert(C->coeff[0] == 0);
   const size_t j = 7; // G-10
   const double omega = .01 * pow(1-.05/3,2);
   const double denom = pow(1-pow(1-.05,j)*.05/3,2) - \
      pow(1-pow(1-.05,j-1)*.05/3,2) - pow(1-pow(1-.05,j),2) + \
      pow(1-(1-.05+.05*.05/3)*pow(1-.05,j-1),2);
   for (int i = 1 ; i <= 10 ; i++) {
      double num =  pow(1-pow(1-.05,j)*.05/3,2) - \
         pow(1-pow(1-.05,j-1)*.05/3,2) - \
         pow(1-pow(1-.05,j)*.05/3-pow(1-.05,i+j-1)*(1-.05/3),2) + \
         pow(1-pow(1-.05,j-1)*.05/3-pow(1-.05,i+j-1)*(1-.05/3),2);
      double target = omega * pow(.99, i-1) * num / denom;
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(C->coeff[i] == 0);
   }

   trunc_pol_t *_C = new_trunc_pol_C(10, 2, YES);
   test_assert_critical(_C != NULL);
   test_assert(_C->mono.deg == 0);
   test_assert(_C->mono.coeff == 0);
   test_assert(_C->coeff[0] == 0);
   const double _omega = .01 * (1 - pow(1-.05/3,2));
   for (int i = 1 ; i <= 10 ; i++) {
      double num =  pow(1-pow(1-.05,j)*.05/3,2) - \
         pow(1-pow(1-.05,j-1)*.05/3,2) - \
         pow(1-pow(1-.05,j)*.05/3-pow(1-.05,i+j-1)*(1-.05/3),2) + \
         pow(1-pow(1-.05,j-1)*.05/3-pow(1-.05,i+j-1)*(1-.05/3),2);
      double target = _omega * pow(.99, i-1) * num / denom;
      test_assert(fabs(_C->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(_C->coeff[i] == 0);
   }

   // Test the special cases N = 0 and N = 1
   trunc_pol_t *C0 = new_trunc_pol_C(50, 0, NO);
   trunc_pol_t *_C0 = new_trunc_pol_C(50, 0, YES);
   trunc_pol_t *C1 = new_trunc_pol_C(50, 1, NO);
   trunc_pol_t *_C1 = new_trunc_pol_C(50, 1, YES);
   test_assert_critical(C0 != NULL);
   test_assert_critical(_C0 != NULL);
   test_assert_critical(C1 != NULL);
   test_assert_critical(_C0 != NULL);
   test_assert(C0->mono.deg == 0);
   test_assert(C0->mono.coeff == 0);
   test_assert(_C0->mono.deg == 0);
   test_assert(_C0->mono.coeff == 0);
   test_assert(C1->mono.deg == 0);
   test_assert(C1->mono.coeff == 0);
   test_assert(_C1->mono.deg == 0);
   test_assert(_C1->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(C0->coeff[i] == 0);
      test_assert(_C0->coeff[i] == 0);
      test_assert(C1->coeff[i] == 0);
      test_assert(_C1->coeff[i] == 0);
   }

   free(C);
   free(_C);
   free(C0);
   free(_C0);
   free(C1);
   free(_C1);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_C
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *C;

   redirect_stderr();
   C = new_trunc_pol_C(0, 0, NO);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   C = new_trunc_pol_C(0, 2, NO);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_C(51, 2, NO);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   C = new_trunc_pol_C(18, 2, NO);
   unredirect_stderr();
   reset_alloc();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_D
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a D polynomial of degree 10 with N = 2.
   trunc_pol_t *D = new_trunc_pol_D(10, 2, NO);
   test_assert_critical(D != NULL);
   test_assert(D->mono.deg == 0);
   test_assert(D->mono.coeff == 0);
   test_assert(D->coeff[0] == 0);
   const double omega = .01 * pow(1-.05/3,2);
   for (int i = 1 ; i <= 10 ; i++) {
      double target = omega * pow(.99, i-1);
      test_assert(fabs(D->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(D->coeff[i] == 0);
   }

   trunc_pol_t *_D = new_trunc_pol_D(10, 2, YES);
   test_assert_critical(_D != NULL);
   test_assert(_D->mono.deg == 0);
   test_assert(_D->mono.coeff == 0);
   test_assert(_D->coeff[0] == 0);
   const double _omega = .01 * (1 - pow(1-.05/3,2));
   for (int i = 1 ; i <= 10 ; i++) {
      double target = _omega * pow(.99, i-1);
      test_assert(fabs(_D->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(_D->coeff[i] == 0);
   }

   // Test the special cases N = 0.
   trunc_pol_t *D0 = new_trunc_pol_D(50, 0, NO);
   trunc_pol_t *_D0 = new_trunc_pol_D(50, 0, YES);
   test_assert_critical(D0 != NULL);
   test_assert_critical(_D0 != NULL);
   test_assert(D0->mono.deg == 0);
   test_assert(D0->mono.coeff == 0);
   test_assert(_D0->mono.deg == 0);
   test_assert(_D0->mono.coeff == 0);

   test_assert(D0->coeff[0] == 0);
   test_assert(_D0->coeff[0] == 0);
   for (int i = 1 ; i <= 50 ; i++) {
      double target = pow(1-.01,i-1) * 0.01;
      test_assert(fabs(D0->coeff[i]-target) < 1e-9);
      test_assert(_D0->coeff[i] == 0);
   }

   free(D);
   free(_D);
   free(D0);
   free(_D0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_D
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *D;

   redirect_stderr();
   D = new_trunc_pol_D(0, 0, NO);
   unredirect_stderr();

   test_assert(D == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   D = new_trunc_pol_D(0, 2, NO);
   unredirect_stderr();

   test_assert(D == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_D(51, 2, NO);
   unredirect_stderr();

   test_assert(D == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   D = new_trunc_pol_D(18, 2, NO);
   unredirect_stderr();
   reset_alloc();

   test_assert(D == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

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
test_new_trunc_pol_w
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a w polynomial of degree 10 with N = 2.
   trunc_pol_t *w = new_trunc_pol_w(10, 2);
   test_assert_critical(w != NULL);

   const double denom = 1-pow(1-.05/3,2);
   const double num = pow(1-pow(1-.05,10)*.05/3,2) - \
      pow(1-pow(1-.05,9)*.05/3,2) - pow(1-pow(1-.05,10),2) + \
      pow(1-(1-.05+.05*.05/3)*pow(1-.05,9),2);
   const double target = pow(1-.01,10) * num / denom;
   test_assert(w->mono.deg == 10);
   test_assert(fabs(w->mono.coeff-target) < 1e-9);
   for (int i = 0 ; i <= 50 ; i++) {
      if (i == 10)
         test_assert(fabs(w->coeff[i]-target) < 1e-9);
      else
         test_assert(w->coeff[i] == 0);
   }

   // Test the special case N = 0.
   trunc_pol_t *w0 = new_trunc_pol_w(1, 0);
   test_assert_critical(w0 != NULL);
   test_assert(w0->mono.deg == 0);
   test_assert(w0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(w0->coeff[i] == 0);
   }
   
   free(w);
   free(w0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_w
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *w;

   redirect_stderr();
   w = new_trunc_pol_w(0, 0);
   unredirect_stderr();

   test_assert(w == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   w = new_trunc_pol_w(0, 2);
   unredirect_stderr();

   test_assert(w == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_w(51, 2);
   unredirect_stderr();

   test_assert(w == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   w = new_trunc_pol_w(10, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(w == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_v
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a v polynomial of degree 10 with N = 2.
   trunc_pol_t *v = new_trunc_pol_v(10, 2);
   test_assert_critical(v != NULL);

   const double num = pow(1-pow(1-.05,10),2) - \
      pow(1-(1-.05+.05*.05/3)*pow(1-.05,9),2);
   const double denom = 1-pow(1-.05/3,2);
   const double target = num / denom * pow(1-.01,10);
   test_assert(v->mono.deg == 10);
   test_assert(fabs(v->mono.coeff-target) < 1e-9);
   for (int i = 0 ; i <= 50 ; i++) {
      if (i == 10)
         test_assert(fabs(v->coeff[i]-target) < 1e-9);
      else
         test_assert(v->coeff[i] == 0);
   }

   // Test the special case N = 0.
   trunc_pol_t *v0 = new_trunc_pol_v(1, 0);
   test_assert_critical(v0 != NULL);
   test_assert(v0->mono.deg == 0);
   test_assert(v0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(v0->coeff[i] == 0);
   }
   
   free(v);
   free(v0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_v
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *v;

   redirect_stderr();
   v = new_trunc_pol_v(0, 0);
   unredirect_stderr();

   test_assert(v == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   v = new_trunc_pol_v(0, 2);
   unredirect_stderr();

   test_assert(v == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_v(51, 2);
   unredirect_stderr();

   test_assert(v == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   v = new_trunc_pol_v(10, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(v == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_y
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a y polynomial of degree 5 with N = 2.
   trunc_pol_t *y = new_trunc_pol_y(4, 5, 2);
   test_assert_critical(y != NULL);

   const double denom =
      pow(1-pow(1-.05,4)*.05/3,2) - \
      pow(1-pow(1-.05,3)*.05/3,2) - \
      pow(1-pow(1-.05,4),2) + \
      pow(1-(1-.05+.05*.05/3)*pow(1-.05,3),2);
   const double num = \
      pow(1-pow(1-.05,4)*.05/3-pow(1-.05,9)*(1-.05/3),2) - \
      pow(1-pow(1-.05,4)*.05/3-pow(1-.05,8)*(1-.05/3),2) - \
      pow(1-pow(1-.05,3)*.05/3-pow(1-.05,9)*(1-.05/3),2) + \
      pow(1-pow(1-.05,3)*.05/3-pow(1-.05,8)*(1-.05/3),2);
   const double target = num / denom * pow(1-.01,5);
   test_assert(y->mono.deg == 5);
   test_assert(fabs(y->mono.coeff-target) < 1e-9);
   for (int i = 0 ; i <= 50 ; i++) {
      if (i == 5)
         test_assert(fabs(y->coeff[i]-target) < 1e-9);
      else
         test_assert(y->coeff[i] == 0);
   }

   // Test the special cases N = 0 and N = 1.
   trunc_pol_t *y0 = new_trunc_pol_y(4, 5, 0);
   test_assert_critical(y0 != NULL);
   test_assert(y0->mono.deg == 0);
   test_assert(y0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(y0->coeff[i] == 0);
   }
   trunc_pol_t *y1 = new_trunc_pol_y(4, 5, 1);
   test_assert_critical(y1 != NULL);
   test_assert(y1->mono.deg == 0);
   test_assert(y1->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(y1->coeff[i] == 0);
   }
   
   free(y);
   free(y0);
   free(y1);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_y
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *y;

   redirect_stderr();
   y = new_trunc_pol_y(10, 0, 2);
   unredirect_stderr();

   test_assert(y == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   y = new_trunc_pol_y(16, 0, 2);
   unredirect_stderr();

   test_assert(y == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   y = new_trunc_pol_y(0, 10, 2);
   unredirect_stderr();

   test_assert(y == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   y = new_trunc_pol_y(0, 16, 2);
   unredirect_stderr();

   test_assert(y == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   y = new_trunc_pol_y(10, 10, 2);
   unredirect_stderr();

   test_assert(y == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   y = new_trunc_pol_y(4, 5, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(y == NULL);
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
   for (int i = 1 ; i <= 50 ; i++) {
      double target = (1-pow(1-pow(1-.05,i)*.05/3,2)) / \
               denom * pow(1-.01,i);
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

   // Test a T_double_down polynomial with N = 2.
   trunc_pol_t *T_ddown = new_trunc_pol_T_double_down(2);
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
   trunc_pol_t *T_ddown0 = new_trunc_pol_T_double_down(0);
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
   T_ddown = new_trunc_pol_T_double_down(2);
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
test_new_trunc_pol_T_sim
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a T_sim polynomial of degree 10 with N = 2.
   trunc_pol_t *T_sim = new_trunc_pol_T_sim(10, 2);
   test_assert_critical(T_sim != NULL);

   test_assert(T_sim->mono.deg == 0);
   test_assert(T_sim->mono.coeff == 0);
   // j = G-10-1 = 6
   const double denom = \
      pow(1-pow(1-.05,6)*.05/3,2) - \
      pow(1-pow(1-.05,5)*.05/3,2) - \
      pow(1-pow(1-.05,6),2) + \
      pow(1-(1-.05+.05*.05/3)*pow(1-.05,5),2);
   test_assert(T_sim->coeff[0] == 1);
   for (int i = 1 ; i <= 10 ; i++) {
      double num = \
         pow(1-pow(1-.05,6)*.05/3,2) - \
         pow(1-pow(1-.05,5)*.05/3,2) - \
         pow(1-pow(1-.05,6)*.05/3-pow(1-.05,6+i)*(1-.05/3),2) + \
         pow(1-pow(1-.05,5)*.05/3-pow(1-.05,6+i)*(1-.05/3),2);
      double target = num / denom * pow(1-.01,i);
      test_assert(fabs(T_sim->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(T_sim->coeff[i] == 0);
   }

   // Test the special cases N = 0.
   trunc_pol_t *T_sim0 = new_trunc_pol_T_sim(10, 0);
   test_assert_critical(T_sim0 != NULL);
   test_assert(T_sim0->mono.deg == 0);
   test_assert(T_sim0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(T_sim0->coeff[i] == 0);
   }
   
   free(T_sim);
   free(T_sim0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_T_sim
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *T_sim;

   redirect_stderr();
   T_sim = new_trunc_pol_T_sim(17, 2);
   unredirect_stderr();

   test_assert(T_sim == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   T_sim = new_trunc_pol_T_sim(10, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(T_sim == NULL);
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

   int success = set_params_mem_prob(17, 100, 0.01, 0.05);
   test_assert_critical(success);

   // Test martrix M with one duplicate because the
   // polynomials are particularly simple in this case.
   matrix_t *M = new_matrix_M(1);
   test_assert_critical(M != NULL);
   
   const int dim = 2*17+1;

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
      double target = pow(.99*.95,G-j-1) * .99*.05;
      test_assert(nspct->mono.deg == G-j);
      test_assert(fabs(nspct->mono.coeff-target) < 1e-9);
      for (int i = 0 ; i <= 100 ; i++) {
         if (i == G-j)
            test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
         else
            test_assert(nspct->coeff[i] == 0);
      }
   }

   // End of second row.
   for (int j = 1 ; j <= G-1 ; j++) {
      test_assert(M->term[dim+G+1+j] == NULL);
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

   // Third row, v terms.
   for (int j = 1 ; j <= G-1 ; j++) {
      test_assert_critical(M->term[2*dim+2+j] != NULL);
      nspct = M->term[2*dim+2+j];
      double target = pow(.99*.95,G-j-1) * .99*.05;
      test_assert(nspct->mono.deg == G-j);
      test_assert(fabs(nspct->mono.coeff-target) < 1e-9);
      for (int i = 0 ; i <= 100 ; i++) {
         if (i == G-j)
            test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
         else
            test_assert(nspct->coeff[i] == 0);
      }
   }

   // Third row, w terms.
   for (int j = 1 ; j <= G-1 ; j++) {
      test_assert_critical(M->term[2*dim+G+1+j] != NULL);
      nspct = M->term[2*dim+G+1+j];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      for (int i = 0 ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }
   }

   // First middle series of rows.
   for (int j = 1 ; j <= G-1 ; j++) {
      // T polynomials.
      test_assert_critical(M->term[(j+2)*dim] != NULL);
      nspct = M->term[(j+2)*dim];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      for (int i = 0 ; i <= j-1 ; i++) {
         double target = pow(.99,i);
         test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
      }
      for (int i = j ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // D polynomials.
      test_assert_critical(M->term[(j+2)*dim+1] != NULL);
      nspct = M->term[(j+2)*dim+1];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      test_assert(nspct->coeff[0] == 0);
      for (int i = 1 ; i <= j ; i++) {
         double target = pow(.99,i-1) * .01*(1-.05/3);
         test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
      }
      for (int i = j+1 ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // Tilde D polynomials.
      test_assert_critical(M->term[(j+2)*dim+2] != NULL);
      nspct = M->term[(j+2)*dim+2];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      test_assert(nspct->coeff[0] == 0);
      for (int i = 1 ; i <= j ; i++) {
         double target = pow(.99,i-1) * .01*.05/3;
         test_assert(fabs(nspct->coeff[i]-target) < 1e-9);
      }
      for (int i = j+1 ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // Rest of the rows.
      for (int i = 1 ; i <= 2*G-2 ; i++) {
         test_assert(M->term[(j+2)*dim+2+i] == NULL);
      }

   }

   // Second middle series of rows.
   for (int j = 1 ; j <= G-1 ; j++) {
      // T polynomials.
      test_assert_critical(M->term[(j+G+1)*dim] != NULL);
      nspct = M->term[(j+G+1)*dim];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      for (int i = 0 ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // C polynomials.
      test_assert_critical(M->term[(j+G+1)*dim+1] != NULL);
      nspct = M->term[(j+G+1)*dim+1];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      for (int i = 0 ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // Tilde C polynomials.
      test_assert_critical(M->term[(j+G+1)*dim+2] != NULL);
      nspct = M->term[(j+G+1)*dim+2];
      test_assert(nspct->mono.deg == 0);
      test_assert(nspct->mono.coeff == 0);
      for (int i = 0 ; i <= 100 ; i++) {
         test_assert(nspct->coeff[i] == 0);
      }

      // Row of y polynomials.
      for (int i = 1 ; i < j ; i++) {
         test_assert_critical(M->term[(j+G+1)*dim+2+i] != NULL);
         nspct = M->term[(j+G+1)*dim+2+i];
         test_assert(nspct->mono.deg == 0);
         test_assert(nspct->mono.coeff == 0);
         for (int i = 0 ; i <= 100 ; i++) {
            test_assert(nspct->coeff[i] == 0);
         }
      }

      // Rest of the rows.
      for (int i = j ; i <= 2*G-2 ; i++) {
         test_assert(M->term[(j+G+1)*dim+2+i] == NULL);
      }

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

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // The first terms can be computed directly.
   for (int i = 0 ; i < 17 ; i++) {
      test_assert(fabs(compute_mem_prob(2,i)-1) < 1e-9);
   }
   double target_17 = 1-pow(.99,17);
   test_assert(fabs(compute_mem_prob(2,17)-target_17) < 1e-9);

   double target_18;
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,2);
   test_assert(fabs(compute_mem_prob(2,18)-target_18) < 1e-9);

   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,3);
   test_assert(fabs(compute_mem_prob(3,18)-target_18) < 1e-9);

   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,4);
   test_assert(fabs(compute_mem_prob(4,18)-target_18) < 1e-9);


   success = set_params_mem_prob(20, 50, 0.02, 0.05);
   test_assert_critical(success);

   // The first terms can be computed directly.
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
   {"new_trunc_pol_D",             test_new_trunc_pol_D},
   {"error_new_trunc_pol_D",       test_error_new_trunc_pol_D},
   {"new_trunc_pol_u",             test_new_trunc_pol_u},
   {"error_new_trunc_pol_u",       test_error_new_trunc_pol_u},
   {"new_trunc_pol_v",             test_new_trunc_pol_v},
   {"error_new_trunc_pol_v",       test_error_new_trunc_pol_v},
   {"new_trunc_pol_w",             test_new_trunc_pol_w},
   {"error_new_trunc_pol_w",       test_error_new_trunc_pol_w},
   {"new_trunc_pol_y",             test_new_trunc_pol_y},
   {"error_new_trunc_pol_y",       test_error_new_trunc_pol_y},
   {"new_trunc_pol_T_down",        test_new_trunc_pol_T_down},
   {"error_new_trunc_pol_T_down",  test_error_new_trunc_pol_T_down},
   {"new_trunc_pol_T_ddown",       test_new_trunc_pol_T_ddown},
   {"error_new_trunc_pol_T_ddown", test_error_new_trunc_pol_T_ddown},
   {"new_trunc_pol_T_up",          test_new_trunc_pol_T_up},
   {"error_new_trunc_pol_T_up",    test_error_new_trunc_pol_T_up},
   {"new_trunc_pol_T_sim",         test_new_trunc_pol_T_sim},
   {"error_new_trunc_pol_T_sim",   test_error_new_trunc_pol_T_sim},
   {"new_null_matrix",             test_new_null_matrix},
   {"error_new_null_matrix",       test_error_new_null_matrix},
   {"new_zero_matrix",             test_new_zero_matrix},
   {"error_new_zero_matrix",       test_error_new_zero_matrix},
   {"new_matrix_M",                test_new_matrix_M},
   {"error_new_matrix_M",          test_error_new_matrix_M},
   {"compute_mem_prob",            test_compute_mem_prob},
   {"error_compute_mem_prob",      test_error_compute_mem_prob},
   {NULL, NULL},
};
