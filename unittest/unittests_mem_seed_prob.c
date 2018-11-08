#include "unittest.h"
#include "mem_seed_prob.c"

#if 0

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
   test_assert_stderr("[mem_seed_prob] error in function `fault_");

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   redirect_stderr();
   // The error is that 'N' is greater than 'MAXN'.
   x = mem_seed_prob(2,1025);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `fault_");

   redirect_stderr();
   // The error is that 'k' is greater than specified value above.
   x = mem_seed_prob(51,2);
   unredirect_stderr();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `fault_");

   set_alloc_failure_countdown_to(0);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   x = mem_seed_prob(10,2);
   unredirect_stderr();
   reset_alloc();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail (somewhere else).
   x = mem_seed_prob(10,2);
   unredirect_stderr();
   reset_alloc();

   test_assert(x != x);
   test_assert_stderr("[mem_seed_prob] error in function `new_n");
   clean_mem_prob();

}


void
test_average_errors
(void)
{

   double array[] = { 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
      0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17,
      0.163141136132321, 0.156113683625966, 0.148917642480934,
      0.141553012697225, 0.13401979427484, 0.126317987213777,
      0.118447591514038, 0.110408607175622, 0.102201034198529,
      0.0938248725827594, 0.085280122328313, 0.0765667834351898,
      0.0676848559033898, 0.0586343397329131, 0.0494152349237596,
      0.0400275414759292, 0.0304712593894221, 0.0278519209369612,
      0.0253482151367323, 0.0229622736484174, 0.0206962281316983,
      0.0185522102462568, 0.0165323516517746, 0.0146387840079337,
      0.0128736389744158, 0.0112390482109028, 0.00973714337707639,
      0.0083700561326185, 0.0071399181372109, 0.00604886105053543,
      0.00509901653227388, 0.00429251624210809, 0.00363149183971987,
      0.00311807498479103, 0.00275439733700339, 0.00242279935474554,
      0.0021218189616605, 0.00184997012315102, 0.00160574284637961,
      0.00138760318026852, 0.00119399321549972, 0.00102333108451495,
      0.00087401096151568, 0.000744403062463109, 0.000632853645078193,
      0.000537685008841622, 0.000457195494993832, 0.000389659486534996,
      0.000333327408225031, 0.000286425726583594, 0.000247156949890084,
      0.000213699628183641, 0.000184208353263147, 0.000158328416353586,
      0.000135721541310787, 0.000116066137064366, 9.90575500606736e-05,
      8.44083167057381e-05, 7.18484158082088e-05, 6.11255210223017e-05,
      5.20052532907434e-05, 4.42714332877157e-05, 3.77263338617996e-05,
      3.21909324789203e-05, 2.7505163665291e-05, 2.35281714503576e-05,
      2.01385618097431e-05, 1.7234655108192e-05, 1.47347385425144e-05,
      1.25773185845308e-05, 1.07213734240163e-05, 9.12958180670945e-06,
      7.76814988713323e-06, 6.60663552787539e-06, 5.61777004532791e-06,
      4.77727740188584e-06, 4.06369084460529e-06, 3.45816699032078e-06,
      2.94429735722175e-06, 2.50791734288845e-06, 2.13691264878695e-06,
      1.82102315122354e-06, 1.55164421875829e-06 };
   
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(fabs(average_errors(17, i ,.01)-array[i]) < 1e-9);
   }

}
#endif


void
test_omega
(void)
{

   // Simple tests for small values.
   test_assert(fabs(omega(0,.05,0)-1.0) < 1e-9);
   test_assert(fabs(omega(1,.05,0)-0.0) < 1e-9);
   test_assert(fabs(omega(2,.05,0)-0.0) < 1e-9);
   test_assert(fabs(omega(0,.05,1)-(1-.05/3)) < 1e-9);
   test_assert(fabs(omega(1,.05,1)-.05/3) < 1e-9);
   test_assert(fabs(omega(2,.05,1)-0.0) < 1e-9);
   test_assert(fabs(omega(0,.05,2)-(1-.05/3)*(1-.05/3)) < 1e-9);
   test_assert(fabs(omega(1,.05,2)-2*(1-.05/3)*.05/3) < 1e-9);
   test_assert(fabs(omega(2,.05,2)-.05/3*.05/3) < 1e-9);

   // Tests for high values checked with R.
   test_assert(fabs(omega(0,.05,9)-0.8596206731) < 1e-9);
   test_assert(fabs(omega(1,.05,9)-0.1311285773) < 1e-9);
   test_assert(fabs(omega(5,.05,9)-1.515016e-07) < 1e-9);
   test_assert(fabs(omega(0,.05,23)-0.6793874326) < 1e-9);
   test_assert(fabs(omega(1,.05,23)-0.2648459483) < 1e-9);
   test_assert(fabs(omega(5,.05,23)-3.197640e-05) < 1e-9);
   test_assert(fabs(omega(0,.05,58)-0.3772629470) < 1e-9);
   test_assert(fabs(omega(1,.05,58)-0.3708686598) < 1e-9);
   test_assert(fabs(omega(5,.05,58)-0.0024179659) < 1e-9);

   // Tests for small values checked with R.
   // Here we use relative error because the values are small.
   test_assert(fabs(omega(9,.05,9)/9.92290301275217213e-17 - 1) < 1e-6);
   test_assert(fabs(omega(18,.05,23)/3.0461652643903103e-28 - 1) < 1e-6);
   test_assert(fabs(omega(23,.05,23)/1.2662551980515106e-41 - 1) < 1e-6);
   test_assert(fabs(omega(18,.05,58)/2.2620043844369284e-18 - 1) < 1e-6);
   test_assert(fabs(omega(34,.05,58)/2.9921066568544365e-45 - 1) < 1e-6);
   test_assert(fabs(omega(44,.05,58)/4.6272133972034210e-66 - 1) < 1e-6);
   test_assert(fabs(omega(58,.05,58)/7.3659281411611084e-104 - 1) < 1e-6);

}


void
test_psi
(void)
{

   test_assert(fabs(psi(1,0,0,0,.05,1)-1.0) < 1e-9);
   test_assert(fabs(psi(17,0,1,1,.05,1)-1.0) < 1e-9);
   test_assert(fabs(psi(15,1,3,7,.05,12)-0.0) < 1e-9);
   test_assert(fabs(psi(15,7,1,3,.05,12)-0.0) < 1e-9);
   test_assert(fabs(psi(30,14,19,15,.05,23)-0.0) < 1e-9);

   // Tests for values checked with R.
   test_assert(fabs(psi(17,1,3,2,.05,5)-2.1637593296) < 1e-9);
   test_assert(fabs(psi(15,1,7,3,.05,12)-14.466081668) < 1e-9);

}


void
test_zeta
(void)
{

   // Simple test for small values.
   test_assert(fabs(zeta(20,0,0,.05,1)-0.0) < 1e-9);
   test_assert(fabs(zeta(17,0,1,.05,1)-pow(.95,17)) < 1e-9);

   // Tests for values checked with R.
   test_assert(fabs(zeta(20,0,1,.05,1)-0.35848592241) < 1e-9);
   test_assert(fabs(zeta(20,4,5,.05,8)-0.65255894306) < 1e-9);
   test_assert(fabs(zeta(20,6,3,.05,8)-0.25509528217) < 1e-9);
   test_assert(fabs(zeta(30,14,19,.05,23)-0.83001525276) < 1e-9);

}

void
test_set_params_mem_prob
(void)
{

   int success;

   // The parameters are 'G' (gamma) the minimum seed size,
   // 'K' the maximum size of the reads, 'P' the probability of
   // read error and 'U' (mu) the divergence rate between
   // the duplicates and the target.
   success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   test_assert(G == 17);
   test_assert(K == 50);
   test_assert(P == 0.01);

   test_assert(KSZ == sizeof(trunc_pol_t) + 51*sizeof(double));

   for (int i = 0 ; i < HSIZE ; i++) {
      test_assert(MEMF[i] == NULL);
      test_assert(TYPI[i] == NULL);
   }

   test_assert_critical(TEMP != NULL);
   test_assert(TEMP->monodeg == 0);
   for (int i = 0 ; i <= K ; i++) {
      test_assert(TEMP->coeff[i] == 0);
   }

   test_assert(PARAMS_INITIALIZED);

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
   success = set_params_mem_prob(17, 50, 0.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   // Case 2.
   redirect_stderr();
   // The error is that 'p' is 1.0.
   success = set_params_mem_prob(17, 50, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   // Case 3.
   redirect_stderr();
   // The error is that 'G' (gamma) is 0.
   success = set_params_mem_prob(0, 50, 0.01);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   // Case 4.
   redirect_stderr();
   // The error is that 'K' is 0.
   success = set_params_mem_prob(17, 0, 0.01);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[mem_seed_prob] error in function `set_");

   // Casae 5. Test memory error.
   set_alloc_failure_rate_to(1.0);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   success = set_params_mem_prob(17, 50, 0.01);
   unredirect_stderr();
   reset_alloc();
   test_assert_stderr("[mem_seed_prob] error in function `new_zero");

   return;

}


void
test_uninitialized_error
(void)
{

   // Do not call 'set_params_mem_prob()'.
   redirect_stderr();
   double x = prob_MEM_failure(5, .05, 20);
   unredirect_stderr();
   test_assert_stderr("[mem_seed_prob] error in function `dynamic_p");
   test_assert(x != x);

   return;

}

void
test_new_zero_trunc_pol
(void)
{

   size_t k = 50;

   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   trunc_pol_t *a = new_zero_trunc_pol();
   test_assert_critical(a);

   test_assert(a->monodeg == 0);
   for (int i = 0 ; i <= k ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   free(a);
   clean_mem_prob();

   return;

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
test_trunc_pol_update_add
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01);
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
test_trunc_pol_mult
(void)
{
   
   // Test multiplications between zero polynomials.

   size_t k = 50;

   set_params_mem_prob(17, k, 0.01);
   trunc_pol_t *a = new_zero_trunc_pol();

   test_assert_critical(a != NULL);

   test_assert(trunc_pol_mult(a, NULL, NULL) == a);

   test_assert_critical(a != NULL);
   test_assert(a->monodeg == 0);
   for (int i = 0 ; i <= k ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   trunc_pol_t *b = new_zero_trunc_pol();
   test_assert_critical(b != NULL);

   test_assert(trunc_pol_mult(a, b, NULL) == a);

   // Same remark as above.
   test_assert_critical(a != NULL);
   test_assert(a->monodeg == 0);
   for (int i = 0 ; i <= k ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   test_assert(trunc_pol_mult(a, NULL, b) == a);

   // Same remark as above.
   test_assert_critical(a != NULL);
   test_assert(a->monodeg == 0);
   for (int i = 0 ; i <= k ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   trunc_pol_t *c = new_zero_trunc_pol();
   test_assert_critical(c != NULL);

   test_assert(trunc_pol_mult(a, b, c) == a);

   // Same remark as above.
   test_assert_critical(a != NULL);
   test_assert(a->monodeg == 0);
   for (int i = 0 ; i <= k ; i++) {
      test_assert(a->coeff[i] == 0);
   }


   // Test multiplications between monomials (b = 5z and c = z^2).
   b->monodeg = 1; b->coeff[1] = 5;
   c->monodeg = 2; c->coeff[2] = 1;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->monodeg == 3);
   for (int i = 0 ; i <= k ; i++) {
      double target = i == 3 ? 5 : 0;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->monodeg == 3);
   for (int i = 0 ; i <= k ; i++) {
      double target = i == 3 ? 5 : 0;
      test_assert(a->coeff[i] == target);
   }

   // Test overflow (multiplying two monomials of sufficient
   // degree yields a zero polynomial).
   bzero(b, KSZ);
   bzero(c, KSZ);
   b->monodeg = k-1; b->coeff[1] = 0; b->coeff[k-1] = 5;
   c->monodeg = k-1; c->coeff[2] = 0; c->coeff[k-1] = 1;

   test_assert(trunc_pol_mult(a, b, c) == NULL);
   test_assert(a != NULL);


   // Test multiplications between a monomial and a
   // polynomial (b = 5z and c = z^2 + 2z^3).
   bzero(b, KSZ);
   bzero(c, KSZ);
   b->monodeg = 1; b->coeff[1] = 5;
   c->monodeg = k+1; c->coeff[2] = 1; c->coeff[3] = 2;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->monodeg > k);
   for (int i = 0 ; i <= k ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 10;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->monodeg > k);
   for (int i = 0 ; i <= k ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 10;
      test_assert(a->coeff[i] == target);
   }

   // Test multiplications between two polynomials
   // (b = 5z + 3z^2 and c = z^2 + 2z^3).
   bzero(b, KSZ);
   bzero(c, KSZ);
   b->monodeg = k+1; b->coeff[1] = 5; b->coeff[2] = 3;
   c->monodeg = k+1; c->coeff[2] = 1; c->coeff[3] = 2;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->monodeg > k);
   for (int i = 0 ; i <= k ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 13;
      if (i == 5) target = 6;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->monodeg > k);
   for (int i = 0 ; i <= k ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 13;
      if (i == 5) target = 6;
      test_assert(a->coeff[i] == target);
   }

   // Test that higher terms overflow.
   // (b = 5z + 3z^2 + z^49 and c = z^2 + 2z^3 + z^49).
   b->coeff[k-1] = 1;
   c->coeff[k-1] = 1;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->monodeg > k);
   for (int i = 0 ; i <= k ; i++) {
      double target = 0;
      if (i == 3)  target = 5;
      if (i == 4)  target = 13;
      if (i == 5)  target = 6;
      if (i == 50) target = 5;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->monodeg > k);
   for (int i = 0 ; i <= k ; i++) {
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
test_new_trunc_pol_A
(void)
{

   size_t k = 50;
   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   trunc_pol_t *A;
   double factor;


   A = new_trunc_pol_A(0,0,.05,1);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg > k);
   factor = .01 * (1-.05/3);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= k ; i++) {
      test_assert(fabs(A->coeff[i]) < 1e-9);
   }

   free(A);
   A = NULL;

   A = new_trunc_pol_A(1,0,.05,1);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg > k);
   factor = .01 * (1-.05/3);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   // The terms below are exactly the same as above. They
   // are separated only to highlight the two terms of the sum.
   for (int i = 18 ; i <= k ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   free(A);
   A = NULL;

   A = new_trunc_pol_A(0,1,.05,1);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg > k);
   factor = .01 * .05/3;
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= k ; i++) {
      double target = pow(.95 * .99, i-1) * .01 * .05/3;
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   free(A);
   A = NULL;

   A = new_trunc_pol_A(1,1,.05,1);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg > k);
   factor = .01 * .05/3;
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   // The terms below are exactly the same as above. They
   // are separated only to highlight the two terms of the sum.
   for (int i = 18 ; i <= k ; i++) {
      double target = factor * pow(.95 * .99, i-1);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   free(A);
   A = NULL;

   // This one was checked with R.
   A = new_trunc_pol_A(3,5,.05,9);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg > k);
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
   for (int i = 18 ; i <= k ; i++) {
      double xi_term = pow(1-pow(.95, i-1), 3);
      double target = (1 - xi_term * (1 - zeta_terms[i-1])) *
            pow(.99, i-1) * factor;
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   free(A);
   A = NULL;

   // Test special case N = 0.
   A = new_trunc_pol_A(0,0,.05,0);
   test_assert_critical(A != NULL);

   test_assert(A->monodeg == 1);
   for (int i = 0 ; i <= k ; i++) {
      double target = 0;
      if (i == 1) target = .01;
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   free(A);

   clean_mem_prob();

}


void
test_error_new_trunc_pol_A
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   trunc_pol_t *A;


   // Test error for 'new_trunc_pol_A()'.

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   A = new_trunc_pol_A(0,0,.05,1);
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

   size_t k = 50;

   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   trunc_pol_t *B = NULL;

   // Test i = 2, N = 1.
   B = new_trunc_pol_B(2,.05,1);
   test_assert_critical(B != NULL);
   test_assert(B->monodeg == 2);
   for (int i = 0 ; i <= k ; i++) {
      double target = 0.0;
      if (i == 2) target = 0.05*(1-0.05) * .99*.99;
      test_assert(fabs(B->coeff[i]-target) < 1e-9);
   }

   free(B);
   B = NULL;

   // Test i = 3, N = 1.
   B = new_trunc_pol_B(3,.05,1);
   test_assert_critical(B != NULL);
   test_assert(B->monodeg == 3);
   for (int i = 0 ; i <= k ; i++) {
      double target = 0.0;
      if (i == 3) target = 0.05*(1-0.05)*(1-0.05) * .99*.99*.99;
      test_assert(fabs(B->coeff[i]-target) < 1e-9);
   }

   free(B);
   B = NULL;

   // Test the special case N = 0.
   
   B = new_trunc_pol_B(2,.05,0);
   test_assert_critical(B != NULL);
   test_assert(B->monodeg == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(B->coeff[i] == 0);
   }

   free(B);
   B = NULL;

   // Only the case i = 1 is non-zero.
   B = new_trunc_pol_B(1,.05,0);
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

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   trunc_pol_t *B;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   B = new_trunc_pol_B(1,.05,2);
   unredirect_stderr();
   reset_alloc();

   test_assert(B == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   redirect_stderr();
   // The error is that 'i = 0'.
   B = new_trunc_pol_B(0,.05,2);
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
   
   size_t k = 50;

   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   trunc_pol_t *C;
   
   C = new_trunc_pol_C(0,.05,1);
   test_assert_critical(C != NULL);
   test_assert(C->monodeg > k);
   for (int i = 0 ; i <= 16 ; i++) {
      double target = pow((1-.05) * .99, i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= k ; i++) {
      test_assert(C->coeff[i] == 0);
   }

   free(C);
   C = NULL;
   
   C = new_trunc_pol_C(1,.05,1);
   test_assert_critical(C != NULL);
   test_assert(C->monodeg > k);
   for (int i = 0 ; i <= k ; i++) {
      double target = pow((1-.05) * .99, i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }

   free(C);
   C = NULL;
   
   C = new_trunc_pol_C(5,.05,10);
   test_assert_critical(C != NULL);
   test_assert(C->monodeg > k);
   for (int i = 0 ; i <= 16 ; i++) {
      double target = (1-pow(1-pow(.95,i),10)) * pow(.99,i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= k ; i++) {
      double target = (1-pow(1-pow(.95,i),5)) * pow(.99,i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }

   free(C);
   C = NULL;

   // Test the special case that 'N' is 0. Do not test
   // any other value than 'm' = 0 (see comment above).
   
   C = new_trunc_pol_C(0,.05,0);
   test_assert_critical(C != NULL);
   test_assert(C->monodeg == 0);
   test_assert(C->coeff[0] == 1.0);
   for (int i = 1 ; i <= k ; i++) {
      test_assert(C->coeff[i] == 0);
   }

   free(C);

   clean_mem_prob();

}


void
test_error_new_trunc_pol_C
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   trunc_pol_t *C;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   C = new_trunc_pol_C(1,.05,3);
   unredirect_stderr();
   reset_alloc();

   test_assert(C == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_D
(void)
{

   size_t k = 50;
   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   trunc_pol_t *D;

   D = new_trunc_pol_D(1,0,.05,1);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg > k);
   test_assert(D->coeff[0] == 0);
   for (int i = 1 ; i <= 16 ; i++) {
      double target = .01 * (1-.05/3) * pow(.99, i-1);
      test_assert(fabs(D->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= k ; i++) {
      test_assert(D->coeff[i] == 0);
   }
   
   free(D);
   D = NULL;

   D = new_trunc_pol_D(16,0,.05,1);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg == 1);
   for (int i = 0 ; i <= k ; i++) {
      double target = 0;
      if (i == 1) target = .01 * (1-.05/3);
      test_assert(fabs(D->coeff[i]-target) < 1e-9);
   }
   
   free(D);
   D = NULL;

   D = new_trunc_pol_D(13,18,.05,20);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg > k);
   test_assert(D->coeff[0] == 0);
   // This is "p times omega sub m" computed with R.
   const double factor = 1.80897521494885666448e-30;
   for (int i = 1 ; i <= 4 ; i++) {
      double target = factor * .01 * pow(.99, i-1);
      // Use the relative error because the numbers are tiny.
      test_assert(fabs(D->coeff[i]/target - 1) < 1e-6);
   }
   for (int i = 5 ; i <= k ; i++) {
      test_assert(D->coeff[i] == 0);
   }
   
   free(D);
   D = NULL;

   // Test the case that 'm' is greater than 'N'.

   D = new_trunc_pol_D(1,10,.05,5);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg == 0);
   for (int i = 0 ; i <= k ; i++) {
      test_assert(D->coeff[i] == 0);
   }

   free(D);
   D = NULL;

   // Test the special cases that 'N' is 0.

   D = new_trunc_pol_D(1,0,.05,0);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg > k);
   test_assert(D->coeff[0] == 0);
   for (int i = 1 ; i <= 16 ; i++) {
      double target = .01 * pow(.99, i-1);
      test_assert(fabs(D->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= k ; i++) {
      test_assert(D->coeff[i] == 0);
   }
   
   free(D);
   D = NULL;

   D = new_trunc_pol_D(16,0,.05,0);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg == 0);
   for (int i = 0 ; i <= k ; i++) {
      test_assert(D->coeff[i] == 0);
   }
   
   free(D);
   D = NULL;

   D = new_trunc_pol_D(1,1,.05,0);
   test_assert_critical(D != NULL);

   test_assert(D->monodeg == 0);
   for (int i = 0 ; i <= k ; i++) {
      test_assert(D->coeff[i] == 0);
   }
   
   free(D);

   clean_mem_prob();

}


void
test_error_new_trunc_pol_D
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   trunc_pol_t *D;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   D = new_trunc_pol_D(5,1,.05,3);
   unredirect_stderr();
   reset_alloc();

   test_assert(D == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   redirect_stderr();
   // The error is that 'j' is greater than G-1
   D = new_trunc_pol_D(17,1,.05,3);
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

   size_t k = 50;

   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   trunc_pol_t *E;

   E = new_trunc_pol_E(16);
   test_assert_critical(E != NULL);

   // E polynomials with j = G-1 are monomials.
   test_assert(E->monodeg == 0);
   test_assert(E->coeff[0] == 1);
   for (int i = 1 ; i <= k ; i++) {
      test_assert(E->coeff[i] == 0);
   }

   free(E);
   E = NULL;

   E = new_trunc_pol_E(1);
   test_assert_critical(E != NULL);

   test_assert(E->monodeg > k);
   for (int i = 0 ; i <= 15 ; i++) {
      double target = pow(1-P,i);
      test_assert(fabs(E->coeff[i]-target) < 1e-9);
   }
   for (int i = 16 ; i <= k ; i++) {
      test_assert(E->coeff[i] == 0);
   }

   free(E);

   clean_mem_prob();

}


void
test_error_new_trunc_pol_E
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01);
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
test_new_trunc_pol_F
(void)
{

   size_t k = 50;

   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   trunc_pol_t *F;

   F = new_trunc_pol_F(0,.05);
   test_assert_critical(F != NULL);

   // F polynomials with j = 0 are monomials.
   test_assert(F->monodeg == 0);
   test_assert(F->coeff[0] == 1);
   for (int i = 1 ; i <= k ; i++) {
      test_assert(F->coeff[i] == 0);
   }

   free(F);
   F = NULL;

   F = new_trunc_pol_F(16,.05);
   test_assert_critical(F != NULL);

   test_assert(F->monodeg > k);
   for (int i = 0 ; i <= 16 ; i++) {
      double target = pow(.99*.95,i);
      test_assert(fabs(F->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= k ; i++) {
      test_assert(F->coeff[i] == 0);
   }

   free(F);

   clean_mem_prob();

}


void
test_error_new_trunc_pol_F
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   trunc_pol_t *F;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   F = new_trunc_pol_F(2,.05);
   unredirect_stderr();
   reset_alloc();

   test_assert(F == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   redirect_stderr();
   // The error is that 'j' is greater than G-1
   F = new_trunc_pol_F(17,.05);
   unredirect_stderr();

   test_assert(F == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_tr");

   clean_mem_prob();

}


void
test_new_trunc_pol_R
(void)
{

   size_t k = 50;

   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   trunc_pol_t *R;

   R = new_trunc_pol_R(0,.05);
   test_assert_critical(R != NULL);

   // R polynomials with j = 0 are monomials.
   test_assert(R->monodeg == 1);
   test_assert(R->coeff[1] == .01*(1-.05/3));
   test_assert(R->coeff[0] == 0);
   for (int i = 2 ; i <= k ; i++) {
      test_assert(R->coeff[i] == 0);
   }

   free(R);
   R = NULL;

   R = new_trunc_pol_R(16,.05);
   test_assert_critical(R != NULL);

   test_assert(R->monodeg > k);
   test_assert(R->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = pow(.99*.95,i-1) * .01*(1-.05/3);
      test_assert(fabs(R->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= k ; i++) {
      test_assert(R->coeff[i] == 0);
   }

   free(R);

   clean_mem_prob();

}


void
test_error_new_trunc_pol_R
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   trunc_pol_t *R;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   R = new_trunc_pol_R(2,.05);
   unredirect_stderr();
   reset_alloc();

   test_assert(R == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   redirect_stderr();
   // The error is that 'j' is greater than G-1
   R = new_trunc_pol_R(17,.05);
   unredirect_stderr();

   test_assert(R == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_tr");

   clean_mem_prob();

}


void
test_new_trunc_pol_r
(void)
{

   size_t k = 50;

   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   trunc_pol_t *r;

   // Test r+ polynomials.
   for (int i = 0 ; i <= 15 ; i++) {
      r = new_trunc_pol_r_plus(i,.05);
      test_assert_critical(r != NULL);

      test_assert(r->monodeg == i+1);
      for (int j = 0 ; j <= k ; j++) {
         double target = 0;
         if (j == i+1) target = pow(.99*.95, i-1) * .99*.05;
         test_assert((r->coeff[j]-target) < 1e-9);
      }

      free(r);
      r = NULL;
   }

   // Test r- polynomials.
   for (int i = 0 ; i <= 15 ; i++) {
      r = new_trunc_pol_r_minus(i,.05);
      test_assert_critical(r != NULL);

      test_assert(r->monodeg == i+1);
      for (int j = 0 ; j <= k ; j++) {
         double target = 0;
         if (j == i+1) target = pow(.99*.95, i-1) * .01*.05/3;
         test_assert((r->coeff[j]-target) < 1e-9);
      }

      free(r);
      r = NULL;
   }

   clean_mem_prob();

}


void
test_error_new_trunc_pol_r
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   trunc_pol_t *r;

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   r = new_trunc_pol_r_plus(2,.05);
   unredirect_stderr();
   reset_alloc();

   test_assert(r == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   r = new_trunc_pol_r_minus(2,.05);
   unredirect_stderr();
   reset_alloc();

   test_assert(r == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   redirect_stderr();
   // The error is that 'j' is greater than G-1
   r = new_trunc_pol_r_plus(16,.05);
   unredirect_stderr();

   test_assert(r == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_tr");

   redirect_stderr();
   // The error is that 'j' is greater than G-1
   r = new_trunc_pol_r_minus(16,.05);
   unredirect_stderr();

   test_assert(r == NULL);
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

   int success = set_params_mem_prob(17, 50, 0.01);
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

   clean_mem_prob();

}


void
test_error_new_zero_matrix
(void)
{

   matrix_t *matrix;

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

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

   clean_mem_prob();

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

   size_t k = 50;

   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   // Test martrix M with 0 duplicate because it is a special case.
   M = new_matrix_M(.05,0);
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
   for (int i = 1 ; i <= k ; i++) {
      test_assert(C->coeff[i] == 0);
   }

   // -- Next 'G-1' rows -- //
   
   for (int i = 1 ; i <= G-1 ; i++) {

      // First term (D polynomial).
      D = M->term[i*dim0];
      test_assert_critical(D != NULL);
      if (i == 1) {
         test_assert(D->monodeg > k);
         test_assert(D->coeff[0] == 0);
         for (int j = 1 ; j <= G-1 ; j++) {
            double target = pow(.99, j-1) * .01;
            test_assert(fabs(D->coeff[j]-target) < 1e-9);
         }
         for (int j = G ; j <= k ; j++) {
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
         test_assert(E->monodeg > k);
         for (int j = 0 ; j <= G-2 ; j++) {
            double target = pow(.99, j);
            test_assert(fabs(E->coeff[j]-target) < 1e-9);
         }
         for (int j = G-1 ; j <= k ; j++) {
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

   // Test martrix M with 'N' = 1 duplicate because the
   // polynomials are particularly simple in this case.
   M = new_matrix_M(.05,1);
   test_assert_critical(M != NULL);

   const int dim1 = 19; // 17+1+1
   test_assert(M->dim == dim1);
   

   // -- First row -- //

   // First two terms (A polynomials).
   A = M->term[0];
   test_assert_critical(A != NULL);
   test_assert(A->monodeg > k);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = pow(.95 * .99, i-1) * .01 * (1-.05/3);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= k ; i++) {
      test_assert(A->coeff[i] == 0);
   }

   A = M->term[1];
   test_assert_critical(A != NULL);
   test_assert(A->monodeg > k);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= k ; i++) {
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
   test_assert(C->monodeg > k);
   for (int i = 0 ; i <= 16 ; i++) {
      double target = pow(.95 * .99, i);
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }
   for (int i = 17 ; i <= k ; i++) {
      test_assert(C->coeff[i] == 0);
   }

   
   // -- Second row -- //

   // First two terms (A polynomials).
   A = M->term[1*dim1+0];
   test_assert_critical(A != NULL);
   test_assert(A->monodeg > k);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= k ; i++) {
      double target = pow(.95 * .99, i-1) * .01 * (1-.05/3);
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }

   A = M->term[1*dim1+1];
   test_assert_critical(A != NULL);
   test_assert(A->monodeg > k);
   test_assert(A->coeff[0] == 0);
   for (int i = 1 ; i <= k ; i++) {
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
   test_assert(C->monodeg > k);
   for (int i = 0 ; i <= k ; i++) {
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
         test_assert(D->monodeg > k);
      test_assert(D->coeff[0] == 0);
      for (int j = 1 ; j <= G-i ; j++) {
         double target = pow(.99, j-1) * .01 * (1-.05/3);
         test_assert(fabs(D->coeff[j]-target) < 1e-9);
      }
      for (int j = G-i+1 ; j <= k ; j++) {
         test_assert(D->coeff[j] == 0);
      }

      D = M->term[(i+1)*dim1+1];
      test_assert_critical(D != NULL);
      if (i == 16)
         test_assert(D->monodeg == 1);
      else
         test_assert(D->monodeg > k);
      test_assert(D->coeff[0] == 0);
      for (int j = 1 ; j <= G-i ; j++) {
         double target = pow(.99, j-1) * .01 * .05/3;
         test_assert(fabs(D->coeff[j]-target) < 1e-9);
      }
      for (int j = G-i+1 ; j <= k ; j++) {
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
         test_assert(E->monodeg > k);
      for (int j = 0 ; j <= G-i-1 ; j++) {
         double target = pow(.99, j);
         test_assert(fabs(E->coeff[j]-target) < 1e-9);
      }
      for (int j = G-i+1 ; j <= k ; j++) {
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

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   matrix_t *M;

   set_alloc_failure_countdown_to(0);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   M = new_matrix_M(.05,1);
   unredirect_stderr();
   reset_alloc();

   test_assert(M == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_n");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail (later).
   M = new_matrix_M(.05,1);
   unredirect_stderr();
   reset_alloc();

   test_assert(M == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_matrix_L
(void)
{

   size_t k = 50;

   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   matrix_t *L;

   trunc_pol_t *F;
   trunc_pol_t *R;
   trunc_pol_t *r;

   // Test martrix M with 0 duplicate because it is a special case.
   L = new_matrix_L(.05);
   test_assert_critical(L != NULL);

   const int dim0 = 34; // 2 x 17
   test_assert(L->dim == dim0);
   
   const double a = .99 * .95;            // (1-p) * (1-u)
   const double b = .01 * .05/3;          // p * u/3;
   const double c = .99 * .05;            // (1-p) * u;
   const double d = .01 * (1-.05/3);      // p * (1-u/3);

   // -- First row -- //

   // First term (R polynomial).
   R = L->term[0];
   test_assert_critical(R != NULL);
   test_assert(R->monodeg > k);
   test_assert(R->coeff[0] == 0);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = pow(a,i-1) * d;
      test_assert(fabs(R->coeff[i]-target) < 1e-9);
   } 
   for (int i = 18 ; i <= k ; i++) {
      test_assert(R->coeff[i] == 0);
   }

   // Next 16 terms are r+ polynomials.
   for (int j = 1 ; j <= 16 ; j++) {
      r = L->term[j];
      test_assert_critical(L != NULL);
      test_assert(r->monodeg == j);
      for (int i = 0 ; i <= k ; i++) {
         double target = 0;
         if (i == j) target = pow(a,i-1) * c;
         test_assert(fabs(r->coeff[i]-target) < 1e-9);
      }
   }
   // Next 16 terms are r- polynomials.
   for (int j = 1 ; j <= 16 ; j++) {
      r = L->term[j+16];
      test_assert_critical(L != NULL);
      test_assert(r->monodeg == j);
      for (int i = 0 ; i <= k ; i++) {
         double target = 0;
         if (i == j) target = pow(a,i-1) * b;
         test_assert(fabs(r->coeff[i]-target) < 1e-9);
      }
   }
   // Last element is F.
   F = L->term[33];
   for (int i = 0 ; i <= 16 ; i++) {
      double target = pow(a,i);
      test_assert(fabs(target-F->coeff[i]) < 1e-9);
   }
   for (int i = 17 ; i <= k ; i++) {
      test_assert(F->coeff[i] == 0);
   }
   
   // Next 16 series of rows.
   for (int j = 1 ; j <= 16 ; j++) {
      // First term is an R polynomial.
      R = L->term[j*dim0];
      test_assert_critical(R != NULL);
      if (j == 16) {
         test_assert(R->monodeg == 1);
         for (int i = 1 ; i <= k ; i++) {
            double target = 0;
            if (i == 1) target = d;
            test_assert(fabs(R->coeff[i]-target) < 1e-9);
         } 
      }
      else {
         test_assert(R->monodeg > k);
         test_assert(R->coeff[0] == 0);
         for (int i = 1 ; i <= 17-j ; i++) {
            double target = pow(a,i-1) * d;
            test_assert(fabs(R->coeff[i]-target) < 1e-9);
         } 
         for (int i = 17-j+1 ; i <= k ; i++) {
            test_assert(R->coeff[i] == 0);
         }
      }
      // Next 16 terms are r+ polynomials.
      for (int i = 1 ; i <= j ; i++) {
         r = L->term[j*dim0+i];
         test_assert(r == NULL);
      }
      for (int i = j+1 ; i <= 16 ; i++) {
         r = L->term[j*dim0+i];
         test_assert_critical(r != NULL);
         test_assert(r->monodeg == i-j);
         for (int l = 0 ; l <= k ; l++) {
            double target = 0;
            if (l == i-j) target = pow(a,i-j-1) * c;
            test_assert(fabs(r->coeff[l]-target) < 1e-9);
         }
      }
      // Next 16 terms are r- polynomials.
      for (int i = 1 ; i <= 17-j ; i++) {
         r = L->term[j*dim0+16+i];
         test_assert_critical(r != NULL);
         test_assert(r->monodeg == i);
         for (int l = 0 ; l <= k ; l++) {
            double target = 0;
            if (l == i) target = pow(a,i-1) * b;
            test_assert(fabs(r->coeff[l]-target) < 1e-9);
         }
      }
      for (int i = 17-j+1 ; i <= 16 ; i++) {
         r = L->term[j*dim0+16+i];
         test_assert(r == NULL);
      }
   }

   // Next 16 series of rows.
   for (int j = 1 ; j <= 16 ; j++) {
      // First term is an R polynomial.
      R = L->term[(j+16)*dim0];
      test_assert_critical(R != NULL);
      if (j == 16) {
         test_assert(R->monodeg == 1);
         for (int i = 1 ; i <= k ; i++) {
            double target = 0;
            if (i == 1) target = d;
            test_assert(fabs(R->coeff[i]-target) < 1e-9);
         } 
      }
      else {
         test_assert(R->monodeg > k);
         test_assert(R->coeff[0] == 0);
         for (int i = 1 ; i <= 17-j ; i++) {
            double target = pow(a,i-1) * d;
            test_assert(fabs(R->coeff[i]-target) < 1e-9);
         } 
         for (int i = 17-j+1 ; i <= k ; i++) {
            test_assert(R->coeff[i] == 0);
         }
      }
      // Next 16 terms are r+ polynomials.
      for (int i = 1 ; i <= 17-j ; i++) {
         r = L->term[(j+16)*dim0+i];
         test_assert_critical(r != NULL);
         test_assert(r->monodeg == i);
         for (int l = 0 ; l <= k ; l++) {
            double target = 0;
            if (l == i) target = pow(a,i-1) * c;
            test_assert(fabs(r->coeff[l]-target) < 1e-9);
         }
      }
      for (int i = 17-j+1 ; i <= 16 ; i++) {
         r = L->term[(j+16)*dim0+i];
         test_assert(r == NULL);
      }
      // Next 16 terms are r- polynomials.
      for (int i = 1 ; i <= j ; i++) {
         r = L->term[(j+16)*dim0+i+16];
         test_assert(r == NULL);
      }
      for (int i = j+1 ; i <= 16 ; i++) {
         r = L->term[(j+16)*dim0+i+16];
         test_assert_critical(r != NULL);
         test_assert(r->monodeg == i-j);
         for (int l = 0 ; l <= k ; l++) {
            double target = 0;
            if (l == i-j) target = pow(a,i-j-1) * b;
            test_assert(fabs(r->coeff[l]-target) < 1e-9);
         }
      }
   }
   
   destroy_mat(L);
   clean_mem_prob();

}


void
test_error_new_matrix_L
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   matrix_t *L;

   set_alloc_failure_countdown_to(0);
   redirect_stderr();
   // The error is that 'malloc()' will fail.
   L = new_matrix_L(.05);
   unredirect_stderr();
   reset_alloc();

   test_assert(L == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_n");

   set_alloc_failure_countdown_to(1);
   redirect_stderr();
   // The error is that 'malloc()' will fail (later).
   L = new_matrix_L(.05);
   unredirect_stderr();
   reset_alloc();

   test_assert(L == NULL);
   test_assert_stderr("[mem_seed_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_matrix_mult
(void)
{

   size_t k = 50;

   int success = set_params_mem_prob(17, k, 0.01);
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
      test_assert(tmp1->term[i]->monodeg > k);
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
      test_assert(tmp2->term[i]->monodeg > k);
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

   int success = set_params_mem_prob(17, 50, 0.01);
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
test_compute_exact_seed_prob
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   trunc_pol_t *w = compute_exact_seed_prob();
   test_assert_critical(w != NULL);

   double array[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      0.157056806616073, 0.148627374682234, 0.140197942748395,
      0.131768510814555, 0.123339078880716, 0.114909646946877,
      0.106480215013038, 0.0980507830791983, 0.0896213511453591,
      0.0811919192115198, 0.0727624872776805, 0.0643330553438413,
      0.055903623410002, 0.0474741914761627, 0.0390447595423234,
      0.0306153276084842, 0.0221858956746449, 0.0137564637408056,
      0.0124325640796893, 0.0111797197413002, 0.0099979307256383,
      0.00888719703270364, 0.00784751866249621, 0.00687889561501601,
      0.00598132789026304, 0.0051548154882373, 0.00439935840893878,
      0.0037149566523675, 0.00310161021852345, 0.00255931910740662,
      0.00208808331901703, 0.00168790285335466, 0.00135877771041952,
      0.00110070789021162, 0.000913693392730939, 0.00079773421797749,
      0.000692934765304654, 0.000598696078705965, 0.000514419202174956,
      0.000439505179705163, 0.000373355055290118, 0.000315369872923355,
      0.000264950676598408, 0.000221498510308811, 0.000184414418048099,
      0.000153099443809803, 0.00012695463158746, 0.000105381025374601,
      8.77796691647619e-05, 7.35516069514757e-05, 6.20978827282764e-05,
      5.28195404886977e-05, 4.51176242262736e-05, 3.83931779345379e-05,
      3.25521314958114e-05, 2.7505463651303e-05, 2.31692020011093e-05,
      1.94644230042148e-05, 1.6317251978492e-05, 1.3658863100701e-05,
      1.14254794064901e-05, 9.55837279039515e-06, 8.00386400584013e-06,
      6.71332266513674e-06, 5.64316723948461e-06, 4.75486505897122e-06,
      4.01493231257193e-06, 3.39493404814997e-06, 2.87148417245643e-06,
      2.42624545113028e-06, 2.04592950869837e-06, 1.72229682857541e-06,
      1.44790085183008e-06, 1.21604541817273e-06, 1.02074220694301e-06,
      8.56668178097524e-07, 7.19123013197522e-07, 6.03986556396534e-07,
      5.07676255428043e-07, 4.27104602593146e-07, 3.59636575748211e-07,
      3.0304707929254e-07, 2.55478385156033e-07, 2.15397573786844e-07};
   
   test_assert(w->monodeg > 50);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(fabs(w->coeff[i]-array[i]) < 1e-9);
   }

   free(w);

   clean_mem_prob();

}


void
test_compute_mem_prob_wgf
(void)
{

   trunc_pol_t *w0  = NULL;
   trunc_pol_t *w1  = NULL;
   trunc_pol_t *w2  = NULL;
   trunc_pol_t *w3  = NULL;
   trunc_pol_t *w4  = NULL;
   trunc_pol_t *w20 = NULL;

   int success = set_params_mem_prob(17, 19, 0.01);
   test_assert_critical(success);

   w0 = compute_mem_prob_wgf(.05,0);
   w1 = compute_mem_prob_wgf(.05,1);
   w2 = compute_mem_prob_wgf(.05,2);
   w3 = compute_mem_prob_wgf(.05,3);
   w4 = compute_mem_prob_wgf(.05,4);
   w20 = compute_mem_prob_wgf(.05,20);

   test_assert_critical(w0 != NULL);
   test_assert_critical(w1 != NULL);
   test_assert_critical(w2 != NULL);
   test_assert_critical(w3 != NULL);
   test_assert_critical(w4 != NULL);
   test_assert_critical(w20 != NULL);

   // The first terms can be computed directly.
   for (int i = 0 ; i < 17 ; i++) {
      test_assert(fabs(w0->coeff[i]-1) < 1e-9);
   }

   for (int i = 0 ; i < 17 ; i++) {
      test_assert(fabs(w1->coeff[i]-1) < 1e-9);
   }

   for (int i = 0 ; i < 17 ; i++) {
      test_assert(fabs(w2->coeff[i]-1) < 1e-9);
   }

   // N = 2, k = 17.
   double target_17 = 1-pow(.99,17);
   test_assert(fabs(w2->coeff[17]-target_17) < 1e-9);


   // N = 2, k = 18.
   double target_18;
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,2);
   test_assert(fabs(w2->coeff[18]-target_18) < 1e-9);

   // N = 3, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,3);
   test_assert(fabs(w3->coeff[18]-target_18) < 1e-9);

   // N = 4, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,4);
   test_assert(fabs(w4->coeff[18]-target_18) < 1e-9);

   // N = 20, k = 18.
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,20);
   test_assert(fabs(w20->coeff[18]-target_18) < 1e-9);

   // N = 1, k = 19.
   double target_19;
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * (1-pow(.95,18)*.05/3) - \
      2*.01*pow(.99,18) * (1-pow(.95,17)*.05/3) - \
      2*.01*.01*pow(.99,17) * (1-pow(.95,17)*.05/3) - \
      .01*.01*pow(.99,17) * (1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3);
   test_assert(fabs(w1->coeff[19]-target_19) < 1e-9);

   // N = 2, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,2) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,2) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,2) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,2);
   test_assert(fabs(w2->coeff[19]-target_19) < 1e-9);

   // N = 3, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,3) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,3) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,3) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,3);
   test_assert(fabs(w3->coeff[19]-target_19) < 1e-9);

   // N = 4, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,4) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,4) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,4) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,4);
   test_assert(fabs(w4->coeff[19]-target_19) < 1e-9);

   // N = 20, k = 19.
   target_19 = 1-pow(.99,19) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,18)*.05/3,20) - \
      2*.01*pow(.99,18) * pow(1-pow(.95,17)*.05/3,20) - \
      2*.01*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,20) - \
      .01*.01*pow(.99,17) * pow(1-pow(.95,17)*(.05/3)*(.05/3) -
            2*pow(.95,17)*(1-.05/3)*.05/3,20);
   test_assert(fabs(w20->coeff[19]-target_19) < 1e-9);

   free(w0);
   free(w1);
   free(w2);
   free(w3);
   free(w4);
   free(w20);

   // Other test case with longer seeds.
   
   success = set_params_mem_prob(20, 21, 0.02);
   test_assert_critical(success);

   // Compute with maximum precision.
   set_mem_prob_max_prcsn();

   w0 = compute_mem_prob_wgf(.05,0);
   w1 = compute_mem_prob_wgf(.05,1);
   w2 = compute_mem_prob_wgf(.05,2);
   w3 = compute_mem_prob_wgf(.05,3);
   w4 = compute_mem_prob_wgf(.05,4);

   test_assert_critical(w0 != NULL);
   test_assert_critical(w1 != NULL);
   test_assert_critical(w2 != NULL);
   test_assert_critical(w3 != NULL);
   test_assert_critical(w4 != NULL);

   for (int i = 0 ; i < 20 ; i++) {
      test_assert(fabs(w0->coeff[i]-1) < 1e-9);
      test_assert(fabs(w1->coeff[i]-1) < 1e-9);
      test_assert(fabs(w2->coeff[i]-1) < 1e-9);
   }

   const double target_20 = 1-pow(.98,20);
   test_assert(fabs(w0->coeff[20]-target_20) < 1e-9);
   test_assert(fabs(w1->coeff[20]-target_20) < 1e-9);
   test_assert(fabs(w2->coeff[20]-target_20) < 1e-9);
   test_assert(fabs(w3->coeff[20]-target_20) < 1e-9);

   double target_21;

   // Special case N = 0.
   target_21 = 1-pow(.98,21) - 2*.02*pow(.98,20);
   test_assert(fabs(w0->coeff[21]-target_21) < 1e-9);

   // Special case N = 1.
   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * (1-pow(.95,20)*.05/3);
   test_assert(fabs(w1->coeff[21]-target_21) < 1e-9);

   // Cases N > 1.
   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,2);
   test_assert(fabs(w2->coeff[21]-target_21) < 1e-9);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,3);
   test_assert(fabs(w3->coeff[21]-target_21) < 1e-9);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,4);
   test_assert(fabs(w4->coeff[21]-target_21) < 1e-9);

   free(w0);
   free(w1);
   free(w2);
   free(w3);
   free(w4);

   unset_mem_prob_max_prcsn();

   // Recompute without maximum precision. Test relative
   // error (results must be within 1% of the target).

   w0 = compute_mem_prob_wgf(.05,0);
   w1 = compute_mem_prob_wgf(.05,1);
   w2 = compute_mem_prob_wgf(.05,2);
   w3 = compute_mem_prob_wgf(.05,3);
   w4 = compute_mem_prob_wgf(.05,4);

   test_assert_critical(w0 != NULL);
   test_assert_critical(w1 != NULL);
   test_assert_critical(w2 != NULL);
   test_assert_critical(w3 != NULL);
   test_assert_critical(w4 != NULL);

   for (int i = 0 ; i < 20 ; i++) {
      test_assert(fabs(w0->coeff[i]-1) < 1e-2);
      test_assert(fabs(w1->coeff[i]-1) < 1e-2);
      test_assert(fabs(w2->coeff[i]-1) < 1e-2);
   }

   test_assert(fabs(w0->coeff[20]/target_20 - 1) < 1e-2);
   test_assert(fabs(w1->coeff[20]/target_20 - 1) < 1e-2);
   test_assert(fabs(w2->coeff[20]/target_20 - 1) < 1e-2);
   test_assert(fabs(w3->coeff[20]/target_20 - 1) < 1e-2);

   // Special case N = 0.
   target_21 = 1-pow(.98,21) - 2*.02*pow(.98,20);
   test_assert(fabs(w0->coeff[21]/target_21 - 1) < 1e-2);

   // Special case N = 1.
   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * (1-pow(.95,20)*.05/3);
   test_assert(fabs(w1->coeff[21]/target_21 - 1) < 1e-2);

   // Cases N > 1.
   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,2);
   test_assert(fabs(w2->coeff[21]/target_21 - 1) < 1e-2);

   // Just make sure that max precision is off (this one
   // requires more iterations to be accurate to within 1e-9).
   test_assert(fabs(w3->coeff[21]-target_21) > 1e-9);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,3);
   test_assert(fabs(w3->coeff[21]/target_21 - 1) < 1e-2);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,4);
   test_assert(fabs(w4->coeff[21]/target_21 - 1) < 1e-2);

   free(w0);
   free(w1);
   free(w2);
   free(w3);
   free(w4);

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
   
   size_t k = 35;
   const double u = 0.05; // Needed to compute 'xi' and 'eta'.

   int success = set_params_mem_prob(17, k, 0.01);
   test_assert_critical(success);

   matrix_t *M  = NULL;

   trunc_pol_t *tmp = tmp = new_zero_trunc_pol();
   test_assert_critical(tmp != NULL);

   for (int N = 1 ; N < 36 ; N++) {

      int dim = 17+N+1;

      // Transfer matrix.
      M = new_matrix_M(.05,N);
      test_assert_critical(M != NULL);
      test_assert(M->dim == dim);

      // State vectors.
      trunc_pol_t **x = malloc(dim * sizeof(trunc_pol_t *));
      trunc_pol_t **y = malloc(dim * sizeof(trunc_pol_t *));
      test_assert_critical(x != NULL);
      test_assert_critical(y != NULL);
      for (int i = 0 ; i < dim ; i++) {
         x[i] = new_zero_trunc_pol();
         y[i] = new_zero_trunc_pol();
         test_assert_critical(x[i] != NULL);
         test_assert_critical(y[i] != NULL);
      }

      // One segment from head to all down/m states (except 'm' = 0).
      for (int j = 1 ; j <= N ; j++) {
         // Store the end state in 'S[j]'.
         trunc_pol_update_add(x[j], M->term[j]);
      }
      // One segment from head to all up/i and
      // one segment to all down/m states.
      for (int i = 1 ; i <= 16 ; i++) {
      for (int j = 1 ; j <= N ; j++) {
         // Store the end state in 'S[j]'.
         trunc_pol_update_add(x[j], trunc_pol_mult(tmp,
            M->term[N+i], M->term[(N+i)*dim+j]));
      }
      }
      // One segment from all down/m states to tail.
      for (int j = 1 ; j <= N ; j++) {
         trunc_pol_update_add(y[j], M->term[j*dim+(dim-1)]);
      }
      // One segment from all down/m to all up/i
      // and one segment to tail.
      for (int i = 1 ; i <= 16 ; i++) {
      for (int j = 1 ; j <= N ; j++) {
         trunc_pol_update_add(y[j], trunc_pol_mult(tmp,
            M->term[j*dim+(N+i)], M->term[(N+i)*dim+(dim-1)]));
      }
      }

      // Combine.
      trunc_pol_t *w = new_zero_trunc_pol();
      for (int j = 1 ; j <= N ; j++) {
         trunc_pol_update_add(w, trunc_pol_mult(tmp, x[j], y[j]));
      }

      double target = 0.0;
      // The position of the error is [j] (1-based).
      for (int j = 1 ; j <= 17 ; j++) {
         double eta_term = 1.0 - pow(1.0-u,k-j) * u/3;
          target += 2 * (1.0 - pow(eta_term,N));
      }
      for (int j = 18 ; j <= k-17 ; j++) {
         double first_eta_term = 1.0 - pow(1.0-u,j-1) * u/3;
         double second_eta_term = 1.0 - pow(1.0-u,k-j) * u/3;
         target +=  1 - ( pow(first_eta_term, N) + 
               pow(second_eta_term, N) -
               pow(1-.05/3 + .05/3 * xi(j-1,u) * xi(k-j,u),N) );
      }
      // Multiply by the probability that there is one error.
      target *= .01 * pow(.99, k-1);
      test_assert(fabs(target - w->coeff[k]) < 1e-12);

      destroy_mat(M);
      for (int i = 0 ; i < dim ; i++) free(x[i]);
      for (int i = 0 ; i < dim ; i++) free(y[i]);
      free(w);
      free(x);
      free(y);

   }

   free(tmp);

   clean_mem_prob();

}


void
test_mcmc_method
(void)
{

   trunc_pol_t *w  = NULL;
   trunc_pol_t *mc = NULL;

   int success = set_params_mem_prob(17, 50, 0.01);
   test_assert_critical(success);

   w = compute_mem_prob_wgf(.05,5);
   test_assert_critical(w != NULL);

   mc = compute_mem_prob_mcmc(.05,5);
   test_assert_critical(mc != NULL);

   test_assert(mc->monodeg > 50);
   // First values are equal to 1 (the precision
   // is not very high for these first values).
   for (int i = 0 ; i < 17 ; i++) {
      test_assert(fabs(mc->coeff[i]-1) < 1e-7);
   }

   // Check that the MCMC esimates are within 2.2 standard deviations of
   // the exact value (note that the Gaussian approximation becomes bad
   // at the tail of the binomial distribution).
   for (int i = 17 ; i <= 50 ; i++) {
      double SD = sqrt(w->coeff[i] * (1-w->coeff[i]) / MCMC_RESAMPLINGS);
      test_assert(fabs(mc->coeff[i] - w->coeff[i]) < 2.2 * SD);
   }

   free(w);
   free(mc);

   clean_mem_prob();

}


// Test cases for export.
const test_case_t test_cases_mem_seed_prob[] = {
   // Math functions
   {"omega",                       test_omega},
   {"psi",                         test_psi},
   {"zeta",                        test_zeta},
   // Initialization functions
   {"set_params_mem_prob",         test_set_params_mem_prob},
   {"error_set_params_mem_prob",   test_error_set_params_mem_prob},
   {"uninitialized_error",         test_uninitialized_error},
   // Basic truncated polynomial constructors
   {"new_zero_trunc_pol",          test_new_zero_trunc_pol},
   {"error_new_zero_trunc_pol",    test_error_new_zero_trunc_pol},
   // Truncated polynomials manipulation
   {"trunc_pol_updated_add",       test_trunc_pol_update_add},
   {"trunc_pol_mult",              test_trunc_pol_mult},
   // Specific polynomials.
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
   {"new_trunc_pol_F",             test_new_trunc_pol_F},
   {"error_new_trunc_pol_F",       test_error_new_trunc_pol_F},
   {"new_trunc_pol_R",             test_new_trunc_pol_R},
   {"error_new_trunc_pol_R",       test_error_new_trunc_pol_R},
   {"new_trunc_pol_r",             test_new_trunc_pol_r},
   {"error_new_trunc_pol_r",       test_error_new_trunc_pol_r},
   // Matrix constructors.
   {"new_null_matrix",             test_new_null_matrix},
   {"error_new_null_matrix",       test_error_new_null_matrix},
   {"new_zero_matrix",             test_new_zero_matrix},
   {"error_new_zero_matrix",       test_error_new_zero_matrix},
   // Specific matrices
   {"new_matrix_M",                test_new_matrix_M},
   {"error_new_matrix_M",          test_error_new_matrix_M},
   {"new_matrix_L",                test_new_matrix_L},
   {"error_new_matrix_L",          test_error_new_matrix_L},
   // Matrix manipulation functions.
   {"matrix_mult",                 test_matrix_mult},
   {"error_matrix_mult",           test_error_matrix_mult},
   // Probabilities.
   {"compute_exact_seed_prob",     test_compute_exact_seed_prob},
   {"compute_mem_prob_wgf",        test_compute_mem_prob_wgf},
   {"misc_correctness",            test_misc_correctness},
   {"mcmc_method",                 test_mcmc_method},
#if 0
   {"average_errors",              test_average_errors},
   {"error_mem_seed_prob",         test_error_mem_seed_prob},
#endif
   {NULL, NULL},
};
