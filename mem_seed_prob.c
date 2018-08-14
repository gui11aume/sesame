#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "random.h"

// MACROS //

#define LIBNAME "mem_seed_prob"
#define VERSION "1.0 08-03-2018"

#define YES 1
#define NO  0

#define SUCCESS 1
#define FAILURE 0

// Maximum allowed number of duplicates.
#define MAXN 1024

// Prob that one of m altnerative threads survives i steps.
#define xi(i)  (      1.0 - pow(1.0-U,(i))            )
#define eta(i) (      1.0 - pow(1.0-U,(i)) * U/3      )

#define iszero(p) (p == NULL || (p->monodeg == 0 && p->coeff[0] == 0))

// Macro to simplify error handling.
#define handle_memory_error(x) do { \
   if ((x) == NULL) { \
      warning("memory_error", __func__, __LINE__); \
      ERRNO = __LINE__; \
      goto in_case_of_failure; \
}} while (0)


//  TYPE DECLARATIONS  //

typedef struct trunc_pol_t  trunc_pol_t;
typedef struct matrix_t     matrix_t;

// The variable 'monodeg' is a number from 0 to 'K' to signify that the
// polynomial is a monomial. If the value is greater than 'K', the
// polynomial is not a monomial. The null polynomial is thus a monomial
// of degree 0 with a coefficient equal to 0.


struct trunc_pol_t {
   size_t monodeg;       // Monomial (if applicable).
   double coeff[];       // Terms of the polynomial.
};

struct matrix_t {
   const size_t    dim;  // Column / row number.
   trunc_pol_t * term[]; // Terms of the matrix.
};


// GLOBAL VARIABLES / CONSTANTS //

static size_t    G = 0;       // Minimum size of MEM seeds.
static size_t    K = 0;       // Max degree (read size).

static double    P = 0.0;     // Probability of a read error.
static double    U = 0.0;     // Divergence rate between duplicates.

static size_t    KSZ = 0;     // Size of the 'trunc_pol_t' struct.

static trunc_pol_t * TEMP = NULL;  // For matrix multipliciation.

static double * XI;     // Precomputed values of 'xi'.
static double * XIc;    // Precomputed values of '1-xi'.
static double * ETA;    // Precomputed values of 'eta'.
static double * ETAc;   // Precomputed values of '1-eta'.

// Computation method.
static enum meth_t { METHOD_AUTO, METHOD_WGF, METHOD_MCMC } METH = 0;

static double * ARRAY[MAXN] = {0};  // Store results indexed by N.

static int ERRNO = 0;

// Error message.
const char internal_error[] =
   "internal error (please contact guillaume.filion@gmail.com)";



// FUNCTION DEFINITIONS //

// IO and error report functions
int  get_mem_prob_error_code (void) { return ERRNO; }   // VISIBLE
void reset_mem_prob_error (void) { ERRNO = 0; }         // VISIBLE
void set_mem_prob_method (int m) { METH = m; }          // VISIBLE


void
warning
(
   const char * msg,
   const char * function,
   const int    line
)
// Print warning message to stderr.
{
   fprintf(stderr, "[%s] error in function `%s' (line %d): %s\n",
         LIBNAME, function, line, msg);
}


// Mathematical functions.

double
omega
(
   size_t m,
   size_t N
)
{
   if (m > N) {
      return 0.0;
   }
   double log_N_choose_m = lgamma(N+1)-lgamma(m+1)-lgamma(N-m+1);
   return exp(log_N_choose_m + (N-m)*log(1-U/3) + m*log(U/3));
}

double
psi
(
   size_t i,
   size_t m,
   size_t n,
   size_t r,
   size_t N
)
{

   if (N < (m+r)) return 0.0;

   // Take out the terms of the sum that do not depend on 'q'.
   const double lC = lgamma(m+1) + lgamma(N-m-r+1);
   const double lx = log(xi(i));
   double val = 0.0;
   size_t topq = n-r < m ? n-r : m;
   for (int q = 0 ; q <= topq ; q++) {
      val += exp((n-r-q) * lx - lgamma(q+1) - lgamma(m-q+1)
                     - lgamma(n-r-q+1) - lgamma(N-m-n+q+1) + lC);
   }
   return val;

}


double
zeta
(
   size_t i,
   size_t m,
   size_t n,
   size_t N
)
{

   if (N < m) return 0.0;

   // Take out the terms of the sum that do not depend on 'r'.
   const double lC = lgamma(N-m+1)+lgamma(n+1)+lgamma(N-n+1)-lgamma(N+1);
   const double lu = i * log(1-U);
   double val = 0.0;
   size_t topr = N-m < n ? N-m : n;
   for (int r = 1 ; r <= topr ; r++) {
      val += psi(i,m,n,r,N) * \
             exp(r*lu - lgamma(r+1) - lgamma(N-m-r+1) + lC);
   }
   return val;
}



// Initialization and clean up.

void
clean_mem_prob // VISIBLE //
(void)
{

   // Set global variables to "fail" values.
   G = 0;
   K = 0;
   P = 0.0;
   U = 0.0;
   KSZ = 0;

   // Free everything.
   free(TEMP); TEMP = NULL;
   free(XI); XI = NULL;
   free(XIc); XIc = NULL;
   free(ETA); ETA = NULL;
   free(ETAc); ETAc = NULL;

   for (int i = 0 ; i < MAXN ; i++) free(ARRAY[i]);
   bzero(ARRAY, MAXN * sizeof(double *));

   ERRNO = 0;

   return;

}


int
set_params_mem_prob // VISIBLE //
(
   size_t g,
   size_t k,
   double p,
   double u
)
// Initialize the global variables from user-defined values.
{

   // Check input
   if (g == 0 || k == 0) {
      ERRNO = __LINE__;
      warning("parameters g and k must greater than 0",
            __func__, __LINE__); 
      goto in_case_of_failure;
   }

   if (p <= 0.0 || p >= 1.0) {
      ERRNO = __LINE__;
      warning("parameter p must be between 0 and 1",
            __func__, __LINE__); 
      goto in_case_of_failure;
   }

   if (u <= 0.0 || u >= 1.0) {
      ERRNO = __LINE__;
      warning("parameter u must be between 0 and 1",
            __func__, __LINE__); 
      goto in_case_of_failure;
   }

   G = g;  // MEM size.
   K = k;  // Read size.
   P = p;  // Sequencing error.
   U = u;  // Divergence rate.

   // Precompute intermediate quantities.
   
   free(XI); XI = malloc((K+1) * sizeof(double));
   free(XIc); XIc = malloc((K+1) * sizeof(double));
   free(ETA); ETA = malloc((K+1) * sizeof(double));
   free(ETAc); ETAc = malloc((K+1) * sizeof(double));
   handle_memory_error(XI);
   handle_memory_error(XIc);
   handle_memory_error(ETA);
   handle_memory_error(ETAc);

   for (int i = 0 ; i <= K ; i++) {
      XI[i]   = 1 - pow(1-U,i);
      XIc[i]  = 1 - XI[i];
      ETA[i]  = 1 - pow(1-U,i) * U/3.0;
      ETAc[i] = 1 - ETA[i];
   }

   // All 'trunc_pol_t' must be congruent.
   KSZ = sizeof(trunc_pol_t) + (K+1) * sizeof(double);

   // Clean previous values (if any).
   for (int i = 0 ; i < MAXN ; i++) free(ARRAY[i]);
   bzero(ARRAY, MAXN * sizeof(double *));

   // Allocate or reallocate 'TEMP'.
   free(TEMP);
   TEMP = calloc(1, KSZ);
   handle_memory_error(TEMP);

   return SUCCESS;

in_case_of_failure:
   clean_mem_prob();
   return FAILURE;

}


// Definitions of weighted generating functions.

trunc_pol_t *
new_zero_trunc_pol
(void)
// IMPORTANT: do not call before set_params_mem_prob().
{

   if (G == 0 || K == 0 || P == 0 || U == 0 || KSZ == 0) {
      warning("parameters unset: call `set_params_mem_prob'",
            __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = calloc(1, KSZ);
   handle_memory_error(new);

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_A
(
   const size_t m,    // Initial state (down).
   const size_t n,    // Final state (down).
   const size_t N     // Number of duplicates.
)
{

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // When N is 0, these polynomials are non-null
   // only if 'm' and 'n' are zero.
   if (N == 0) {
      if (m == 0 && n == 0) {
         new->monodeg = 1;
         new->coeff[1] = P;
      }
      return new;
   }

   // This is not a monomial.
   new->monodeg = K+1;

   double omega_p_pow_of_q = omega(n,N) * P;
   for (int i = 1 ; i <= G ; i++) {
      new->coeff[i] = (1 - pow(xi(i-1),N)) * omega_p_pow_of_q;
      omega_p_pow_of_q *= (1.0-P);
   }
   for (int i = G+1 ; i <= K ; i++) {
      new->coeff[i] = (1 - pow(xi(i-1),m) *
            (1- zeta(i-1,m,n,N))) * omega_p_pow_of_q;
      omega_p_pow_of_q *= (1.0-P);
   }
   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_B
(
   const size_t i,    // Final state (up).
   const size_t N     // Number of duplicates.
)
{

   if (i < 1) {
      warning(internal_error, __func__, __LINE__);
      ERRNO = __LINE__;
      goto in_case_of_failure;
   }   

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   if (N == 0) {
      // Zero if i is not 1.
      if (i != 1)
         return new;
         new->monodeg = 1;
         new->coeff[1] = 1-P;
      return new;
   }

   // This is a monomial.
   new->monodeg = i;
   new->coeff[i] = (pow(xi(i),N) - pow(xi(i-1),N)) * pow(1-P,i);
   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_C
(
   const size_t m,    // Initial state (down).
   const size_t N     // Number of duplicates.
)
// Note: the polynomials are defined in the case m > N, but
// they have no meaning in the present context. Here they are
// simply not forbidden, but also not used (and not tested).
{

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // Special case that 'N' is 0. Assuming that 'm' is 0 (see
   // remark above), this is a monomial equal to 1.
   if (N == 0) {
      new->monodeg = 0;
      new->coeff[0] = 1.0;
      return new;
   }

   // This is not a monomial (except in the case that
   // 'G' is 1 and 'm' is 1, but this is too rare to justify
   // the extra code (it won't trigger a mistake anyway).
   new->monodeg = K+1;

   double pow_of_q = 1.0;
   for (int i = 0 ; i <= G-1 ; i++) {
      new->coeff[i] = (1 - pow(xi(i),N)) * pow_of_q;
      pow_of_q *= (1.0-P);
   }
   for (int i = G ; i <= K ; i++) {
      new->coeff[i] = (1 - pow(xi(i),m)) * pow_of_q;
      pow_of_q *= (1.0-P);
   }
   return new;

in_case_of_failure:
   return NULL;

}



trunc_pol_t *
new_trunc_pol_D
(
   const size_t j,    // Initial state (up).
   const size_t m,    // Final state (down).
   const size_t N     // Number of duplicates.
)
{

   if (j > G-1) {
      warning(internal_error, __func__, __LINE__);
      ERRNO = __LINE__;
      goto in_case_of_failure;
   }   

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // In the case that 'N' is 0, the only polynomial that is
   // defined is D(1,0,0) -- the others are set to 0. In the
   // case that 'm' is greater than 'N' the polynomial is 0.
   if ((N == 0 && !(j == 1 && m == 0)) || m > N) {
      return new;
   }

   // This is a monomial only if j is equal to G-1.
   new->monodeg = j == G-1 ? 1 : K+1;

   double omega_p_pow_of_q = omega(m,N) * P;
   for (int i = 1 ; i <= G-j ; i++) {
      new->coeff[i] = omega_p_pow_of_q;
      omega_p_pow_of_q *= (1.0-P);
   }

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_E
(
   const size_t j    // Initial state (up).
)
{

   if (j > G-1) {
      warning(internal_error, __func__, __LINE__);
      ERRNO = __LINE__;
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // This is a monomial only if j is equal to G-1
   new->monodeg = j == G-1 ? 0 : K+1;

   double pow_of_q = 1.0;
   for (int i = 0 ; i <= G-j-1 ; i++) {
      new->coeff[i] = pow_of_q;
      pow_of_q *= (1.0-P);
   }

   return new;

in_case_of_failure:
   return NULL;

}



matrix_t *
new_null_matrix
(
   const size_t dim
)
// Create a matrix where all truncated polynomials
// (trunc_pol_t) are set NULL.
{

   // Initialize to zero.
   size_t sz = sizeof(matrix_t) + dim*dim * sizeof(trunc_pol_t *);
   matrix_t *new = calloc(1, sz);
   handle_memory_error(new);

   // The dimension is set upon creation
   // and must never change afterwards.
   *(size_t *)&new->dim = dim;

   return new;

in_case_of_failure:
   return NULL;

}


void
destroy_mat
(
   matrix_t * mat
)
{

   // Do not try anything on NULL.
   if (mat == NULL) return;

   size_t nterms = (mat->dim)*(mat->dim);
   for (size_t i = 0 ; i < nterms ; i++)
      free(mat->term[i]);
   free(mat);

}



matrix_t *
new_zero_matrix
(
   const size_t dim
)
// Create a matrix where all terms
// are 'trunc_pol_t' struct set to zero.
{

   matrix_t *new = new_null_matrix(dim);
   handle_memory_error(new);

   for (int i = 0 ; i < dim*dim ; i++) {
      new->term[i] = new_zero_trunc_pol();
      handle_memory_error(new->term[i]);
   }

   return new;

in_case_of_failure:
   destroy_mat(new);
   return NULL;

}



matrix_t *
new_matrix_M
(
   const size_t N    // Number of duplicates.
)
{

   const size_t dim = G+N+1;
   matrix_t *M = new_null_matrix(dim);
   handle_memory_error(M);

   // First N+1 series of rows.
   for (int j = 0 ; j <= N ; j++) {
      for (int i = 0 ; i <= N ; i++) {
         M->term[j*dim+i] = new_trunc_pol_A(j,i,N);
         handle_memory_error(M->term[j*dim+i]);
      }
      for (int i = N+1 ; i <= N+G-1 ; i++) {
         M->term[j*dim+i] = new_trunc_pol_B(i-N,N);
         handle_memory_error(M->term[j*dim+i]);
      }
      M->term[j*dim+(N+G)] = new_trunc_pol_C(j,N);
   }

   // Next G-1 rows.
   for (int j = N+1 ; j <= N+G-1 ; j++) {
      for (int i = 0 ; i <= N ; i++) {
         M->term[j*dim+i] = new_trunc_pol_D(j-N,i,N);
         handle_memory_error(M->term[j*dim+i]);
      }
      M->term[j*dim+(N+G)] = new_trunc_pol_E(j-N);
   }

   // Last row is null (nothing to do).

   return M;

in_case_of_failure:
   destroy_mat(M);
   return NULL;

}


trunc_pol_t *
trunc_pol_mult
(
         trunc_pol_t * dest,
   const trunc_pol_t * a,
   const trunc_pol_t * b
)
{

   // Erase destination.
   bzero(dest, KSZ);

   // If any of the two k-polynomials is zero, keep 'dest' as
   // zero and return 'NULL' (not dest). See macro 'iszero'.
   if (iszero(a) || iszero(b)) return NULL;

   if (a->monodeg < K+1 && b->monodeg < K+1) {
      // Both are monomials, just do one multiplication.
      // If degree is too high, all coefficients are zero.
      if (a->monodeg + b->monodeg > K) return NULL;
      // Otherwise do the multiplication.
      dest->monodeg = a->monodeg + b->monodeg;
      dest->coeff[dest->monodeg] =
         a->coeff[a->monodeg] * b->coeff[b->monodeg];
      return dest;
   }

   // The result cannot be a monomial.
   dest->monodeg = K+1;
   if (a->monodeg < K+1) {
      // 'a' is a monomial, do one "row" of multiplications.
      for (int i = a->monodeg ; i <= K ; i++)
         dest->coeff[i] = a->coeff[a->monodeg] * b->coeff[i-a->monodeg];
   }
   else if (b->monodeg < K+1) {
      // 'b' is a monomial, do one "row" of multiplications.
      for (int i = b->monodeg ; i <= K ; i++)
         dest->coeff[i] = b->coeff[b->monodeg] * a->coeff[i-b->monodeg];
   }
   else {
      // Standard convolution product.
      for (int i = 0 ; i <= K ; i++) {
         dest->coeff[i] = a->coeff[0] * b->coeff[i];
         for (int j = 1 ; j <= i ; j++) {
            dest->coeff[i] += a->coeff[j] * b->coeff[i-j];
         }
      }
   }

   return dest;

}


void
trunc_pol_update_add
(
         trunc_pol_t * dest,
   const trunc_pol_t * a
)
{

   // No update when adding zero.
   if (iszero(a)) return;

   // Only adding two monomials of same degree returns
   // a monomial (the case of adding a zero polynomial
   // was taken care of above.
   if (a->monodeg != dest->monodeg) {
      dest->monodeg = K+1;
   }

   for (int i = 0 ; i <= K ; i++) {
      dest->coeff[i] += a->coeff[i];
   }

}


matrix_t *
matrix_mult
(
         matrix_t * dest,
   const matrix_t * a,
   const matrix_t * b
)
{

   if (a->dim != dest->dim || b->dim != dest->dim) {
      warning(internal_error, __func__, __LINE__);
      ERRNO = __LINE__;
      goto in_case_of_failure;
   }

   size_t dim = dest->dim;

   for (int i = 0 ; i < dim ; i++) {
   for (int j = 0 ; j < dim ; j++) {
      // Erase destination entry.
      bzero(dest->term[i*dim+j], KSZ);
      for (int m = 0 ; m < dim ; m++) {
         trunc_pol_update_add(dest->term[i*dim+j],
            trunc_pol_mult(TEMP, a->term[i*dim+m], b->term[m*dim+j]));
      }
   }
   }

   return dest;

in_case_of_failure:
   return NULL;

}


int
fault_in_params
(
   const size_t k,    // Segment or read size.
   const size_t N     // Number of duplicates.
)
{

   // Check if sequencing parameters were set. Otherwise
   // warn user and fail gracefully (return nan).
   if (G == 0 || K == 0 || P == 0 || U == 0 || KSZ == 0) {
      warning("parameters unset: call `set_params_mem_prob'",
            __func__, __LINE__);
      return YES;
   }

   // Check input.
   if (N > MAXN-1) {
      char msg[128];
      snprintf(msg, 128, "argument N greater than %d", MAXN);
      warning(msg, __func__, __LINE__);
      ERRNO = __LINE__;
      return YES;
   }

   if (k > K) {
      char msg[128];
      snprintf(msg, 128, "argument k greater than set value (%ld)", K);
      warning(msg, __func__, __LINE__);
      ERRNO = __LINE__;
      return YES;
   }

   return NO;

}


// Need this snippet to compute a bound of the numerical imprecision.
double HH(double x, double y)
         {  return x*log(x/y)+(1-x)*log((1-x)/(1-y));  }

trunc_pol_t *
compute_mem_prob_wgf
(
   const size_t N    // Number of duplicates.
)
{

   // Assume parameters were checked by the caller.

   // The truncated polynomial 'w' stands the weighted generating
   // function of the reads without on-target MEM seed for set
   // parameters and 'N' duplicates.
   trunc_pol_t *w = new_zero_trunc_pol();

   // The matrix 'M' is the transfer matrix of reads without MEM
   // seed. The row and column of the head state have been deleted
   // because the reads can be assumed to start in state "double-
   // down". The entries of 'M' are truncated polynomials that
   // stand for weighted generating functions of segments joining
   // different states. The m-th power of 'M' contains the weighted
   // generating functions of sequences of m segments joining the
   // dfferent states. We thus update 'w' with the entry of M^m
   // that joins the initial "double-down" state to the tail state.
   matrix_t *M = new_matrix_M(N);

   // Create two temporary matrices 'powM1' and 'powM2' to compute
   // the powers of M. On first iteration, powM1 = M*M, and later
   // perform the operations powM2 = M*powM1 and powM1 = M*powM2,
   // so that powM1 is M^2m at iteration m. Using two matrices
   // allows the operations to be performed without requesting more
   // memory. The temporary matrices are implicitly erased in the
   // course of the multiplication by 'matrix_mult()'.
   matrix_t *powM1 = new_zero_matrix(G+N+1);
   matrix_t *powM2 = new_zero_matrix(G+N+1);

   handle_memory_error(w);
   handle_memory_error(M);
   handle_memory_error(powM1);
   handle_memory_error(powM2);

   // Update weighted generating function with
   // one-segment reads (i.e. tail only).
   trunc_pol_update_add(w, M->term[G+N]);

   matrix_mult(powM1, M, M);

   // Update weighted generating function with two-segment reads.
   trunc_pol_update_add(w, powM1->term[G+N]);

   // There is at least one sequencing error for every three
   // segments. We bound the probability that a read of size k has
   // at least m/3 errors by a formula for the binomial distribution,
   // where m is the number of segments, i.e. the power of matirx M.
   // https://en.wikipedia.org/wiki/Binomial_distribution#Tail_Bounds
   for (int m = 2 ; m < K ; m += 2) {
      // Increase the number of segments and update
      // the weighted generating function accordingly.
      matrix_mult(powM2, M, powM1);
      trunc_pol_update_add(w, powM2->term[G+N]);
      matrix_mult(powM1, M, powM2);
      trunc_pol_update_add(w, powM1->term[G+N]);
#ifdef MAX_PRECISION_MODE
      // In max precision debug mode, get all possible digits.
      continue;
#endif
      // Otherwise, stop when reaching 1% precision.
      double x = floor((m+2)/3) / ((double) K);
      double bound_on_imprecision = exp(-HH(x, P)*K);
      if (bound_on_imprecision / w->coeff[K] < 1e-2) break;
   }

   // Clean temporary variables.
   destroy_mat(powM1);
   destroy_mat(powM2);
   destroy_mat(M);

   return w;

in_case_of_failure:
   // Clean everything.
   destroy_mat(powM1);
   destroy_mat(powM2);
   destroy_mat(M);
   free(w);
   return NULL;

}


size_t
rgeom
(void)
{
   return ceil( log(runifMT()) / log(1-P) );
}


size_t
rpos
(
   const size_t m,
   const size_t i
)
{
   if (m == 0) return 1;
   const double mth_root_of_unif = pow(runifMT(), 1.0/m);
   return ceil( log(1-XI[i]*mth_root_of_unif) / log(1-U) );
}


void
one_mcmc
(
   const size_t   N,
         double * pos
)
{

   size_t i;     // Size of the error-free segment.
   size_t m = 0; // Number of hard masks.

   // Note: stop if last nucleotide is an error.
   for (int sz = K ; sz > 0 ; sz -= (i+1)) {

      // Get size of the error-free segment.
      i = rgeom() - 1;

      if (i >= sz) {
         // If the error-free segment is longer than the read, resize
         // the error-free segment to 'sz' and ignore soft masks.
         
         // If the size of the read is shorter than the minimum
         // seed size 'g', then there is no on-target MEM seed.
         int long_enough = sz >= G;
         // Otherwise, there is an on-target MEM seed if there is
         // no survivng hard mask.
         int unmasked = rbinom(m, XIc[sz]) == 0;
         if (unmasked && long_enough) {
            // The segment is unmasked and long enough: there is a global
            // MEM seed. We need to find its position.
            size_t failpos = rpos(m, sz);
            size_t from = failpos < G ? K-sz + G : K-sz + failpos;
            for (int j = from ; j <= K ; j++) pos[j]++;
         }
         return;
      }
      else {
         // The next error is before the end of the read.
         // Get the number of hard masks that survive i
         // nucleotides -- each survives with probability (1-u)^i.
         size_t hsurv = rbinom(m, XIc[i]);

         // Get the number of soft masks that survive i+1
         // nucleotides -- each survives with probability (1-u)^i*u/3.
         size_t ssurv = rbinom(N-m, ETAc[i]);

         // If the segment has size greater than the
         // minimum seed length and there are no surviving
         // threads then there is an on-target MEM seed.
         int unmasked = hsurv == 0;
         int long_enough = i >= G;
         if (unmasked && long_enough) {
            int all_soft_masks_failed = ssurv == 0;
            size_t failpos = rpos(m, i);
            size_t from = failpos < G ? K-sz + G : K-sz + failpos;
            if (all_soft_masks_failed) {
               for (int j = from ; j <= K ; j++) pos[j]++;
               return;
            }
            else {
               size_t to = K-sz+i;
               for (int j = from ; j <= to ; j++) pos[j]++;
            }
         }

         // Otherwise, get the number of strictly masking threads (m).
         // All the surviving soft masks become strictly masking. The
         // hard masks have a probability u/3 of matching the error.
         // The soft masks that did not survive have a probability
         // pp = u/3 * xi / eta of matching the error. 
         double pp = 1 - (1 - U/3.0) / ETA[i];
         // The new state is down/m.
         m = ssurv + rbinom(m, U/3.0) + rbinom(N-m-ssurv, pp);

      }

   }

}


double *
compute_mem_prob_mcmc
(
   const size_t N    // Number of duplicates.
)
{

   // Assume parameters were checked by the caller.
   
   // Perform 10 million resamplings.
   // TODO: allow user to change this.
   const size_t R = 10000000;
   
   double * seeds = malloc((K+1) * sizeof(double));
   handle_memory_error(seeds);

   for (int i = 0 ; i < R ; i++) {
      one_mcmc(N, seeds);
   }

   for (int i = 0 ; i <= K ; i++) {
      seeds[i] = 1.0 - seeds[i] / R;
   }

   return seeds;

in_case_of_failure:
   return NULL;

}
   

double
mem_seed_prob
(
   const size_t k,       // Segment or read size.
   const size_t N        // Number of duplicates.
)
{
   
   // Check parameters.
   if (fault_in_params(k,N)) {
      goto in_case_of_failure;
   }

   // If results were not computed for given number of duplicates (N),
   // need to compute truncated genearting function and store the
   // results for future use.
   if (ARRAY[N] == NULL) {

      // Choose method. If N > 20 use MCMC.
      int use_method_wgf = METH == METHOD_WGF ||
                              (METH == METHOD_AUTO && N < 21);
      
      if (use_method_wgf) {
         trunc_pol_t *w = compute_mem_prob_wgf(N);
         if (w == NULL)
            goto in_case_of_failure;
         // Copy results to ARRAY.
         ARRAY[N] = malloc((K+1) * sizeof(double));
         handle_memory_error(ARRAY[N]);
         memcpy(ARRAY[N], w->coeff, (K+1) * sizeof(double));
         free(w);
      }
      else {
         ARRAY[N] = compute_mem_prob_mcmc(N);
         if (ARRAY[N] == NULL)
            goto in_case_of_failure;
      }

   }

   return ARRAY[N][k];

in_case_of_failure:
   return 0.0/0.0;  //  nan

}
