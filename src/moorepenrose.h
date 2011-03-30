#include <gsl/gsl_matrix.h>

#define CHOL_TOL            (double) 1.e-12

int moore_penrose(gsl_matrix* pA, double tol);

// not needed after gsl 1.12
//int gsl_linalg_cholesky_invert (gsl_matrix * cholesky);
void gslprint(gsl_matrix* pmymat);
void gslprint(gsl_vector* pmymat);
