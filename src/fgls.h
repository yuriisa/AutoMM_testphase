#include <new>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <R.h>
#include <Rinternals.h>



extern "C" {

	double fgls(gsl_vector * Y, gsl_matrix * X,
			unsigned short int test, // 1 = wald, 2 = score, 3 = robust
			gsl_matrix * W,  // can be NULL
			gsl_matrix * tXWfixed, 	// can be NULL, otherwise contains pre-computed t(X)%*%W,
			// where only elements t(X[,!whichtest]) %*% W are fixed and
			// pre-computed
			unsigned int WhichToTest, // at what column of X 'fixed' part starts
			double scorevar,   // can be null, otherwise residual variance used in score test
			// ON RETURN
			gsl_vector * beta, // length ncolX, estimates
			gsl_matrix * V
	);

	SEXP fgls_caller (SEXP RinY, SEXP RinX,
			SEXP RintXWfixed, SEXP RinW,
			SEXP WTC, SEXP RinTest, SEXP RinScoreVar);


	// end extern "C"

}
