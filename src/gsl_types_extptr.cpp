#include "Logger.h"
#include "fgls.h"
#include "gsl_types_extptr.h"

#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif

	void check_gslMatrix_Pointer(SEXP s) {
		if (TYPEOF(s) != EXTPTRSXP) {
			errorLog << "Pointer is not EXTPTRSXP" << endl << errorExit;
		}
		if (R_ExternalPtrTag(s) != install("gslMatrix"))
		{
			errorLog << "R_ExternalPtrTag(s) = " << (void*)R_ExternalPtrTag(s) << endl;
			errorLog << "Pointer is not gslMatrix" << endl << errorExit;
		}
	}

	void check_gslVector_Pointer(SEXP s) {
		if (TYPEOF(s) != EXTPTRSXP) {
			errorLog << "Pointer is not EXTPTRSXP" << endl << errorExit;
		}
		if (R_ExternalPtrTag(s) != install("gslVector"))
		{
			errorLog << "R_ExternalPtrTag(s) = " << (void*)R_ExternalPtrTag(s) << endl;
			errorLog << "Pointer is not gslVector" << endl << errorExit;
		}
	}

	static void gslMatrixRFinalizer(SEXP x) {
		check_gslMatrix_Pointer(x);
		if (x == R_NilValue) return;
		gsl_matrix * p = (gsl_matrix *) EXTPTR_PTR(x);
		if (p == NULL) return;
		wrapperLog << "Finalizing gslMatrix: "<< (void *)p << endl;
		gsl_matrix_free (p);
	}

	static void gslVectorRFinalizer(SEXP x) {
		check_gslVector_Pointer(x);
		if (x == R_NilValue) return;
		gsl_vector * p = (gsl_vector *) EXTPTR_PTR(x);
		if (p == NULL) return;
		wrapperLog << "Finalizing gslMatrix: "<< (void *)p << endl;
		gsl_vector_free (p);
	}

	gsl_matrix * gsl_matrix_from_EXTPTRSEXP(SEXP s){
		check_gslMatrix_Pointer(s);
		if (TYPEOF(s) == EXTPTRSXP) {
			return  (gsl_matrix*) R_ExternalPtrAddr(s);
		}
		errorLog << "External pointer not valid!" << endl << errorExit ;
		return NULL;
	}

	gsl_vector * gsl_vector_from_EXTPTRSEXP(SEXP s){
		check_gslVector_Pointer(s);
		if (TYPEOF(s) == EXTPTRSXP) {
			return  (gsl_vector*) R_ExternalPtrAddr(s);
		}
		errorLog << "External pointer not valid!" << endl << errorExit ;
		return NULL;
	}

	SEXP gslMatrix_from_gsl_matrix(gsl_matrix * X) {
		SEXP val = R_MakeExternalPtr(X, install("gslMatrix"), R_NilValue);
		return val;
	}

	SEXP gslMatrix_from_SEXP(SEXP RinX) {
		unsigned int nrow = (INTEGER(getAttrib(RinX,R_DimSymbol))[0]);
		unsigned int ncol = (INTEGER(getAttrib(RinX,R_DimSymbol))[1]);
		double * indata = REAL(RinX);
		gsl_matrix * X = gsl_matrix_alloc(nrow,ncol);

		if (X == NULL) {
			errorLog << "Error creating gslMatrix" << endl;
			return R_NilValue;
		}

		for (unsigned int i=0;i<nrow;i++)
			for (unsigned int j=0;j<ncol;j++)
				gsl_matrix_set(X,i,j,indata[nrow*j+i]);

		SEXP val = R_MakeExternalPtr(X, install("gslMatrix"), R_NilValue);
		R_RegisterCFinalizerEx(val, gslMatrixRFinalizer, (Rboolean) TRUE);
		return val;
	}

	SEXP gslVector_from_SEXP(SEXP RinY) {
		unsigned int length = LENGTH(RinY);
		double * indata = REAL(RinY);
		gsl_vector * Y = gsl_vector_alloc(length);

		if (Y == NULL) {
			errorLog << "Error creating gslVector" << endl;
			return R_NilValue;
		}

		for (unsigned int i=0;i<length;i++)
				gsl_vector_set(Y,i,indata[i]);

		SEXP val = R_MakeExternalPtr(Y, install("gslVector"), R_NilValue);
		R_RegisterCFinalizerEx(val, gslVectorRFinalizer, (Rboolean) TRUE);
		return val;
	}



#ifdef __cplusplus
}
#endif
