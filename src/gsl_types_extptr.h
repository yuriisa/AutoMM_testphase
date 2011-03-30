#ifdef __cplusplus
extern "C" {
#endif


	gsl_matrix * gsl_matrix_from_EXTPTRSEXP(SEXP s);
	gsl_vector * gsl_vector_from_EXTPTRSEXP(SEXP s);
	SEXP gslMatrix_from_gsl_matrix(gsl_matrix * X);
	SEXP gslMatrix_from_SEXP(SEXP RinX);


#ifdef __cplusplus
}
#endif
