
#include "Rstaff.h"
#include "iterator_functions.h"
#include "iterator.h"
#include "fgls.h"
#include "gsl_types_extptr.h"


#ifdef __cplusplus
extern "C" {
#endif



	void fglsWrapper(double *indata, unsigned long int indataHeight,
			unsigned long int indataWidth, double *outdata,
			unsigned long int &outdataNcol, unsigned long int &outdataNrow,
			unsigned int narg,
			SEXP *argList
			//SEXP RinY, [0]
			//SEXP RinX, [1]
			//SEXP RintXWfixed, [2]
			//SEXP RinW, [3]
			//SEXP RinTest, [4]
			//SEXP WhichAreGT, [5]
			//// c("additive","dominantB","recessiveB","overdominant","genotypic")
			//SEXP modelGT, [6]
			//// c("typed","dose","probability")
			//SEXP codingGT, [7]
			//SEXP mincall, [8]
			//SEXP RinScoreVar [9]
			//SEXP minimal [10]
	)
	{

		gsl_matrix * tmp_gsl_X = gsl_matrix_from_EXTPTRSEXP(argList[1]);
		unsigned int ncol = tmp_gsl_X->size2;
		unsigned int outsize_Ncol = 1+2*ncol+ncol*(ncol+1)/2;

		//unsigned int ncol = (INTEGER(getAttrib(RinX,R_DimSymbol))[1]);

		if (indata) {
			//SEXP X;
			unsigned int nrow = tmp_gsl_X->size1;
			gsl_matrix * gsl_X = gsl_matrix_alloc(nrow,ncol);
			gsl_matrix_memcpy(gsl_X,tmp_gsl_X);
			//gsl_matrix_fprintf(stdout,gsl_X,"gsl_X=%f");
			//unsigned int nrow = INTEGER(getAttrib(RinX,R_DimSymbol))[0];
			//PROTECT(X = allocMatrix(REALSXP, nrow, ncol));
			//double * pnew = REAL(X);
			//double * pold = REAL(RinX);
			int model_gt = INTEGER(argList[6])[0];
			int coding_gt = INTEGER(argList[7])[0];

			// impute iterated data
			for (unsigned int i = 0; i < indataWidth; i++) {
				double mean = 0.0;
				double fullN = 0.0;
				bool anynan = false;
				for (unsigned int j = 0; j < nrow; j++) {
					double tmp = indata[i*nrow + j];
					if (!isnan(tmp)) {mean+=tmp;fullN+=1.0;} else anynan=true;
				}
				if (anynan) {
					double SNPcall = fullN/((double) nrow);
					if (SNPcall < REAL(argList[8])[0]) {
						double zero = 0;
						Rprintf("Warning: call < mincall\n");
						for (unsigned int j = 0; j < outsize_Ncol; j++) outdata[j++] = 0/zero;
						return;
					}
					// HERE MUST TREAT 'impute' and 'exact'
					mean /= fullN;
					for (unsigned int j = 0; j < nrow; j++)
						if (isnan(indata[i*nrow + j])) // take care of NaNs later on...
							if (coding_gt != 1 || model_gt == 1)
								indata[i*nrow + j] = mean;
							else // take care of NANs later on
								indata[i*nrow + j] = -1.0; //mean;
				}
			}
			// arrange gt
			/**
			if (gtcoding == "dose" && model.SNP != "additive")
				stop("'dose' type only possible with 'additive' model")
			 **/
			// c("additive","dominantB","recessiveB","overdominant","genotypic")
			// c("typed","dose","probability")

			int nGTcols = 1; if (model_gt == 5) nGTcols = 2;
			double * processedGT = new (nothrow) double [nrow*nGTcols];
			if (coding_gt == 2 && model_gt != 1) {
				Rprintf("'dose' type only possible with 'additive' model\n");
				delete [] processedGT;
				return;
			}
			/**
			if (gtcoding == "probability") {
				if (model.SNP=="additive") return(matrix(as(genodata[,csnp],"double")*2+as(genodata[,csnp],"double")*1,ncol=1))
				GGG <- matrix(c(as(genodata[,csnp+1],"double"),as(genodata[,csnp],"double")),ncol=2)
			 **/
			if (coding_gt == 3) {
				if (model_gt == 1) { // additive
					for (unsigned int i = 0; i < nrow; i++) {
						processedGT[i] = indata[i]*2.0+indata[nrow+i];
					}
				} else if (model_gt == 2) { // dominantB
					for (unsigned int i = 0; i < nrow; i++) {
						processedGT[i] = indata[i]+indata[nrow+i];
					}
				} else if (model_gt == 3) { // recessiveB
					for (unsigned int i = 0; i < nrow; i++) {
						processedGT[i] = indata[i];
					}
				} else if (model_gt == 4) { // overdominant
					for (unsigned int i = 0; i < nrow; i++) {
						processedGT[i] = indata[nrow+i];
					}
				} else if (model_gt == 5) { // genotypic
					for (unsigned int i = 0; i < nrow; i++) {
						processedGT[i] = indata[nrow+i];
						processedGT[nrow+i] = indata[i];
					}
				} else {
					Rprintf("unrecognised model_gt\n");
					delete [] processedGT;
					return;
				}
			} else if (coding_gt == 1) {
				/**
			} else if (gtcoding == "typed") {
				if (model.SNP=="additive") return(matrix(as(genodata[,csnp],"double"),ncol=1))
				tmp <- as(genodata[,csnp],"double")
				tmp <- as.integer(tmp)
				GGG <- matrix(c(1*(tmp==1),1*(tmp==2)),ncol=2)
				 **/
				if (model_gt == 1) { // additive
					for (unsigned int i = 0; i < nrow; i++) {
						processedGT[i] = indata[i];
					}
				} else if (model_gt == 2) { // dominantB
					double nrel = 0.0, nnn = 0.0, mn = 0.0;
					for (unsigned int i = 0; i < nrow; i++) {
						if (abs(indata[i]-1.0) < 1.e-8 || abs(indata[i]-2.0) < 1.e-8) {
							processedGT[i] = 1.0;
							nrel+=1.0;
							nnn+=1.0;
						} else if (abs(indata[i])<1.e-8) {
							processedGT[i] = 0.0;
							nnn+=1.0;
						} else processedGT[i] = -1.0;
					}
					mn = nrel/nnn;
					for (unsigned int i = 0; i < nrow; i++) {
						if (processedGT[i] < (-0.1)) processedGT[i] = mn;
					}
				} else if (model_gt == 3) { // recessiveB
					double nrel = 0.0, nnn = 0.0, mn = 0.0;
					for (unsigned int i = 0; i < nrow; i++) {
						if (indata[i] > (2.0-1.e-8)) {
							processedGT[i] = 1.0;
							nrel++;
							nnn++;
						} else if (abs(indata[i]-1.0) < 1.e-8 || abs(indata[i]) < 1.e-8) {
							processedGT[i] = 0.0;
							nnn++;
						} else processedGT[i] = -1.0;
					}
					mn = nrel/nnn;
					for (unsigned int i = 0; i < nrow; i++)
						if (processedGT[i] < (-0.1)) processedGT[i] = mn;
				} else if (model_gt == 4) { // overdominant
					double nrel = 0.0, nnn = 0.0, mn = 0.0;
					for (unsigned int i = 0; i < nrow; i++) {
						if (abs(indata[i] - 1.0) < 1.e-8) {
							processedGT[i] = 1.0;
							nrel++;
							nnn++;
						} else if (abs(indata[i]-2.0) < 1.e-8 || abs(indata[i]) < 1.e-8) {
							processedGT[i] = 0.0;
							nnn++;
						} else processedGT[i] = -1.0;
					}
					mn = nrel/nnn;
					for (unsigned int i = 0; i < nrow; i++)
						if (processedGT[i] < (-0.1)) processedGT[i] = mn;
				} else if (model_gt == 5) { // genotypic
					double nrel1 = 0.0, nrel2 = 0.0, nnn = 0.0, mn1 = 0.0, mn2 = 0.0;
					for (unsigned int i = 0; i < nrow; i++) {
						if (abs(indata[i] - 1.0) < 1.e-8) {
							processedGT[i] = 1.0;
							processedGT[nrow+i] = 0.0;
							nrel1++;
							nnn++;
						} else if (abs(indata[i] - 2.0) < 1.e-8) {
							processedGT[i] = 0.0;
							processedGT[nrow+i] = 1.0;
							nrel2++;
							nnn++;
						} else if (abs(indata[i]) < 1.e-8) {
							processedGT[i] = 0.0;
							processedGT[nrow+i] = 0.0;
							nnn++;
						} else processedGT[i] = -1.0;
					}
					mn1 = nrel1/nnn;
					mn2 = nrel2/nnn;
					for (unsigned int i = 0; i < nrow; i++)
						if (processedGT[i] < (-0.1))
						{
							processedGT[i] = mn1;
							processedGT[nrow+i] = mn2;
						}
				} else {
					Rprintf("unrecognised model_gt\n");
					delete [] processedGT;
					return;
				}
			} else if (coding_gt == 2) {
				/**
			} else if (gtcoding == "dose") {
				if (model.SNP=="additive") return(matrix(as(genodata[,csnp],"double"),ncol=1))
				else stop("uhaha")
			} else {
				stop("brrrr")
			}
			if (model.SNP == "genotypic") return(GGG)
			else if (model.SNP == "dominantB") return(GGG[,1,drop=FALSE]+GGG[,2,drop=FALSE])
			else if (model.SNP == "recessiveB") return(GGG[,2,drop=FALSE])
			else if (model.SNP == "overdominant") return(GGG[,1,drop=FALSE])
				 **/
				if (model_gt == 1) { // additive
					for (unsigned int i = 0; i < nrow; i++) {
						processedGT[i] = indata[i];
					}
				} else {
					Rprintf("'dose' type only possible with 'additive' model (o2)\n");
					delete [] processedGT;
					return;
				}
			} else {
				Rprintf("unrecognised coding_gt\n");
				delete [] processedGT;
				return;
			}

			// copy constant data
			/**
			for (unsigned int i = 0; i < ncol; i++)
				for (unsigned int j = 0; j < nrow; j++) {
					if (INTEGER(WhichAreGT)[i] == 0)
						pnew[j + nrow*i] = pold[nrow*i + j];
					else if (INTEGER(WhichAreGT)[i] == 1)
						pnew[j + nrow*i] = pold[nrow*i + j]*processedGT[j];
					else if (INTEGER(WhichAreGT)[i] == 2)
						pnew[j + nrow*i] = pold[nrow*i + j]*processedGT[nrow+j];
					else {
						Rprintf("WhichAreGT[%d] == %d; not recognised\n",
								i,INTEGER(WhichAreGT)[i]);
						delete [] processedGT;
						return;
					}
				}
			 **/

			// copy iterated data
			for (unsigned int i = 0; i < ncol; i++)
				for (unsigned int j = 0; j < nrow; j++) {
					if (INTEGER(argList[5])[i] == 0) {}
					else if (INTEGER(argList[5])[i] == 1) {
						double inVal = gsl_matrix_get(gsl_X,j,i)*processedGT[j];
						gsl_matrix_set(gsl_X,j,i,inVal);
					} else if (INTEGER(argList[5])[i] == 2) {
						double inVal = gsl_matrix_get(gsl_X,j,i)*processedGT[nrow+j];
						gsl_matrix_set(gsl_X,j,i,inVal);
					} else {
						Rprintf("WhichAreGT[%d] == %d; not recognised\n",
								i,INTEGER(argList[5])[i]);
						delete [] processedGT;
						return;
					}
				}
			delete [] processedGT;

			// copy iterated data
			//for (unsigned int i = noldcol; i < ncol; i++)
			//	for (unsigned int j = 0; j < nrow; j++)
			//		pnew[i*nrow + j] = indata[(i-noldcol)*nrow + j];
			//		pnew[i*nrow + j] = indata[(i-noldcol) + indataWidth*j];

			// run function
			unsigned int WhichToTest = 0, kkk = 0;
			while (WhichToTest==0) {
				if (INTEGER(argList[5])[kkk] != 0) WhichToTest = kkk+1;
				kkk++;
			}
			//SEXP WTC;
			//PROTECT(WTC = allocVector(INTSXP, 1));
			//INTEGER(WTC)[0] = intWTC;
			//SEXP tmpout;
			//PROTECT(tmpout = allocVector(REALSXP, outsize_Ncol ));
			//SEXP in_gsl_X;
			//PROTECT(in_gsl_X = allocSExp(EXTPTRSXP) );
			//in_gsl_X = gslMatrix_from_gsl_matrix(gsl_X);

			//tmpout = fgls_caller ( argList[0],
			//		in_gsl_X, argList[2], argList[3], WTC, argList[4], argList[9]);

			// BEGIN part which was in fgsl_caller

			//unsigned int WhichToTest = INTEGER(WTC)[0];

			//gsl_matrix * X = gsl_matrix_from_EXTPTRSEXP(SEXP_ptr_to_gsl_X);
			//unsigned int nrowX = X->size1;
			//unsigned int ncolX = X->size2;

			//double * inY = REAL(RinY);
			//gsl_vector * Y = gsl_vector_alloc(nrowX);
			//for (unsigned int i=0;i<nrowX;i++)
			//	gsl_vector_set(Y,i,inY[i]);
			gsl_vector * Y = gsl_vector_from_EXTPTRSEXP(argList[0]);

			unsigned short int test = INTEGER(argList[4])[0];   // 1 = wald, 2 = score, 3 = robust

			gsl_matrix * W = NULL;      // can be NULL
			if (argList[3] != R_NilValue)
				W = gsl_matrix_from_EXTPTRSEXP(argList[3]);
			else
				W = NULL;

			gsl_matrix * tXWfixed = NULL;   // can be NULL, otherwise contains pre-computed t(X)%*%W,
			if (argList[2] != R_NilValue)
				tXWfixed = gsl_matrix_from_EXTPTRSEXP(argList[2]);
			else
				tXWfixed = NULL;

			double scorevar = NULL;      // should provide residual variance used in score test (if score used)
			if (argList[9] != R_NilValue) scorevar = REAL(argList[9])[0];

			unsigned short int minimal = INTEGER(argList[10])[0];   // only betas
			//Rprintf("minimal = %d\n", minimal);

			// ON RETURN
			gsl_vector * beta = gsl_vector_alloc(ncol);     // length ncolX, estimates
			gsl_matrix * V = gsl_matrix_alloc(ncol,ncol);        // ncolX x ncolX, var-cov matrix
			double T2 = -2.0;       // chi2 test statistic

			//gsl_matrix_fprintf(stdout,gsl_X,"before fgsl gsl_X=%f");
			T2 = fgls(Y, gsl_X, test, W, tXWfixed, WhichToTest, scorevar, beta, V, minimal);

			outdata[0] = T2;
			unsigned int k = 1;
			if (T2 < (-0.1)) {
				double zero = 0;
				//incorrect, set to NA
				for (unsigned int i = 0; i < beta->size; i++) {
					outdata[k++] = 0/zero;
				}

				for (unsigned int i = 0; i < beta->size; i++)
					for (unsigned int j = i; j < beta->size; j++) {
						outdata[k++] = 0/zero;
					}
			} else {
				//correct
				// return betas
				for (unsigned int i = 0; i < beta->size; i++) {
					outdata[k++] = gsl_vector_get(beta,i);
				}
				//return var-cov diag and off-diag
				for (unsigned int i = 0; i < beta->size; i++)
					for (unsigned int j = i; j < beta->size; j++) {
						outdata[k++] = gsl_matrix_get(V,i,j);
					}
			}
			// always return mean
			for (unsigned int j = 0; j < gsl_X->size2; j++) {
				double mn = 0;
				for (unsigned int i = 0; i < gsl_X->size1; i++)
					mn += gsl_matrix_get(gsl_X,i,j);
				mn /= (double) gsl_X->size1;
				outdata[k++] = mn;
			}
			//Rprintf("k = %d\n",k);

			//gsl_vector_free(Y);
			//gsl_matrix_free(X);
			//if (inW != NULL) gsl_matrix_free(W);
			//if (intXWfixed != NULL) gsl_matrix_free(tXWfixed);
			gsl_vector_free(beta);
			gsl_matrix_free(V);

			// END part which was in fgsl_caller



			//double * tmpout_double_ptr = REAL(tmpout);
			//Rprintf("outsize_Ncol REALSXP in fake_fglsWrapper = %d\n",outsize_Ncol);
			//for (unsigned int i = 0; i< outsize_Ncol ; i++ )
			//{
			//	outdata[i] = tmpout_double_ptr[i];
			//	Rprintf("out %d %f ",i,outdata[i]);
			//}
			//UNPROTECT(3);
			gsl_matrix_free(gsl_X);
		} else {
			outdataNcol = outsize_Ncol;
			outdataNrow = 1;
		}

	}



#ifdef __cplusplus
}
#endif
