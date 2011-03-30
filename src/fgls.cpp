#include "fgls.h"
#include "gsl_types_extptr.h"

using namespace std;


/**

C++ code similar to the code of R's FGLS

#' Feasible GLS
#'
#' Feasible Generalised Least Squares
#'
#' @param Y dependent variable
#' @param X design matrix (including intercept, if necessary)
#' @param test test to be applied, one of 'wald', 'score' or 'robust'
#' @param whichtest which independent variables to be tested (set to 'TRUE')
#' @param W for GLS, inverse variance-covariance mRinXÄatrix, as such returned by
#' GenABEL's polygenic(...)$InvSigma, or NULL for LS
#'
#' @return List with elements 'beta' -- estimates for the regression
#' coefficients; 'V' -- variance covariance matrix for parameters
#' estimates; 'T2' -- test statistics (distributed as Chi-squared under
#' the null) for the testing of whichtest parameters; 'df' -- the number of
#' degrees of freedom of the T2 test; 'tested' -- which parameters were
#' tested with T2; 'meanColX' -- mean value of variable in columns of X;
#' 'n' -- length of Y (== height of X)
#'
#' @author Yurii Aulchenko
#'
FGLS <- function(Y,X,test="wald",whichtest = c(FALSE,rep(TRUE,dim(X)[2]-1)),W=NULL)

 **/

extern "C" {


	int invert_by_LU (gsl_matrix * intoinvert, gsl_matrix * inverted) {
		gsl_matrix * toinvert = gsl_matrix_alloc(intoinvert->size1,intoinvert->size2);
		gsl_matrix_memcpy(toinvert,intoinvert);
		gsl_permutation * perm = gsl_permutation_alloc(toinvert->size2);
		int s;
		//cout << "before gsl_linalg_LU_decomp\n";
		gsl_linalg_LU_decomp(toinvert,perm,&s);
		//cout << "before gsl_linalg_LU_invert\n";
		gsl_set_error_handler_off();
		int status = gsl_linalg_LU_invert(toinvert,perm,inverted);
		double det = gsl_linalg_LU_det(toinvert,1);
		gsl_set_error_handler (NULL);
		//cout << "invert done\n";
		gsl_matrix_free(toinvert);
		gsl_permutation_free(perm);
		//Rprintf("invert status = %d; lndet = %e\n",status,lndet);
		//if (det<0) det=det*(-1.0);
		if (status || !isfinite(det) || det==0) { // ln(1e-6)=-13.81551
			gsl_matrix_set(inverted,0,0,-999.999);
			return 1;
		} else {
			return 0;
		}
	}

	/*  SCCS @(#)cholesky2.c	5.2 10/27/98
	 ** subroutine to do Cholesky decompostion on a matrix: C = FDF'
	 **   where F is lower triangular with 1's on the diagonal, and D is diagonal
	 **
	 ** arguments are:
	 **     n         the size of the matrix to be factored
	 **     **matrix  a ragged array containing an n by n submatrix to be factored
	 **     toler     the threshold value for detecting "singularity"
	 **
	 **  The factorization is returned in the lower triangle, D occupies the
	 **    diagonal and the upper triangle is left undisturbed.
	 **    The lower triangle need not be filled in at the start.
	 **
	 **  Return value:  the rank of the matrix (non-negative definite), or -rank
	 **     it not SPD or NND
	 **
	 **  If a column is deemed to be redundant, then that diagonal is set to zero.
	 **
	 **   Terry Therneau
	 **
	 **   gsl modification -- Yurii Aulchenko
	 **
	 */

	int cholesky2_rank(gsl_matrix *matrix, double toler)
	{
		int n = matrix->size1;
		double internal_matrix[n][n];
		double temp;
		int  i,j,k;
		double eps, pivot;
		int rank;
		int nonneg;

		for (i=0;i<n;i++) for (j=0;j<n;j++)
			internal_matrix[i][j] = gsl_matrix_get(matrix,i,j);

		nonneg=1;
		eps =0;
		for (i=0; i<n; i++) {
			if (internal_matrix[i][i] > eps)  eps = internal_matrix[i][i];
			for (j=(i+1); j<n; j++)  internal_matrix[j][i] = internal_matrix[i][j];
		}
		eps *= toler;

		rank =0;
		for (i=0; i<n; i++) {
			pivot = internal_matrix[i][i];
			if (pivot < eps) {
				internal_matrix[i][i] =0;
				if (pivot < -8*eps) nonneg= -1;
			}
			else  {
				rank++;
				for (j=(i+1); j<n; j++) {
					temp = internal_matrix[j][i]/pivot;
					internal_matrix[j][i] = temp;
					internal_matrix[j][j] -= temp*temp*pivot;
					for (k=(j+1); k<n; k++) internal_matrix[k][j] -= temp*internal_matrix[k][i];
				}
			}
		}
		return(rank * nonneg);
	}


	int invert_by_chol (gsl_matrix * intoinvert, gsl_matrix * inverted) {
		gsl_matrix_memcpy(inverted,intoinvert);
		//cout << "before gsl_linalg_cholesky_decomp\n";
		//cout << "before gsl_linalg_LU_invert\n";
		gsl_set_error_handler_off();
		int status = gsl_linalg_cholesky_decomp(inverted);
		gsl_set_error_handler (NULL);
		gsl_linalg_cholesky_invert(inverted);
		/**
		int rank = cholesky2_rank(intoinvert,1e-6); //,1.5e-12);
		if (rank < intoinvert->size1) {
			gsl_matrix_set(inverted,0,0,-999.999);
			return 1;
		} else {
			return 0;
		}
		 **/
		return(status);
	}


	// returns chi2 test statistic
	double fgls(gsl_vector * Y, gsl_matrix * X,
			unsigned short int test, // 1 = wald, 2 = score, 3 = robust
			gsl_matrix * W,  // can be NULL
			gsl_matrix * tXWfixed,   // can be NULL, otherwise contains pre-computed t(X)%*%W,
			// where only elements t(X[,!whichtest]) %*% W are fixed and
			// pre-computed
			unsigned int WhichToTest, // at what column of X 'fixed' part starts (the 1st column is 1!)
			double scorevar,   // can be negative, otherwise residual variance used in score test
			// ON RETURN
			gsl_vector * beta, // length ncolX, estimates
			gsl_matrix * V,    // ncolX x ncolX, var-cov matrix
			unsigned short int minimal // only betas to be returned
	)
	{
		//gsl_vector_fprintf(stdout,Y,"Y=%f");
		//gsl_matrix_fprintf(stdout,X,"X=%f");
		int tmpStatus = 0, Status = 0;
		double T2 = -2.0; // return chi2 test stat value or -2.0 if failed; -1.0 means could not invert


		//cout << "minimal = " << minimal << "\n";
		/**
  if (any(is.na(Y)) || any(is.na(X))) {
    warning("missing data points in Y or X, dropping")
    IScomplete <- !is.na(Y) & complete.cases(X)
    Y <- Y[IScomplete]
    X <- X[IScomplete,]
    if (!is.null(W)) W <- W[IScomplete,IScomplete]
  }

  NANs have to be dealt with before the call to the function
  because GSL does not support NaN

		 **/


		/**
  # precompute tXW
  if (is.null(W)) {
    tXW <- t(X)
  } else {
    tXW <- t(X) %*% W
  }
		 **/
		gsl_matrix * tXW = gsl_matrix_alloc(X->size2,X->size1);

		if (W != NULL) {
			if (tXWfixed == NULL) {
				/**
        gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB,
            double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
        These functions compute the matrix-matrix product and sum
        C = \alpha op(A) op(B) + \beta C
        where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans
        and similarly for the parameter TransB
        t(X) %*% W
				 **/
				gsl_blas_dgemm (CblasTrans, CblasNoTrans,
						1.0, X, W, 0.0, tXW);
				//Rprintf("X(2-1,1-1)=%f\n",gsl_matrix_get(X,1,0));
				//Rprintf("W(2-1,1-1)=%f\n",gsl_matrix_get(W,1,0));
				//Rprintf("tXW(2-1,1-1)=%f\n",gsl_matrix_get(tXW,1,0));
				//cout << "*** tXW:\n"; gsl_matrix_fprintf(stdout,tXW,"%10f");

			} else {
				//cout << "before\n";
				gsl_matrix_memcpy(tXW,tXWfixed);
				//cout << "tXWfixed in:\n";gsl_matrix_fprintf(stdout,tXWfixed,"%f");
				//cout << "1:\n";gsl_matrix_fprintf(stdout,tXW,"%f");
				//Rprintf("%d %d %d %d\n",WhichToTest-1,0,(tXW->size1-WhichToTest+1),tXW->size2);
				gsl_matrix_view tmp_tXW_unfixed =
						gsl_matrix_submatrix(tXW,WhichToTest-1,0,(tXW->size1-WhichToTest+1),tXW->size2);
				//cout << "after\n";
				gsl_matrix_view tmp_X_unfixed =
						gsl_matrix_submatrix(X,0,WhichToTest-1,X->size1,(X->size2-WhichToTest+1));
				gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0,
						&(tmp_X_unfixed.matrix), W, 0.0, &(tmp_tXW_unfixed.matrix));
				//cout << "2:\n";gsl_matrix_fprintf(stdout,tXW,"%f");
			}
		} else {
			//if (tXWfixed == NULL) {
			gsl_vector * tmp = gsl_vector_alloc(X->size1);
			for (unsigned int j = 0; j < (X->size2); j++) {
				gsl_matrix_get_col (tmp, X, j);
				gsl_matrix_set_row (tXW, j, tmp);
			}
			gsl_vector_free(tmp);
			//} else {

			//}
		}


		//Rprintf ("*** tXW (%d,%d) :\n",tXW->size1,tXW->size2); gsl_matrix_fprintf(stdout,tXW,"%f");

		/**
  # estimate beta
  XpWXm1 <- ginv(tXW %*% X)
		 **/

		gsl_matrix * XpWX = gsl_matrix_alloc(X->size2,X->size2);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
				1.0, tXW, X, 0.0, XpWX);

		gsl_matrix * XpWXm1 = gsl_matrix_alloc(X->size2,X->size2);
		tmpStatus = invert_by_LU (XpWX,XpWXm1);
		if (tmpStatus) {
			//cout << "tmpStatus " << tmpStatus << " returned by invert XpWx\n";
			Status = 1;
		}

		gsl_matrix_free(XpWX);
		/**
    gsl_permutation * perm = gsl_permutation_alloc(X->size2);
    int s;
    gsl_linalg_LU_decomp(XpWX,perm,&s);
    gsl_linalg_LU_invert(XpWX,perm,XpWXm1);
		 **/

		/**
  tXWY <- tXW %*% Y
  betaW <- XpWXm1 %*% tXWY
		 **/
		gsl_vector * tXWY = gsl_vector_alloc(X->size2);
		gsl_blas_dgemv (CblasNoTrans, 1.0, tXW, Y, 0.0, tXWY);
		gsl_blas_dgemv (CblasNoTrans, 1.0, XpWXm1, tXWY, 0.0, beta);

		gsl_vector_free(tXWY);
		//Rprintf("beta");
		//gsl_vector_fprintf(stdout,beta,"%f");


		gsl_matrix_memcpy(V,XpWXm1);

		if (minimal) {
			//cout << "oooo";
			return(0.01);
		}

		/**
  # estimate V
  if (test=="wald") {
    YmXB <- Y - X %*% betaW
    if (is.null(W)) {
      sigma2 <- ((t(YmXB) %*% YmXB)/(dim(X)[1] - dim(X)[2]))[1,1]
    } else {
      sigma2 <- ((t(YmXB) %*% W %*% YmXB)/(dim(X)[1] - dim(X)[2]))[1,1]
    }
    V <- sigma2*XpWXm1
		 **/
		if (test == 1) {
			gsl_vector * YmXB = gsl_vector_alloc(X->size1);
			gsl_vector_memcpy(YmXB,Y);
			gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, -1.0, YmXB);
			//gsl_vector_fprintf(stdout,YmXB,"%f");
			double sigma2 = 0;
			if (W == NULL) {
				gsl_vector_mul (YmXB,YmXB);
			} else {
				gsl_vector * WYmXB = gsl_vector_alloc(X->size1);
				gsl_vector_memcpy(WYmXB,YmXB);
				gsl_blas_dgemv (CblasNoTrans, 1.0, W, YmXB, 0.0, WYmXB);
				gsl_vector_mul (YmXB,WYmXB);
				gsl_vector_free(WYmXB);
			}
			for (unsigned int i = 0; i < YmXB->size; i++) sigma2 += gsl_vector_get(YmXB,i);

			gsl_vector_free(YmXB);

			sigma2 /= (X->size1 - X->size2);
			//Rprintf("YmXB");
			//gsl_vector_fprintf(stdout,YmXB,"%f");
			gsl_matrix_scale(V,sigma2);
			//Rprintf("sigma2=%f; V = \n",sigma2);
			//gsl_matrix_fprintf(stdout,V,"%f");
			/**
  } else if (test=="score") {
    Xm <- X[,!whichtest,drop=FALSE]
    if (is.null(W)) {
      bt <- ginv(t(Xm) %*% Xm) %*% (t(Xm) %*% Y)
    } else {
      bt <- ginv(t(Xm) %*% W %*% Xm) %*% (t(Xm) %*% W %*% Y)
    }
    sigma2 <- sum(as(Y - Xm %*% bt,"vector")^2,na.rm=T)/(dim(Xm)[1]-dim(Xm)[2])
    V <- sigma2*XpWXm1
			 **/
		} else if (test == 2) {

			if (scorevar > 0 ) {
				gsl_matrix_scale(V,scorevar);
			} else {
				// here should compute score var
				gsl_matrix_free(XpWXm1);
				gsl_matrix_free(tXW);
				cout << "score test, but missing scorevar\n";
				return(-2.0);
			}

			/**
  } else if (test=="robust") {
		YmXB <- Y - X %*% betaW
		if (is.null(W)) {
			sigma2vec <- as.vector(YmXB * YmXB)
		} else {
			sigma2vec <- as.vector(as.vector(t(YmXB) %*% W) * YmXB)
		}
		V <- XpWXm1 %*% (tXW %*% diag(sigma2vec) %*% X) %*% XpWXm1
			 **/
		} else if (test == 3) {

			// the same code as in 'wald'
			gsl_vector * YmXB = gsl_vector_alloc(Y->size);
			gsl_vector_memcpy(YmXB,Y);
			gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, -1.0, YmXB);

			if (W == NULL) {
				gsl_vector_mul (YmXB,YmXB);
			} else {
				gsl_vector * WYmXB = gsl_vector_alloc(Y->size);
				gsl_vector_memcpy(WYmXB,YmXB);
				gsl_blas_dgemv (CblasNoTrans, 1.0, W, YmXB, 0.0, WYmXB);
				gsl_vector_mul (YmXB,WYmXB);
				gsl_vector_free(WYmXB);
			}

			// now the difference with 'wald' comes
			// V <- XpWXm1 %*% (tXW %*% diag(sigma2vec) %*% X) %*% XpWXm1
			gsl_matrix * tXW_diag_sigma2vec = gsl_matrix_alloc(tXW->size1,tXW->size2);
			gsl_vector * tmp = gsl_vector_alloc(tXW->size2);
			gsl_matrix * tXW_diag_sigma2vec_X = gsl_matrix_alloc(tXW->size1,tXW->size1);
			gsl_matrix * XpWXm1_tXW_diag_sigma2vec_X = gsl_matrix_alloc(tXW->size1,tXW->size1);
			for (unsigned int i = 0; i < tXW->size1; i++) {
				gsl_matrix_get_row(tmp,tXW,i);
				gsl_vector_mul(tmp,YmXB);
				gsl_matrix_set_row(tXW_diag_sigma2vec,i,tmp);
			}
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
					1.0, tXW_diag_sigma2vec, X, 0.0, tXW_diag_sigma2vec_X);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
					1.0, XpWXm1, tXW_diag_sigma2vec_X, 0.0, XpWXm1_tXW_diag_sigma2vec_X);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
					1.0, XpWXm1_tXW_diag_sigma2vec_X, XpWXm1, 0.0, V);
			gsl_matrix_transpose(V);

			//gsl_matrix_memcpy(V,tXW_diag_sigma2vec_X);

			// here should compute all the staff

			gsl_matrix_free(tXW_diag_sigma2vec);
			gsl_matrix_free(XpWXm1_tXW_diag_sigma2vec_X);
			gsl_matrix_free(tXW_diag_sigma2vec_X);
			gsl_vector_free(YmXB);
			gsl_vector_free(tmp);

			//gsl_matrix_fprintf(stdout,V," V = %f");
			// on error
			// gsl_matrix_free(XpWXm1);
			// gsl_matrix_free(tXW);
			// return(-2.0);

			/**
  } else {
    stop("test not recognised")
  }
			 **/
		} else {
			// option not recognised
			gsl_matrix_free(XpWXm1);
			gsl_matrix_free(tXW);
			cout << "test not recognised\n";
			return(-2.0);
		}



		/**
  # do the test
  if (sum(whichtest)>1) {
    Vinv <- ginv(V)
    Vbg <- Vinv[whichtest,whichtest]-Vinv[whichtest,!whichtest] %*%
           ginv(Vinv[!whichtest,!whichtest]) %*% Vinv[!whichtest,whichtest]
		 **/
		//fprintf(stdout,"WhichToTest=%d\n",WhichToTest);
		if (WhichToTest < X->size2) {
			//Rprintf("VV = \n");
			//gsl_matrix_fprintf(stdout,V,"%f");
			gsl_matrix * Vm1 = gsl_matrix_alloc(X->size2,X->size2);
			tmpStatus = invert_by_LU (V,Vm1);
			if (tmpStatus) {
				//cout << "tmpStatus " << tmpStatus << " returned by invert V\n";
				//gsl_matrix_fprintf(stdout,V,"V=%f");
				//gsl_matrix_fprintf(stdout,Vm1,"Vm1=%f");
				Status = 1;
			}

			gsl_matrix_view Vm1ww = gsl_matrix_submatrix(Vm1,
					WhichToTest-1,WhichToTest-1,
					(Vm1->size1-WhichToTest+1),(Vm1->size1-WhichToTest+1));
			gsl_matrix_view Vm1wNw = gsl_matrix_submatrix(Vm1,
					WhichToTest-1,0,
					(Vm1->size1-WhichToTest+1),(WhichToTest-1));
			gsl_matrix_view Vm1NwNw = gsl_matrix_submatrix(Vm1,
					0,0,
					(WhichToTest-1),(WhichToTest-1));
			gsl_matrix * Vm1NwNwm1 = gsl_matrix_alloc((WhichToTest-1),(WhichToTest-1));
			tmpStatus = invert_by_LU (&(Vm1NwNw.matrix),Vm1NwNwm1);
			if (tmpStatus) {
				//cout << "tmpStatus " << tmpStatus << " returned by invert &(Vm1NwNw.matrix)\n";
				Status = 1;
			}
			gsl_matrix_view Vm1Nww = gsl_matrix_submatrix(Vm1,
					0,WhichToTest-1,
					(WhichToTest-1),(Vm1->size1-WhichToTest+1));

			//Rprintf("AAA\n");
			gsl_matrix * res = gsl_matrix_alloc((Vm1->size1-WhichToTest+1),(WhichToTest-1));
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0,
					&(Vm1wNw.matrix), Vm1NwNwm1, 0.0, res);
			//Rprintf("AAA %d %d\n",Vm1NwNwm1->size1,Vm1NwNwm1->size2);
			//gsl_matrix_fprintf(stdout,res,"%f");
			gsl_matrix * res1 = gsl_matrix_alloc((Vm1->size1-WhichToTest+1),(Vm1->size1-WhichToTest+1));
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0,
					res, &(Vm1Nww.matrix), 0.0, res1);
			//gsl_matrix_fprintf(stdout,res1,"%f");
			gsl_matrix_sub(res1,&(Vm1ww.matrix));
			gsl_matrix_scale(res1,-1.0);

			/**
			 *     T2 <- t(betaW[whichtest]) %*% Vbg %*% betaW[whichtest]
			 */
			gsl_vector * Vbgbeta = gsl_vector_alloc(res1->size1);
			gsl_vector_view testBeta = gsl_vector_subvector (beta,(WhichToTest-1),res1->size1);
			gsl_blas_dgemv (CblasNoTrans, 1.0, res1, &(testBeta.vector), 0.0, Vbgbeta);
			T2 = 0;
			for (unsigned int i = 0; i< (Vbgbeta->size); i++)
				T2 += gsl_vector_get(Vbgbeta,i)*gsl_vector_get(&(testBeta.vector), i );

			gsl_matrix_free(Vm1);
			gsl_matrix_free(Vm1NwNwm1);
			gsl_matrix_free(res);
			gsl_matrix_free(res1);
			gsl_vector_free(Vbgbeta);


			/**
  } else if (sum(whichtest) == 1) {
    T2 <- betaW[whichtest]^2/diag(V)[whichtest]
			 **/
		} else if (WhichToTest == X->size2) {
			T2 =
					gsl_vector_get(beta,WhichToTest-1)*gsl_vector_get(beta,WhichToTest-1)
					/gsl_matrix_get(V,WhichToTest-1,WhichToTest-1);
			/**
  } else {
    stop("unreachable statement")
  }
			 **/
		} else {
			gsl_matrix_free(XpWXm1);
			gsl_matrix_free(tXW);
			cout << "WhichToTest > X->size2\n";
			return(-2.0);
		}

		/**

  rownames(betaW) <- colnames(X)
  dimnames(V) <- list(colnames(X),colnames(X))
  out <- list(beta=betaW,V=V,T2=T2,df=sum(whichtest),tested=whichtest,
      meanColX = apply(X,FUN=mean,MAR=2), n = dim(X)[1])
  class(out) <- "FGLS"
  out
		 **/

		gsl_matrix_free(XpWXm1);
		gsl_matrix_free(tXW);

		//cout << "tmpStatus = " << tmpStatus << endl;
		if (Status) {
			//cout << "Status -1.0\n";
			return(-1.0);
		} else {
			return(T2);
		}

	}



	// end 'extern "C"'

}
