/**





 **/

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif


	void rwth_example (
			// print the data out?
			short int * Loud,
			// dimensions of the data, explained later
			unsigned int * Length_y, unsigned int * NXRs,
			unsigned int * WidthXR, unsigned int * WidthXL,
			// parameters and data to compute M
			double * Phi, double * H2, double * Sigma2,
			// vector of length length_y
			double * y,
			// matrix of height length_y and width widthXL
			double * XL,
			// matrix of height length_y and width nXRs*widthXR
			// it contains nXRs matrices XR
			double * XRs,
			// this is output
			// 'beta's length is (widthXL + widthXR)
			// and we get a nXRs betas (so make a matrix)
			// for time-being, just fill it with mean of
			// columns of X
			double * beta //,
			// 'vcbeta' is square matrix with
			// width/height = (widthXL + widthXR)
			// for now, will skip it
			// double * vcbeta
	)
	{
		short int loud = (*Loud);
		unsigned int length_y = (*Length_y);
		unsigned int nXRs = (*NXRs);
		unsigned int widthXR = (*WidthXR);
		unsigned int widthXL = (*WidthXL);
		double h2 = (*H2);
		double sigma2 = (*Sigma2);

		Rprintf("loud = %d\n",loud);

		// M = Sigma2*(h2*phi + (1-h2)*I)
		// note that matrix is transposed in R compared to C++; i.e.
		// for R phi[2] is the element from 2nd row, 1 column,
		// while in C++ phi[2] will be second column, 1st row
		// print phi, the same view as in R
		if (loud) {
			Rprintf("phi:\n");
			for (int i = 0; i < length_y; i++) {
				for (int j = 0; j < length_y; j++) Rprintf(" %f",Phi[i*length_y+j]);
				Rprintf("\n");
			}
		}
		// access XL, print and
		// compute the column-means for XL
		if (loud) Rprintf("XL:\n");
		double mean_XL[widthXL];
		for (int i=0;i<widthXL;i++) {
			double sum = 0.0;
			for (int j=0;j<length_y;j++) {
				double data = XL[i*length_y+j];
				if (loud) Rprintf(" %f",data);
				sum += data;
			}
			if (loud) Rprintf("\n");
			mean_XL[i] = sum / (double) length_y;
		}
		// access XRs one after another,
		// print them and compute means;
		// fill in the beta
		double mean_XR[widthXR];
		unsigned int length_XR = length_y * widthXR;
		unsigned int length_beta = widthXL + widthXR;
		for (int i = 0; i < nXRs; i++) {
			if (loud) Rprintf("XR number %d:\n",i+1);
			for (int k = 0; k < widthXR; k++) {
				double sum = 0.0;
				for (int j = 0; j < length_y; j++) {
					double data = XRs[length_XR*i + k*length_y+j];
					if (loud) Rprintf(" %f",data);
					sum += data;
				}
				if (loud) Rprintf("\n");
				mean_XR[k] = sum / (double) length_y;
			}
			// fill in beta with means
			for (int j = 0; j<widthXL; j++) beta[j*nXRs + i] = mean_XL[j];
			for (int j = 0; j<widthXR; j++) beta[(j + widthXL)*nXRs + i] = mean_XR[j];
		}
	}



#ifdef __cplusplus
}
#endif


