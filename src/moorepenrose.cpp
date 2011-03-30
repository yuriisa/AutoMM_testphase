#include "moorepenrose.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include <math.h>

// for debug
#include <iostream>
int moore_penrose(gsl_matrix* pA_mat, double tol)
 {
   // compute the SVD of the variance-covariance structure
   gsl_matrix* psvd_V_mat;
   gsl_vector* psvd_S_vec;
   gsl_matrix* pUrest_mat;
   unsigned int thresh;
   
   gsl_matrix_view svd_Vrest_matview;
   unsigned int it;
   unsigned int status;
   
   if(pA_mat==0)
     GSL_ERROR ("Null input", GSL_EINVAL);
   if((pA_mat->size1==0)||(pA_mat->size1!=pA_mat->size2))
     GSL_ERROR ("Null input", GSL_EINVAL);
   
   psvd_V_mat=gsl_matrix_alloc(pA_mat->size1,pA_mat->size1);
   if(psvd_V_mat==0)
     GSL_ERROR ("Allocation failure", GSL_ENOMEM);
   
   psvd_S_vec=gsl_vector_alloc(pA_mat->size1);
   if(psvd_S_vec==0)
     {
       gsl_matrix_free(psvd_V_mat);
       GSL_ERROR ("Allocation failure", GSL_ENOMEM);
     }
   // compute svd of input   
   status=gsl_linalg_SV_decomp_jacobi(pA_mat, psvd_V_mat, psvd_S_vec);
   if(status!=GSL_SUCCESS)
     { 
       gsl_matrix_free(psvd_V_mat);
       gsl_vector_free(psvd_S_vec);
       GSL_ERROR ("Singular value decomposition error", status);
     }
   
   thresh=0;
   while((thresh<pA_mat->size2)&&(gsl_vector_get(psvd_S_vec,thresh)>CHOL_TOL))
     thresh++;
   
   // if A==0 just return A!
   if(thresh==0){
     gsl_matrix_free(psvd_V_mat);
     gsl_vector_free(psvd_S_vec);
     // std::cout<<"were are returning a null!";
     return GSL_SUCCESS;
   }
   for(it=0;it<thresh;it++)
     gsl_vector_scale(&(gsl_matrix_column(psvd_V_mat,it).vector),1.0/gsl_vector_get(psvd_S_vec,it));
   
   pUrest_mat=gsl_matrix_alloc(pA_mat->size1, thresh);
   gsl_matrix_memcpy(pUrest_mat,&(gsl_matrix_submatrix(pA_mat,0,0,pA_mat->size1, thresh).matrix));
   
   
   svd_Vrest_matview=gsl_matrix_submatrix(psvd_V_mat, 0, 0, pA_mat->size1, thresh);
   
   gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &(svd_Vrest_matview.matrix), pUrest_mat,0.0, pA_mat);
   
   gsl_matrix_free(psvd_V_mat);
   gsl_vector_free(psvd_S_vec);
   gsl_matrix_free(pUrest_mat);
   return GSL_SUCCESS;
}


// for debug
void gslprint(gsl_matrix* pmymat)
{ 
  unsigned int i,j; 
  std::cout << std::endl;
  for(i=0;i<pmymat->size1;i++)
    {
      for(j=0;j<pmymat->size2;j++)
	std::cout << gsl_matrix_get(pmymat,i,j) << ", "; 
      std::cout << std::endl;
    }
 std::cout << std::endl;
}

/*
int gsl_linalg_cholesky_invert (gsl_matrix * cholesky)
{
  unsigned int col;
  gsl_vector_view acolview;
  gsl_matrix* inverse;
 
  inverse=gsl_matrix_alloc(cholesky->size1, cholesky->size2);
  gsl_matrix_set_identity(inverse);
  for(col=0;col<cholesky->size2;col++)
    {      
      acolview=gsl_matrix_column(inverse,col);
      gsl_linalg_cholesky_svx(cholesky, &(acolview.vector));
    }
  gsl_matrix_memcpy(cholesky, inverse);
 
  gsl_matrix_free(inverse);
  return GSL_SUCCESS;
  }*/


// for debug
void gslprint(gsl_vector* pmyvec)
{
  unsigned int i;
  std::cout << std::endl;
  for(i=0;i<pmyvec->size;i++)
    std::cout << gsl_vector_get(pmyvec,i) << ", ";
  std::cout << std::endl;
}
