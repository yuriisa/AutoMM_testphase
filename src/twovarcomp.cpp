#include "twovarcomp.h"
#include "exception.h"
#include "moorepenrose.h"
#include <iostream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <gsl/gsl_min.h>

// should this be a friend of the class?
struct dev_pars
{
  TwoVarCompModel* pTVCM;
  // why does DevWorkspace need to be a subclass of TVCM?
  TwoVarCompModel::DevWorkspace* pdev_work;
};


TwoVarCompModel::TwoVarCompModel(gsl_vector* presponse_vec, gsl_matrix* pcovar_mat, gsl_matrix* pvar1_mat, gsl_matrix* pvar2_mat, gsl_matrix* pincid1_mat, gsl_matrix* pincid2_mat, double* pnu_naught, double* pgammasq_naught)
{
  AmITheDaddy=1;
  
  // Check for unexpected null pointer (other arguments are checked by respective 'Set' functions)
  if(presponse_vec==0)
    THROW_INVALID_ARGUMENT(EXC_CSTR_RESPNULL, "unexpected null pointer.");
  
  // Get sample size from response vector 
  n=(unsigned int) presponse_vec->size;

  // Check we have a sample of at least 1
  if(n==0)
    THROW_INVALID_ARGUMENT(EXC_CSTR_NO_SAMPLE, "sample size of zero.");
  
  // Allocate memory for member data
  py_vec=gsl_vector_alloc(n); 
  // Check allocation went ok
  if(this->py_vec==0)
    THROW_BAD_ALLOC(EXC_CSTR_RESPALLOC, "failure to allocate py_vec on the heap.");
  
  // Set member pointers to null to indicate memory needs allocating 
  this->pw_mat=0;
  this->ppP_mat[0]=0;
  this->ppP_mat[1]=0;
  this->ppvar_mat[0]=0;
  this->ppvar_mat[1]=0;
  this->ppincid_mat[0]=0;
  this->ppincid_mat[1]=0;
  this->ppfullvar_mat[0]=0;
  this->ppfullvar_mat[1]=0;
  this->ppD_vec[0]=0;
  this->ppD_vec[1]=0;
  this->ppvR_mat[0]=0;
  this->ppvR_mat[1]=0;
  this->px_mat=0;  
  this->ppPVP_mat[0]=0;
  this->ppPVP_mat[1]=0;
  this->ppPVPy_vec[0]=0;
  this->ppPVPy_vec[1]=0;
  this->ppRvP_mat[0]=0;
  this->ppRvP_mat[1]=0;
  this->ppRvPy_vec[0]=0;
  this->ppRvPy_vec[1]=0;
  this->ppDminus1_vec[0]=0;
  this->ppDminus1_vec[1]=0;
  this->ppLambda_vec[0]=0;
  this->ppLambda_vec[1]=0;
  this->ppRvPx_mat[0]=0;
  this->ppRvPx_mat[1]=0;
  this->ppPVPx_mat[0]=0;
  this->ppPVPx_mat[1]=0;
  this->ppxPVPx_mat[0]=0;
  this->ppxPVPx_mat[1]=0;
  this->ppxPVPy_vec[0]=0;
  this->ppxPVPy_vec[1]=0;
 
  

  // Copy data into object
  gsl_vector_memcpy(this->py_vec,presponse_vec);
  
 // freed
  ppPVP_mat[0]=gsl_matrix_alloc(n,n);
   //freed
  ppPVPy_vec[0]=gsl_vector_alloc(n); 
   //freed
  ppRvP_mat[0]=gsl_matrix_alloc(n,n);
   //freed
  ppRvPy_vec[0]=gsl_vector_alloc(n);
  // freed
  ppvR_mat[0]=gsl_matrix_alloc(n,n);
  //freed
  ppLambda_vec[0]=gsl_vector_alloc(n);
  // freed
  ppD_vec[0]=gsl_vector_alloc(n);
  // need we store D at all?
  //freed
  ppDminus1_vec[0]=gsl_vector_alloc(n);
 
  SetCovariates(pcovar_mat);
 
  SetVar(pvar1_mat, pincid1_mat,0);
  SetVar(pvar2_mat, pincid2_mat,1);
  SetScalarVarPrior(*pnu_naught, *pgammasq_naught);

  
  TakeYourMarks(0);
  
  // this needs to be computed somewhere
  rankCovar=pw_mat->size2;
  
}

// Copy constructor should this be const?
TwoVarCompModel::TwoVarCompModel(TwoVarCompModel & daddy)
{
  AmITheDaddy=0;
  
  // cant i just copy the object in one go instead of all these pointers?
  n=daddy.n;
  // Set member pointers to null to indicate memory needs allocating 
  this->pw_mat=daddy.pw_mat;
  this->ppvar_mat[0]=daddy.ppvar_mat[0];
  this->ppvar_mat[1]=daddy.ppvar_mat[1];
  this->ppincid_mat[0]=daddy.ppincid_mat[0];
  this->ppincid_mat[1]=daddy.ppincid_mat[1];
  this->ppfullvar_mat[0]=daddy.ppfullvar_mat[0];
  this->ppfullvar_mat[1]=daddy.ppfullvar_mat[1];
  this->ppD_vec[0]=daddy.ppD_vec[0];
  this->ppD_vec[1]=daddy.ppD_vec[1];
  this->ppvR_mat[0]=daddy.ppvR_mat[0];
  this->ppvR_mat[1]=daddy.ppvR_mat[1];
  this->py_vec=daddy.py_vec;
  this->ppPVP_mat[0]=daddy.ppPVP_mat[0];
  this->ppPVP_mat[1]=daddy.ppPVP_mat[1];
  this->ppPVPy_vec[0]=daddy.ppPVPy_vec[0];
  this->ppPVPy_vec[1]=daddy.ppPVPy_vec[1];
  this->ppRvPy_vec[0]=daddy.ppRvPy_vec[0];
  this->ppRvPy_vec[1]=daddy.ppRvPy_vec[1];
  this->ppDminus1_vec[0]=daddy.ppDminus1_vec[0];
  this->ppDminus1_vec[1]=daddy.ppDminus1_vec[1];
  this->ppP_mat[0]=daddy.ppP_mat[0];
  this->ppP_mat[1]=daddy.ppP_mat[1];
  this->ppLambda_vec[0]=daddy.ppLambda_vec[0];
  this->ppLambda_vec[1]=daddy.ppLambda_vec[1];
  this->ppRvP_mat[0]=daddy.ppRvP_mat[0];
  this->ppRvP_mat[1]=daddy.ppRvP_mat[1];
  // note here actually copying values is this what we want?
  this->pyPVPy[0]=daddy.pyPVPy[0];
  this->pyPVPy[1]=daddy.pyPVPy[1];

  this->gamma0sq=daddy.gamma0sq;
  this->nu0=daddy.nu0;
   
  ppRvPx_mat[0]=0; 
  ppRvPx_mat[1]=0;
  ppPVPx_mat[0]=0; 
  ppPVPx_mat[1]=0;
  ppxPVPx_mat[0]=0;
  ppxPVPx_mat[1]=0;
  ppxPVPy_vec[0]=0;  
  ppxPVPy_vec[1]=0;
  

  // this needs to be computed somewhere
  rankCovar=daddy.rankCovar;
}


void TwoVarCompModel::SetCovariates(gsl_matrix* pcovar_mat)
{
  if(pcovar_mat==0)
    {
      ClearUp();
      THROW_INVALID_ARGUMENT(EXC_SETCOV_COVARNULL, "unexpected null pointer.");
    }

  if(pcovar_mat->size1!=n)
    {
      ClearUp();
      THROW_INVALID_ARGUMENT(EXC_SETCOV_COVARDIM, "vector or matrix of unexpected dimension.");
    }
  
  pw_mat=pcovar_mat;
  /*// Check for unexpected null pointers
    

  if(pcovar_mat->size2>0)
    {
      pw_mat=gsl_matrix_alloc(n,pcovar_mat->size2);
      // check if allocation went okay.
      if(pw_mat==0)
	{
	  ClearUp();
	  THROW_BAD_ALLOC(EXC_SETCOV_COVARALLOC, "failure to allocate pw_mat on the heap.");
	}
	gsl_matrix_memcpy(pw_mat, pcovar_mat);
	}*/
}

void TwoVarCompModel::SetVar1(gsl_matrix* pvar_mat, gsl_matrix* pincid_mat)
{
  SetVar(pvar_mat, pincid_mat, 0);
  TakeYourMarks(0);
}

void TwoVarCompModel::SetVar2(gsl_matrix* pvar_mat, gsl_matrix* pincid_mat)
{
  SetVar(pvar_mat, pincid_mat, 1);
  TakeYourMarks(1);
}

void TwoVarCompModel::SetVar(gsl_matrix* pvar_mat, gsl_matrix* pincid_mat, unsigned int k)
{
  // Check for unexpected null pointers
  if((pvar_mat==0)||(pincid_mat==0))
    {
      ClearUp();
      THROW_INVALID_ARGUMENT(EXC_SETVAR_NULLPTR, "unexpected null pointer.");
    }
  // Check dimensions are consistent 
  if((pincid_mat->size1!=n)||(pvar_mat->size1!=pvar_mat->size2)||(pincid_mat->size2!=pvar_mat->size1))
    {
      ClearUp();
      THROW_INVALID_ARGUMENT(EXC_SETVAR_DIM, "vector or matrix of unexpected dimension.");
    }

  // Delete any previous variance objects
  if(ppvar_mat[k]!=0)
    {
      delete ppP_mat[k];
      delete ppvar_mat[k];
      delete ppincid_mat[k];
      delete ppfullvar_mat[k];
    }

  ppvar_mat[k]=gsl_matrix_alloc(pvar_mat->size1, pvar_mat->size1);
  ppincid_mat[k]=gsl_matrix_alloc(n,pvar_mat->size1);
  ppfullvar_mat[k]=gsl_matrix_alloc(n,n);

  // check if allocation went okay.
  if((ppvar_mat[k]==0)||(ppincid_mat[k]==0)||(ppfullvar_mat[k]==0))
    {
      ClearUp();
      THROW_BAD_ALLOC(EXC_SETVAR_ALLOC, "failure to allocate one of ppvar_mat[k], ppincid_mat[k], ppfullvar_mat[k] on the heap.");
    }
 
  gsl_matrix_memcpy(ppincid_mat[k], pincid_mat);
  gsl_matrix_memcpy(ppvar_mat[k], pvar_mat);
  
  gsl_matrix* pa_mat;
  
  // Compute variance-covariance structure from the two matrices
  pa_mat=gsl_matrix_alloc(n,pvar_mat->size1);
  
  if(pa_mat==0)
    {
      ClearUp();
      THROW_BAD_ALLOC(EXC_SETVAR_AMATALLOC, "failure to allocate pa_mat on the heap.");
    }
  
  gsl_blas_dsymm (CblasLeft, CblasUpper, 1.0, ppincid_mat[k], ppvar_mat[k], 0.0, pa_mat);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, pa_mat, ppincid_mat[k], 0.0, ppfullvar_mat[k]);
  gsl_matrix_free(pa_mat);
    // this temporary to fix wKw =0 for k=1
  if(k==0)
    {
	  
      gsl_matrix* pK_mat;
      pK_mat=gsl_matrix_alloc(n,n);
      
      if(pK_mat==0)
	{
	  ClearUp();
	  THROW_BAD_ALLOC(EXC_SETVAR_KMATALLOC, "failure to allocate pK_mat on the heap.");
	}
      
      gsl_matrix_memcpy(pK_mat, ppfullvar_mat[k]);
   
      if(moore_penrose(pK_mat, CHOL_TOL)!=GSL_SUCCESS)
	{
	  ClearUp();
	  gsl_matrix_free(pK_mat); 
	  THROW(EXC_SETVAR_PSEUDOERROR, "GSL error, failure to compute a pseudo-inverse.");
	}
      
      ppP_mat[k]=gsl_matrix_alloc(n,n);
      if(ppP_mat[k]==0)
	{
	  ClearUp();
	  gsl_matrix_free(pK_mat); 
	  THROW_BAD_ALLOC(EXC_SETVAR_PMATALLOC, "failure to allocate ppP_mat[k] on the heap.");
	} 
      gsl_matrix_memcpy(ppP_mat[k], pK_mat);
      
      
      // if there are covariates compute projection part
      if(pw_mat)
	{ 
	  // Calculate 
	  // what if wK=0?
	  gsl_matrix* pwK_mat, *pwKw_mat, *pKwIN_mat;
	  pwK_mat=gsl_matrix_alloc(pw_mat->size2,n);
	  pwKw_mat=gsl_matrix_alloc(pw_mat->size2,pw_mat->size2);
	  pKwIN_mat=gsl_matrix_alloc(n,pw_mat->size2);
	  
	  if((pwK_mat==0)||(pwKw_mat==0)||(pKwIN_mat==0))
	    {
	      if(pwK_mat)
		gsl_matrix_free(pwK_mat);
	      if(pwKw_mat)
		gsl_matrix_free(pwKw_mat);
	      if(pKwIN_mat)
		gsl_matrix_free(pKwIN_mat);
	      ClearUp();
	      gsl_matrix_free(pK_mat); 
	      THROW_BAD_ALLOC(EXC_SETVAR_KWMATALLOC, "failure to allocate pwK, pwKw_mat or pKwIN_mat on the heap.");
	    } 
	  
	  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, pw_mat, pK_mat,0.0, pwK_mat);
	  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, pwK_mat, pw_mat,0.0, pwKw_mat);
	  
	 
	  // !!we need to check if wKw is invertible here
	  gsl_linalg_cholesky_decomp(pwKw_mat);
	  gsl_linalg_cholesky_invert(pwKw_mat);
	


	  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, pwK_mat, pwKw_mat,0.0, pKwIN_mat);
	  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, pKwIN_mat, pwK_mat, 1.0, ppP_mat[k]);
	 
	  gsl_matrix_free(pwK_mat); 
	  gsl_matrix_free(pwKw_mat);
	  gsl_matrix_free(pKwIN_mat);
	}
      gsl_matrix_free(pK_mat); 
    }
}


void TwoVarCompModel::TakeYourMarks(unsigned int k)
{

  unsigned it;
  //need to add throws 
  // need to check whether stuff needs deleting, lambda etc
  // ppvR_mat[k]=gsl_matrix_alloc(n,n);
  //ppD_vec[k]=gsl_vector_alloc(n);
  gsl_matrix* pvLambdaSqrt_mat=gsl_matrix_alloc(n,n);
  gsl_matrix* pBDB_mat=gsl_matrix_alloc(n,n);
  gsl_matrix* pB_mat=gsl_matrix_alloc(n,n);
  gsl_matrix* pPvLambdaSqrt_mat=gsl_matrix_alloc(n,n);
  gsl_vector_view avec_view;
  gsl_eigen_symmv_workspace* peig_workspace=gsl_eigen_symmv_alloc(n);
  
  //ppvR_mat[k]=gsl_matrix_alloc(n,n);
  //ppLambda_vec[k]=gsl_vector_alloc(n);
  //ppD_vec[k]=gsl_vector_alloc(n);
  // need we store D at all?
  // ppDminus1_vec[k]=gsl_vector_alloc(n);
 
  gsl_eigen_symmv(ppfullvar_mat[1-k], ppLambda_vec[k], pvLambdaSqrt_mat, peig_workspace);
  /*std::cout<<"v"<<std::endl;
  gslprint(pvLambdaSqrt_mat);
  std::cout<<"ppLambda_vec[k]"<<std::endl;
  gslprint(ppLambda_vec[k]);*/
  for(it=0;it<n;it++)
    {
      if(gsl_vector_get(ppLambda_vec[k],it)<0.0)
	gsl_vector_set(ppLambda_vec[k],it, 0.0);
    }
  for(it=0;it<n;it++)
    {
      // abs to prevent taking sqrt of small -ve numbers 
      avec_view=gsl_matrix_column(pvLambdaSqrt_mat, it);
      gsl_vector_scale(&(avec_view.vector), sqrt(gsl_vector_get(ppLambda_vec[k],it)));
    }
  // CALCULATING Pv twice (here and below) - need to fix this.
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ppP_mat[k], pvLambdaSqrt_mat, 0.0,pPvLambdaSqrt_mat);
  
 
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, pvLambdaSqrt_mat, pPvLambdaSqrt_mat, 0.0,pBDB_mat);
 
  /*
   std::cout<<"pBDB"<<std::endl;
   gslprint(pBDB_mat);
  */
 
  gsl_eigen_symmv(pBDB_mat, ppD_vec[k],  pB_mat, peig_workspace);
  gsl_eigen_symmv_free(peig_workspace); 
  /*std::cout<<"pP_mat"<<std::endl;
  gslprint(ppP_mat[k]);
  */
  
  // do we need to sort them? or does this waste time
  gsl_eigen_symmv_sort (ppD_vec[k], pB_mat, GSL_EIGEN_SORT_VAL_DESC);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, pvLambdaSqrt_mat, pB_mat, 0.0,ppvR_mat[k]);
 
  /*
  std::cout<<"pvLambdaSqrt_mat"<<std::endl;
   gslprint(pvLambdaSqrt_mat);
   std::cout<<"pB_mat"<<std::endl;
   gslprint(pB_mat);
   std::cout<<"ppvR_mat[k]"<<std::endl;
   gslprint(ppvR_mat[k]);
   std::cout<<"ppD_vec[k]"<<std::endl;
   gslprint(ppD_vec[k]);
  */
  gsl_matrix_free(pvLambdaSqrt_mat);
  gsl_matrix_free(pPvLambdaSqrt_mat);
  gsl_matrix_free(pBDB_mat);
   gsl_matrix_free(pB_mat);
 
 
  // do the precalculations on the 
  gsl_vector_memcpy(ppDminus1_vec[k],ppD_vec[k]);
  gsl_vector_add_constant(ppDminus1_vec[k], -1.0);
  // do we need to store R in the class etc? need to check everything is essential
  //   gslprint(ppDminus1_vec[k]);

  
  gsl_matrix* pVP_mat=gsl_matrix_alloc(n,n);
  //  ppPVP_mat[k]=gsl_matrix_alloc(n,n);
  //ppPVPy_vec[k]=gsl_vector_alloc(n); 
  //ppRvP_mat[k]=gsl_matrix_alloc(n,n);
  //ppRvPy_vec[k]=gsl_vector_alloc(n);
  

  // use symmetric routines faster?
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ppfullvar_mat[k], ppP_mat[k], 0.0,pVP_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ppP_mat[k], pVP_mat, 0.0, ppPVP_mat[k]);
  gsl_blas_dsymv(CblasUpper, 1.0, ppPVP_mat[k], py_vec, 0.0, ppPVPy_vec[k]);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, ppvR_mat[k], ppP_mat[k], 0.0, ppRvP_mat[k]);
  gsl_blas_dgemv(CblasNoTrans, 1.0, ppRvP_mat[k], py_vec, 0.0, ppRvPy_vec[k]);
  gsl_blas_ddot(ppPVPy_vec[k], py_vec, &(pyPVPy[k]));
 gsl_matrix_free(pVP_mat);
  /*
  std::cout<<"ppPVP_mat[k]"<<std::endl;
  gslprint(ppPVP_mat[k]);
  std::cout<<"ppPVPy_vec[k]"<<std::endl;
  gslprint(ppPVPy_vec[k]);
  std::cout<<"ppRvP_mat[k]"<<std::endl;
  gslprint(ppRvP_mat[k]);
  std::cout<<"ppRvPy_vec[k]"<<std::endl;
  gslprint(ppRvPy_vec[k]);
  */

}

void TwoVarCompModel::SetExplan(gsl_matrix* pexplan_mat)
{
  // should really check if we need to reallocate by checking for size match.
  px_mat=pexplan_mat;//gsl_matrix_alloc(pexplan_mat->size1,pexplan_mat->size2);
  //gsl_matrix_memcpy(px_mat, pexplan_mat);
  UpdatePostNewExplan(0);
  
}

void TwoVarCompModel::SetExplan(gsl_vector* pexplan_vec)
{
  /* should we use gsl_matrix_const_view_vector ?*/
  x_matview=gsl_matrix_view_vector(pexplan_vec, pexplan_vec->size, 1);
  px_mat=&(x_matview.matrix);
  UpdatePostNewExplan(0);
}
 
void TwoVarCompModel::UpdatePostNewExplan(unsigned int k)
{  
  // need to check the sizes are the same before we reallocate
  if(!ppRvPx_mat[k])
    { 
      ppRvPx_mat[k]=gsl_matrix_alloc(n,GetNoExplan());
      ppPVPx_mat[k]=gsl_matrix_alloc(n,GetNoExplan());
      ppxPVPx_mat[k]=gsl_matrix_alloc(GetNoExplan(),GetNoExplan());
      ppxPVPy_vec[k]=gsl_vector_alloc(GetNoExplan());
    }
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, ppRvP_mat[k], px_mat, 0.0, ppRvPx_mat[k]); 
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, ppPVP_mat[k], px_mat, 0.0, ppPVPx_mat[k]);
  gsl_blas_dgemv(CblasTrans, 1.0, ppPVPx_mat[k], py_vec, 0.0, ppxPVPy_vec[k]); 
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, ppPVPx_mat[k], px_mat, 0.0, ppxPVPx_mat[k]); 
  
}

TwoVarCompModel::~TwoVarCompModel()
{
  ClearUp();
}

double TwoVarCompModel::Deviance(double sigma)
{
  // note 2nd par is k - we might want to try both ways.
  DevWorkspace dev_worker(this, false,0);
    
  return deviance(sigma, &dev_worker);
}


double TwoVarCompModel::NullDeviance(double sigma)
{
  DevWorkspace dev_worker(this, true, 0);
   
  return deviance(sigma, &dev_worker);
}

double TwoVarCompModel::MinimiseNullDeviance(double* pminimand)
{
  return MinimiseADeviance(true, pminimand);
}
double TwoVarCompModel::MinimiseDeviance(double *pminimand)
{ 
  return MinimiseADeviance(false, pminimand);
}


double TwoVarCompModel::MinimiseNullDeviance()
{ 
  return MinimiseADeviance(true, 0);
}
double TwoVarCompModel::MinimiseDeviance()
{
  return MinimiseADeviance(false, 0);
}




double TwoVarCompModel::MinimiseADeviance(bool null, double* pminimand)
{
  // this needs to be improved.
  int status;
  int iter = 0, max_iter = 1000;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m;//, m_expected = 0.5;
  // should really think about top limit
  double a = 0.0, b = 0.9999;
  gsl_function F;
  if(pminimand)
    m=*pminimand;
  else
    m=(a+b)/2.0;
  DevWorkspace dev_worker(this, null, 0);
  dev_pars Fpars;
  /* double deva= deviance(a, &dev_worker);
  double devm= deviance(m, &dev_worker);
  double devb= deviance(b, &dev_worker);*/
  // need a check for case that Kin=I or seem to get crash
  /* if((deva<=devm)||(devb<=devm))
    {
      if(deva<devb)
	{
	  // m=b;
    	  while((deviance(m, &dev_worker)>deva)&&(fabs(m-a)>0.00000001))
	    {
	      m=(m+a)/2.0;
	    }
	  if(!((fabs(m-a)>0.00000001)))
	    {
	      *pminimand=a;
	      return deva;
	    }
	}
      
      if(devb<deva)
	{
	 // m=a;
    	  while((deviance(m, &dev_worker)>devb)&&(fabs(m-b)>0.00000001))
	    {
	      m=(m+b)/2.0;
	    }
	  if(!((fabs(m-b)>0.00000001)))
	    {
	      *pminimand=b;
	      return devb;
	    }
	}
	}*/
  
  Fpars.pdev_work=&dev_worker;
  Fpars.pTVCM=this;
   
  F.function =  TwoVarCompModel_deviance_wrap;
  F.params = &Fpars;
  
  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, a, b);
  
  /*  printf ("using %s method\n",
	  gsl_min_fminimizer_name (s));
  
  printf ("%5s [%9s, %9s] %9s %10s %9s\n",
	  "iter", "lower", "upper", "min",
	  "err", "err(est)");
  
  printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
	  iter, a, b,
	  m, m - m_expected, b - a);
  */
  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);
     
      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);
    

      status = gsl_min_test_interval (a, b, 0.001, 0.0);
      
      /* if (status == GSL_SUCCESS)
	 printf ("Converged:\n");*/
     
      /* printf ("%5d [%.7f, %.7f] "
	      "%.7f %+.7f %.7f\n",
	      iter, a, b,
	      m, m - m_expected, b - a);*/
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  if(status!=GSL_SUCCESS)
    std::cout<<status;

  gsl_min_fminimizer_free (s);
  if(pminimand)
    *pminimand=m;
  return deviance(m, &dev_worker);
  
}
  
double TwoVarCompModel_deviance_wrap(double sigma, void* ppars)
{
  return ((dev_pars*) ppars)->pTVCM->deviance(sigma,((dev_pars*) ppars)->pdev_work);
}

double TwoVarCompModel::deviance(double sigma, DevWorkspace* pdev_work)
{

  // this assumes jeffreys prior at the moment
  // also assumes linear independence of covariates
  gsl_vector_memcpy(pdev_work->pa_vec, ppDminus1_vec[pdev_work->k]);
  //gslprint(ppDminus1_vec[pdev_work->k]);
  //gslprint(pdev_work->pa_vec);
 
   
  gsl_vector_scale(pdev_work->pa_vec, sigma);
  gsl_vector_add_constant(pdev_work->pa_vec, 1.0);
  
  gsl_vector_set_all(pdev_work->pa2_vec, sigma);
  gsl_vector_div(pdev_work->pa2_vec, pdev_work->pa_vec); 
  gsl_vector_memcpy(pdev_work->pa4_vec, pdev_work->pa2_vec);
  gsl_vector_mul(pdev_work->pa4_vec,ppRvPy_vec[pdev_work->k]);
  gsl_blas_ddot(pdev_work->pa4_vec, ppRvPy_vec[pdev_work->k], &(pdev_work->yPhiy));
 
  // this check to avoid 0/0 - best check?
  //if(fabs(pdev_work->yPhiy)>0)
    pdev_work->yPhiy*=-1.0;
    //else
    //  std::cout<<"yPhiy="<<pdev_work->yPhiy<<std::endl;

  pdev_work->yPhiy+=pyPVPy[pdev_work->k];

  // ADD ON THE DETERMINANT PART //
  pdev_work->deviance=0;
  for(pdev_work->it=0;pdev_work->it<n;pdev_work->it++)
    pdev_work->deviance+=log(gsl_vector_get(pdev_work->pa_vec,pdev_work->it));
  pdev_work->deviance-=rankCovar*log(1.0-sigma);
  
  if(pdev_work->null)
    {
      pdev_work->deviance+=(nu0+n-rankCovar)*log(pdev_work->yPhiy/(1-sigma)+(nu0*gamma0sq)); 
      return pdev_work->deviance;
    }
  
  gsl_matrix_memcpy(pdev_work->pa_mat, ppRvPx_mat[pdev_work->k]);
  for(pdev_work->it=0;pdev_work->it<n;pdev_work->it++)
    gsl_vector_scale(&(gsl_matrix_row(pdev_work->pa_mat,pdev_work->it).vector), gsl_vector_get(pdev_work->pa2_vec,pdev_work->it));
  
  gsl_vector_memcpy(pdev_work->pxPhiy_vec, ppxPVPy_vec[pdev_work->k]);
  gsl_matrix_memcpy(pdev_work->pxPhix_mat, ppxPVPx_mat[pdev_work->k]);

  gsl_blas_dgemv(CblasTrans, -1.0, pdev_work->pa_mat, ppRvPy_vec[pdev_work->k], 1.0, pdev_work->pxPhiy_vec);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, -1.0, pdev_work->pa_mat, ppRvPx_mat[pdev_work->k], 1.0, pdev_work->pxPhix_mat);

  // std::cout<<"xPhiy="<<std::endl;
  // gslprint(pdev_work->pxPhiy_vec);

   // _must be able to do faster than this! 
  //really ought to call a function to give hat{beta}
  // std::cout<<"xPhix=";
  // gslprint(pdev_work->pxPhix_mat);
  if(gsl_matrix_get(pdev_work->pxPhix_mat,0,0)>0.000000000000001)
    {
       moore_penrose(pdev_work->pxPhix_mat,CHOL_TOL);
      //  std::cout<<"xPhiy=";
      // gslprint(pdev_work->pxPhiy_vec);

      gsl_blas_dgemv(CblasTrans, 1.0, pdev_work->pxPhix_mat, pdev_work->pxPhiy_vec, 0.0, pdev_work->pa3_vec);
  
      gsl_blas_ddot(pdev_work->pxPhiy_vec,pdev_work->pa3_vec, &(pdev_work->SS));
      pdev_work->SS*=-1;
    }
  else
    { 
      std::cout<<"*";
      pdev_work->SS=0.0;
    }   

  // std::cout<<"yy/(1-sigma)="<<pdev_work->yPhiy/(1-sigma)<<std::endl;
  // std::cout<<"SS="<<pdev_work->SS<<std::endl;
  
  //  std::cout<<"SS"<<pdev_work->SS<<std::endl;
  // we add this for the null seperately, cant we do it better?
  pdev_work->SS+=pdev_work->yPhiy;
  //  std::cout<<"result/(1-sigma)="<<pdev_work->SS/(1-sigma)<<std::endl;
  //std::cout<<"SS+yy"<<pdev_work->SS<<std::endl;
  // std::cout<<nu0<<std::endl;
  //std::cout<<gamma0sq<<std::endl;
pdev_work->deviance+=(nu0+n-rankCovar)*log(pdev_work->SS/(1-sigma)+(nu0*gamma0sq));
  return (pdev_work->deviance);
}

unsigned int TwoVarCompModel::GetNoExplan()
{
  if(px_mat)
    return px_mat->size2;
  else
    return 0;
}

unsigned int TwoVarCompModel::GetSampleSize()
{
  return n;
}

void TwoVarCompModel::ClearUp()
{
   // free private members from the heap
  if(AmITheDaddy==1)
    {
      gsl_vector_free(py_vec);
      if(ppP_mat[0])
	gsl_matrix_free(ppP_mat[0]);
      if(ppP_mat[1])
	gsl_matrix_free(ppP_mat[1]);
      if(ppvar_mat[0])
	gsl_matrix_free(ppvar_mat[0]);
      if(ppvar_mat[1])
	gsl_matrix_free(ppvar_mat[1]);
      if(ppincid_mat[0])
	gsl_matrix_free(ppincid_mat[0]);
      if(ppincid_mat[1])
	gsl_matrix_free(ppincid_mat[1]);
      if(ppfullvar_mat[0])
	gsl_matrix_free(ppfullvar_mat[0]);
      if(ppfullvar_mat[1])
	gsl_matrix_free(ppfullvar_mat[1]);
      if(ppvR_mat[0])
	gsl_matrix_free(ppvR_mat[0]);
      if(ppvR_mat[1])
	gsl_matrix_free(ppvR_mat[1]);
      if(ppD_vec[0])
	gsl_vector_free(ppD_vec[0]);
      if(ppD_vec[1])
	gsl_vector_free(ppD_vec[1]);
      if(ppLambda_vec[0])
	gsl_vector_free(ppLambda_vec[0]);
      if(ppLambda_vec[1])
	gsl_vector_free(ppLambda_vec[1]);
      if(ppDminus1_vec[0])
	gsl_vector_free(ppDminus1_vec[0]);
      if(ppDminus1_vec[1])
	gsl_vector_free(ppDminus1_vec[1]);
      if(ppPVP_mat[0])
	gsl_matrix_free(ppPVP_mat[0]);
      if(ppPVP_mat[1])
	gsl_matrix_free(ppPVP_mat[1]);
      if(ppPVPy_vec[0])
	gsl_vector_free(ppPVPy_vec[0]);
      if(ppPVPy_vec[1])
	gsl_vector_free(ppPVPy_vec[1]);
      if(ppRvP_mat[0])
	gsl_matrix_free(ppRvP_mat[0]);
      if(ppRvP_mat[1])
	gsl_matrix_free(ppRvP_mat[1]);
      if(ppRvPy_vec[0])
	gsl_vector_free(ppRvPy_vec[0]);
      if(ppRvPy_vec[1])
      gsl_vector_free(ppRvPy_vec[1]);
    }
 
  if(ppRvPx_mat[0])
    gsl_matrix_free(ppRvPx_mat[0]);
  if(ppRvPx_mat[1])
    gsl_matrix_free(ppRvPx_mat[1]);
  if(ppPVPx_mat[0])
    gsl_matrix_free(ppPVPx_mat[0]);
  if(ppPVPx_mat[1])
    gsl_matrix_free(ppPVPx_mat[1]);
  if(ppxPVPx_mat[0])
    gsl_matrix_free(ppxPVPx_mat[0]);
  if(ppxPVPx_mat[1])
    gsl_matrix_free(ppxPVPx_mat[1]);
  if(ppxPVPy_vec[0])
    gsl_vector_free(ppxPVPy_vec[0]);
  if(ppxPVPy_vec[1])
    gsl_vector_free(ppxPVPy_vec[1]);
  
}

void TwoVarCompModel::SetScalarVarPrior(double nu0, double gamma0sq)
{
  this->gamma0sq=gamma0sq;
  this->nu0=nu0;
}


// best way of doing this?
void TwoVarCompModel::GetBeta(gsl_vector* pbeta_vec, double sigma)
{
  // note 2nd par is k - we might want to try both ways.
  DevWorkspace dev_worker(this, false,0);
    
  get_beta(pbeta_vec, sigma, &dev_worker);
  return;
}



void TwoVarCompModel::get_beta(gsl_vector* pbeta_vec, double sigma, DevWorkspace* pdev_work)
{
  // this assumes jeffreys prior at the moment
  // also assumes linear independence of covariates
  gsl_vector_memcpy(pdev_work->pa_vec, ppDminus1_vec[pdev_work->k]);

 
  gsl_vector_scale(pdev_work->pa_vec, sigma);
  gsl_vector_add_constant(pdev_work->pa_vec, 1.0);
  
  gsl_vector_set_all(pdev_work->pa2_vec, sigma);
   
  gsl_vector_div(pdev_work->pa2_vec, pdev_work->pa_vec); 
 
  gsl_vector_memcpy(pdev_work->pa4_vec, pdev_work->pa2_vec);

  gsl_vector_mul(pdev_work->pa4_vec,ppRvPy_vec[pdev_work->k]);
 
  

  gsl_matrix_memcpy(pdev_work->pa_mat, ppRvPx_mat[pdev_work->k]);

  for(pdev_work->it=0;pdev_work->it<n;pdev_work->it++)
    gsl_vector_scale(&(gsl_matrix_row(pdev_work->pa_mat,pdev_work->it).vector), gsl_vector_get(pdev_work->pa2_vec,pdev_work->it));

   
  gsl_vector_memcpy(pdev_work->pxPhiy_vec, ppxPVPy_vec[pdev_work->k]);
  gsl_matrix_memcpy(pdev_work->pxPhix_mat, ppxPVPx_mat[pdev_work->k]);

  gsl_blas_dgemv(CblasTrans, -1.0, pdev_work->pa_mat, ppRvPy_vec[pdev_work->k], 1.0, pdev_work->pxPhiy_vec);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, -1.0, pdev_work->pa_mat, ppRvPx_mat[pdev_work->k], 1.0, pdev_work->pxPhix_mat);

  // _must be able to do faster than this! 
  //really ought to call a function to give hat{beta}
  moore_penrose(pdev_work->pxPhix_mat,CHOL_TOL);
 

  gsl_blas_dgemv(CblasTrans, 1.0, pdev_work->pxPhix_mat, pdev_work->pxPhiy_vec, 0.0, pbeta_vec);
  
  return;
}
