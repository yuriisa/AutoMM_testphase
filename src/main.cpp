#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "twovarcomp.h"

#ifdef _OPENMP
#include <omp.h>
#define OMP_GET_MAX_THREADS omp_get_max_threads()
#define OMP_GET_THREAD_NUM omp_get_thread_num()
#else
#define OMP_GET_MAX_THREADS 1
#define OMP_GET_THREAD_NUM 0
#endif // _OPENMP

// for debug
#include "moorepenrose.h"
extern "C" {

#include "main.h"


SEXP rint_flmm(SEXP pexplan_sexp, SEXP presp_sexp, SEXP pn_sexp, SEXP pp_sexp, SEXP pcovar_sexp, SEXP pp_covar_sexp, SEXP pVar2_sexp, SEXP nu_naught_sexp, SEXP gamma_naught_sexp)
{

  double *pexplan, *presp,  *pnu_naught, *pgamma_naught, *pcovar, *pVar2;
  double* pchisq;
  double* pherit;
  unsigned int *pn, *pp_covar, *pp;
 

  char pret_names[][100]={"coefs", "chi.sq", "herit", "null.herit"};

  SEXP preturn_list_SEXP, preturn_names_SEXP, paname_SEXP;
    
  SEXP pchisq_SEXP;
  SEXP pherit_SEXP;
  SEXP pbeta_SEXP;
  SEXP pnullherit_SEXP;

  gsl_matrix* pvar1_mat, *pvar2_mat;
  gsl_matrix* pcovar_mat;
 
  gsl_vector* presponse_vec;

  double* pnullherit;
  double* pbeta;
  

  // really must check all gsl returns
  gsl_set_error_handler_off();

  // C side pointers to R objects
  pexplan=(double*) REAL(pexplan_sexp);
  presp=(double*) REAL(presp_sexp);
  pn=(unsigned int*) INTEGER(pn_sexp);
  pp=(unsigned int*) INTEGER(pp_sexp);
  pcovar=(double*) REAL(pcovar_sexp);
  pp_covar=(unsigned int*) INTEGER(pp_covar_sexp);
  pVar2=(double*) REAL(pVar2_sexp);
  pnu_naught=(double*) REAL(nu_naught_sexp);
  pgamma_naught=(double*) REAL(gamma_naught_sexp);
  
  
  /* gsl_vector_view response_vecview=gsl_vector_view_array(presp, *pn);
     presponse_vec=&(response_vecview.vector);*/
  gsl_matrix_view var2_matview=gsl_matrix_view_array(pVar2, *pn, *pn);
  pvar2_mat=&(var2_matview.matrix);
  // freed
  pvar1_mat=gsl_matrix_alloc(*pn, *pn);
  /* gsl_matrix_view covar_matview=gsl_matrix_view_array(pcovar, *pn, *pp_covar);
     pcovar_mat=&(covar_matview.matrix);*/

  // sort out missing in response
  // better to bulk copy then iterate?
  unsigned int it;
  unsigned int nonzerocount=0;
  // freed
  presponse_vec=gsl_vector_alloc(*pn);
  double meanval=0.0;
  for(it=0;it<*pn;it++)
    {
      if(!ISNA(presp[it]))
	{
	  meanval+=presp[it];
	  nonzerocount+=1;
	}
    }
  meanval/=(double) nonzerocount;
  
  for(it=0;it<*pn;it++)
    {
      if(ISNA(presp[it]))
	gsl_vector_set(presponse_vec, it, meanval);
      else
	gsl_vector_set(presponse_vec, it, presp[it]);  
    }

  // freed
  pcovar_mat=gsl_matrix_alloc( *pn, *pp_covar);
  unsigned it2;
  for(it2=0;it2<*pp_covar;it2++)
    {
      meanval=0.0;
      nonzerocount=0;
      for(it=0;it<*pn;it++)
	{
	  if(!ISNA(pcovar[it*(*pp_covar)+it2]))
	    {
	      meanval+=pcovar[it*(*pp_covar)+it2];
	      nonzerocount+=1;
	    }
	}
      meanval/=(double) nonzerocount;
      for(it=0;it<*pn;it++)
	{
	  if(ISNA(pcovar[it*(*pp_covar)+it2]))
	    gsl_matrix_set(pcovar_mat, it, it2, meanval);
	  else
	    gsl_matrix_set(pcovar_mat, it, it2, pcovar[it*(*pp_covar)+it2]);  
	}
    }

  
  /*std::cout<<"cov="<<pcovar[0]<<","<<pcovar[1]<<","<<pcovar[2]<<std::endl;
    std::cout<<"pcovar_mat";
  gslprint(pcovar_mat);*/
  

  gsl_matrix* pincid1_mat, *pincid2_mat;
 
  // freed
  pincid1_mat=gsl_matrix_alloc(*pn,*pn);
  //freed
  pincid2_mat=gsl_matrix_alloc(*pn,*pn);

  gsl_matrix_set_identity(pvar1_mat);
  // gsl_matrix_set_identity(pvar2_mat); 
  gsl_matrix_set_identity(pincid1_mat);
  gsl_matrix_set_identity(pincid2_mat); 

  PROTECT(pbeta_SEXP=NEW_NUMERIC(*pp));
  pbeta=NUMERIC_POINTER(pbeta_SEXP);
  PROTECT(pchisq_SEXP=NEW_NUMERIC(*pp));
  pchisq=NUMERIC_POINTER(pchisq_SEXP);
  PROTECT(pherit_SEXP=NEW_NUMERIC(*pp));
  pherit=NUMERIC_POINTER(pherit_SEXP);
  PROTECT(pnullherit_SEXP=NEW_NUMERIC(1));
  pnullherit=NUMERIC_POINTER(pnullherit_SEXP);

  
  TwoVarCompModel DaddyTwoVarCompModel(presponse_vec, pcovar_mat, pvar1_mat, pvar2_mat, pincid1_mat, pincid2_mat, pnu_naught, pgamma_naught);  
  double nullminimand=0.5;
  double altminimand;
  double nulldev=DaddyTwoVarCompModel.MinimiseNullDeviance(&nullminimand);
  *pnullherit=nullminimand;
  
  /*std::cout<<"si==0.2"<<std::endl<<DaddyTwoVarCompModel.NullDeviance(0.2)<<std::endl;	 
  std::cout<<"si==0.4"<<std::endl<<DaddyTwoVarCompModel.NullDeviance(0.4)<<std::endl;	 
  std::cout<<"si==0.6"<<std::endl<<DaddyTwoVarCompModel.NullDeviance(0.6)<<std::endl;
  std::cout<<"si==0.8"<<std::endl<<DaddyTwoVarCompModel.NullDeviance(0.8)<<std::endl;
  */
  
 
  pgsl_vector* ppexplantemp_vec = new pgsl_vector[OMP_GET_MAX_THREADS];
  pgsl_vector* ppbeta_vec=new pgsl_vector[OMP_GET_MAX_THREADS];
  for(it=0;it<OMP_GET_MAX_THREADS;it++)
    {
      ppexplantemp_vec[it]=gsl_vector_alloc(*pn);
      ppbeta_vec[it]=gsl_vector_alloc(1);

    }
  #pragma omp parallel for shared(pexplan, pp, pn, pchisq, pherit, nulldev, nullminimand, ppexplantemp_vec, pbeta, ppbeta_vec) private(altminimand, it2, meanval, nonzerocount)
  for(it=0;it<*pp;it++)
    {
      TwoVarCompModel ChildTwoVarCompModel(DaddyTwoVarCompModel);
      std::cout<<".";
      
      meanval=0.0;
      nonzerocount=0;
      for(it2=0;it2<*pn;it2++)
	{
	  if(!ISNA(pexplan[it2+(*pn)*it]))
		{
		  meanval+=pexplan[it2+(*pn)*it];
		  nonzerocount+=1;
		}
	}
      meanval/=(double) nonzerocount;
      
      for(it2=0;it2<*pn;it2++)
	{
	  if(ISNA(pexplan[it2+(*pn)*it]))
	    gsl_vector_set(ppexplantemp_vec[OMP_GET_THREAD_NUM], it2, meanval);
	  else
	    gsl_vector_set(ppexplantemp_vec[OMP_GET_THREAD_NUM], it2, pexplan[it2+(*pn)*it]);  
	}
      ChildTwoVarCompModel.SetExplan(ppexplantemp_vec[OMP_GET_THREAD_NUM]);
      altminimand=nullminimand;
    
      pchisq[it]=nulldev-ChildTwoVarCompModel.MinimiseDeviance(&altminimand);
      
      pherit[it]=altminimand;
      ChildTwoVarCompModel.GetBeta(ppbeta_vec[OMP_GET_THREAD_NUM], altminimand);
      pbeta[it]=gsl_vector_get(ppbeta_vec[OMP_GET_THREAD_NUM], 0);
      
    }
  for(it=0;it<OMP_GET_MAX_THREADS;it++)
    {
      gsl_vector_free(ppexplantemp_vec[it]);
      gsl_vector_free(ppbeta_vec[it]);
    }
  delete[] ppexplantemp_vec;
  delete[] ppbeta_vec;
  
  gsl_matrix_free(pvar1_mat);
  gsl_vector_free(presponse_vec);
  gsl_matrix_free(pcovar_mat);
  gsl_matrix_free(pincid1_mat);
  gsl_matrix_free(pincid2_mat);
  
  
  PROTECT(preturn_list_SEXP=allocVector(VECSXP,4));
  SET_VECTOR_ELT(preturn_list_SEXP, 0,pbeta_SEXP);
  SET_VECTOR_ELT(preturn_list_SEXP, 1,pchisq_SEXP);
  SET_VECTOR_ELT(preturn_list_SEXP, 2,pherit_SEXP);
  SET_VECTOR_ELT(preturn_list_SEXP, 3,pnullherit_SEXP);
 
 
  PROTECT(preturn_names_SEXP=allocVector(STRSXP,4));

  
  for(int it=0;it<4;it++)
    {
      
      PROTECT(paname_SEXP=Rf_mkChar(pret_names[it]));
      SET_STRING_ELT(preturn_names_SEXP,it,paname_SEXP);
    }
  setAttrib(preturn_list_SEXP, R_NamesSymbol,preturn_names_SEXP);
  
  UNPROTECT(10);

  return preturn_list_SEXP;
}

}

