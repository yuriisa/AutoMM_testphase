#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// need to understand const and make changes

typedef gsl_matrix* pgsl_matrix;
typedef gsl_vector* pgsl_vector;


double TwoVarCompModel_deviance_wrap(double sigma, void* ppars);

// just needed for debug
#include <iostream>

class TwoVarCompModel
{    
  friend double TwoVarCompModel_deviance_wrap(double sigma, void* ppars);
 public:
  class DevWorkspace{
 
  public:
    DevWorkspace(TwoVarCompModel* pparent, bool null, unsigned int k)
      {
	// all allocs are freed
	n=pparent->GetSampleSize();
	NoExplan=pparent->GetNoExplan();
	
	this->null=null;
	this->k=k;
	pa_vec=gsl_vector_alloc(n);
	pa2_vec=gsl_vector_alloc(n);		
	pa4_vec=gsl_vector_alloc(n);
	
	if((!null)&&(NoExplan!=0))
	  {

	    pa3_vec=gsl_vector_alloc(NoExplan);
	    pa_mat=gsl_matrix_alloc(n,NoExplan);
	
	    pxPhix_mat=gsl_matrix_alloc(NoExplan,NoExplan);
	    pxPhiy_vec=gsl_vector_alloc(NoExplan);
	  }
	
      }; 

    ~DevWorkspace()
      {
	gsl_vector_free(pa4_vec);
	gsl_vector_free(pa_vec);
	gsl_vector_free(pa2_vec);
	if((!null)&&(NoExplan!=0))
	  {
	    gsl_vector_free(pa3_vec);
	    gsl_matrix_free(pa_mat);
	    gsl_matrix_free(pxPhix_mat);
	    gsl_vector_free(pxPhiy_vec);
	  }
      };
   
    gsl_vector* pa_vec, *pa2_vec, *pa3_vec, *pa4_vec; 
    gsl_matrix* pa_mat;
    gsl_matrix* pxPhix_mat;
    gsl_vector* pxPhiy_vec;
    bool null;
    double SS;
    bool k;
    unsigned int n;
    unsigned int NoExplan;
    double deviance;
    double yPhiy;
    unsigned int it;
    
  }; private:
 bool AmITheDaddy;
 gsl_matrix_view x_matview;
 double MinimiseADeviance(bool, double*); 
  // should we use V rathar than fullvar?
  pgsl_matrix ppxPVPx_mat[2];
  pgsl_vector ppxPVPy_vec[2];
  pgsl_matrix ppRvPx_mat[2];
  pgsl_matrix ppPVPx_mat[2];

  pgsl_matrix ppPVP_mat[2];
  pgsl_vector ppPVPy_vec[2]; 
  double pyPVPy[2];
  pgsl_matrix ppRvP_mat[2];
  pgsl_vector ppRvPy_vec[2];

  pgsl_vector ppDminus1_vec[2];
 

  // var //
  // Eigenvalues of the covariance matrices
  pgsl_vector ppLambda_vec[2];
  // Basis for simultaneous diagonalisation
  pgsl_matrix ppvR_mat[2];
  // Projection 
  pgsl_matrix ppP_mat[2];
  // Covariates
  pgsl_matrix pw_mat;
  
  //Need to document what these do;
  pgsl_vector ppD_vec[2];

  double deviance(double sigmasq, DevWorkspace*);

  void get_beta(gsl_vector* pbeta_vec, double sigmasq, DevWorkspace*);

  pgsl_matrix ppfullvar_mat[2];
  
  // var that user might want to access 
  // sample size
  unsigned long int n;
  // reponse
  gsl_vector *py_vec;
  // explanatory variables
  gsl_matrix *px_mat;
  // variance components matrices
  pgsl_matrix ppvar_mat[2];
  bool var2isI;
  // variance components incidence matrices
  pgsl_matrix ppincid_mat[2];
  bool incidisI[2];
  // heritability parameters MLE
  double hMLE;
  // deviance
  //  double deviance;
  // regression coefficients
  gsl_vector *pregcoef_vec;

  // Clear up data members in preparation for object destruction.
  void ClearUp();


  // Setup the variance structure
  void SetVar(gsl_matrix* pvar1_mat, gsl_matrix* pincid1_mat, unsigned int k);
  
  void UpdatePostNewExplan(unsigned int k);
  // Compute the data members required to calculate likelihood
  void TakeYourMarks(unsigned int);  
  
  unsigned int rankCovar;
  double gamma0sq, nu0;
 public:
  // Constructor
  TwoVarCompModel(gsl_vector* presponse_vec, gsl_matrix* pcovar_mat, gsl_matrix* pvar1_mat, gsl_matrix* pvar2_mat, gsl_matrix* pincid1_mat, gsl_matrix* pincid2_mat, double* pnu_naught, double* pgammasq_naught);
  TwoVarCompModel(TwoVarCompModel &);
  void SetCovariates(gsl_matrix* pcovar_mat);
  void SetVar1(gsl_matrix* pvar1_mat, gsl_matrix* pincid1_mat);
  void SetVar2(gsl_matrix* pvar2_mat, gsl_matrix* pincid2_mat);
  void SetScalarVarPrior(double nu0, double gamma0sq);

  //Destructor
  ~TwoVarCompModel();
  // func //
  void FitModel(unsigned char //how to calculate likelihood: use dual rep or not
		);
  void SetExplan(gsl_matrix* pexplan_mat);
  void SetExplan(gsl_vector* pexplan_vec);
  
  void GetBeta(gsl_vector* pbeta_vec, double sigma);

  unsigned int GetSampleSize();
  unsigned int GetNoExplan();
  
  double MinimiseDeviance(); 
  double MinimiseNullDeviance();
  double MinimiseDeviance(double*); 
  double MinimiseNullDeviance(double*);
  double Deviance(double sigma);
  double NullDeviance(double sigma);
  
};

// Exceptions

// Constructor
#define EXC_CSTR_RESPNULL        1001            // py_vec argument of constructor is null
#define EXC_CSTR_NO_SAMPLE       1002            // Trying to crete object with no sample
#define EXC_CSTR_RESPALLOC       1003            // Allocation of this->py_vec failed in constructor

// ::SetCovariates
#define EXC_SETCOV_COVARNULL     2001            // pcovar_mat argument of SetCovariates is null;
#define EXC_SETCOV_COVARDIM      2002            // pcovar_mat argument of SetCovariates has incorrect first dimension
#define EXC_SETCOV_COVARALLOC    2003            // Allocation of this->pcovar_mat argument failed 

// ::SetVar
#define EXC_SETVAR_NULLPTR       3001            // pvar_mat or pincid_mat argument of SetVar is null
#define EXC_SETVAR_DIM           3002            // A vector or matrix argument of SetVar has incorrect dimension
#define EXC_SETVAR_ALLOC         3003            // Allocation of one of ppvar_mat[k], ppincid_mat[k], ppfullvar_mat[k] on heap failed 
#define EXC_SETVAR_AMATALLOC     3004            // Allocation of pa_mat on the heap failed
#define EXC_SETVAR_PSEUDOERROR   3005            // Failure to compute a moore-penrose inverse
#define EXC_SETVAR_PMATALLOC     3006            // Failure to allocate projection matrix ppP_mat[k] on the heap
#define EXC_SETVAR_KMATALLOC     3007            // Failure to allocate pK_mat on the heap.
#define EXC_SETVAR_KWMATALLOC    3008            // Failure to allocate pwK_mat or pwKw_mat on the heap.
