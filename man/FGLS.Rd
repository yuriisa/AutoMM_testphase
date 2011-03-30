\name{FGLS}
\alias{FGLS}
\title{Feasible GLS...}
\usage{FGLS(Y, X, W, test="wald", whichtest=c(FALSE, rep(TRUE, dim(X)[2] -
    1)))}
\description{Feasible GLS}
\details{Feasible Generalised Least Squares}
\value{List with elements 'beta' -- estimates fo the regression 
coefficients; 'V' -- variance covariance matrix for parameters 
estimates; 'T2' -- test statistics (distributed as Chi-squared under 
the null) for the testing of whichtest parameters; 'df' -- the number of 
degrees of freedom of the T2 test; 'tested' -- which parameters were 
tested with T2; 'meanColX' -- mean value of variable in columns of X;
'n' -- length of Y (== height of X)}
\author{Yurii Aulchenko}
\arguments{\item{Y}{dependent variable}
\item{X}{design matrix (including intercept, if necessary)}
\item{test}{test to be applied, one of 'wald', 'score' or 'robust'}
\item{whichtest}{which independent variables to be tested (set to 'TRUE')}
\item{W}{for GLS, inverse variance-covariance matrix, as such returned by 
GenABEL's polygenic(...)$InvSigma, or NULL for LS}}
