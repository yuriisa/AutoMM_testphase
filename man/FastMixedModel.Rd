\name{FastMixedModel}
\alias{FastMixedModel}
\title{fast mixed models...}
\usage{FastMixedModel(Response, Explan, Kin, Covariates, nu_naught=0,
    gamma_naught=0)}
\description{fast mixed models}
\details{fast mixed models -- BETA VERSION
If compiled against OMP this library can exploit multi-core parallelism 
Does not cope with missing data at present}
\value{a list with values ...}
\seealso{mmscore}
\references{reference to fill in}
\author{William Astle \email{fio@where}}
\arguments{\item{Response}{is an n dimensional vector of Responses}
\item{Explan}{is an n*p matrix of Explanatory variables, each to be tested marginally (SNPS)}
\item{Kin}{is the n*n kinship matrix}
\item{Covariates}{is an n*k matrix of Covariates}
\item{nu_naught}{(and gamma_naught) are hyperparameters which control the heaviness 
of the tails of the test distribution (recommend leave them unchanged).}
\item{gamma_naught}{(and nu_naught) are hyperparameters which control the heaviness 
of the tails of the test distribution (recommend leave them unchanged).}}
\examples{require(mvtnorm)
data(ge03d2.clean)
df <- ge03d2.clean[1:250,autosomal(ge03d2.clean)]
NSNPS <- nsnps(df)
modh2 <- 0.8
gkin <- ibs(df[,autosomal(df)],w="freq")

ngkin <- gkin
ngkin[upper.tri(ngkin)] <- t(ngkin)[upper.tri(ngkin)] 
ngkin[1:5,1:5]
mysig <- (modh2*2*ngkin+(1.-modh2)*diag(dim(ngkin)[1]))
mysig[1:5,1:5]
mytra <- as.vector(rmvnorm(1,sigma=mysig)) + phdata(df)$sex*0.05 + phdata(df)$age*0.002
mytra[1:10]
df@phdata$mytra <- mytra
df@phdata[1:5,]

time0.h2 <- proc.time()
h2 <- polygenic(mytra~sex+age,data=df,kin=gkin)
time.h2 <- proc.time() - time0.h2

time0.mms <- proc.time()
mms <- mmscore(h2,data=df)
time.mms <- proc.time() - time0.mms

time0.grs <- proc.time()
grs <- qtscore(h2$pgres,data=df)
time.grs <- proc.time() - time0.grs

res <- mytra
summary(res)
expl <- as.numeric(df[,1:NSNPS])
summary(res)
covariates <- matrix(c(phdata(df)$sex,phdata(df)$age),ncol=2)
summary(covariates) 

time0.fmm <- proc.time()
fmm <- FastMixedModel(Response=res,
Explan=expl,
Kin = gkin,
Cov=covariates)
time.fmm <- proc.time() - time0.fmm

time.h2
time.h2+time.grs
time.h2+time.mms
time.fmm

h2$h2an
#mms$effB
#mms$chi2.1df
fmm$null.herit

cor(mms[,"chi2.1df"],fmm$chi.sq)^2
plot(mms[,"chi2.1df"],fmm$chi.sq)}
