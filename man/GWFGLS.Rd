\name{GWFGLS}
\alias{GWFGLS}
\title{Genome-wide FGLS...}
\usage{GWFGLS(formula, data, subset, weights, na.action, contrasts, offset, W,
    inverse=TRUE, na.SNP="impute", mincall=0.95, residuals=FALSE,
    test="wald", model.SNP="additive", genodata, gtcoding="typed",
    verbosity=1, varcov=FALSE, include.means=TRUE, singular="ignore",
    with.lm=FALSE, old=FALSE)}
\description{Genome-wide FGLS}
\details{Genome-wide Feasible GLS}
\author{Yurii Aulchenko}
\arguments{\item{formula}{analysis formula; should contain special 'SNP' term}
\item{data}{phenotypic data frame, or 'gwaa.data-class' object}
\item{subset}{subset of data (individuals)}
\item{weights}{RESERVED FOR FURTURE USE}
\item{na.action}{RESERVED FOR FURTURE USE}
\item{contrasts}{RESERVED FOR FURTURE USE}
\item{offset}{RESERVED FOR FURTURE USE}
\item{W}{for GLS, (inverse of) variance-covariance matrix, as such returned by 
GenABEL's polygenic(...)\$InvSigma, or NULL for LS}
\item{inverse}{wheather W is already inverted}
\item{na.SNP}{how to deal with missing SNP data; 'impute' -- substitute 
with mean, 'drop' -- drop rows with missing genotypes}
\item{mincall}{minimall call rate for a SNP (if less, the SNP is dropped)}
\item{residuals}{use residuals for analysis? (faster, but less precise)}
\item{test}{test to be applied, one of 'wald', 'score' or 'robust'}
\item{model.SNP}{SNP model to apply, one of 
c("additive","dominantB","recessiveB","overdominant","genotypic")}
\item{genodata}{genotypic data; can be missing if data is of 'gwaa.data-class'.
Otherwise can be regular matrix or 'databel' matrix}
\item{gtcoding}{one of c("typed","dose","probability")
'typed' -- coded with NA, 0, 1, 2}
\item{verbosity}{what to report; 0 -- only test stats, 
1 -- test stats and estimates concerning 'whichtest' parameters, 
2 -- test stats and estimates of all parameters}
\item{varcov}{wheather var-cov matrix for estimated parameters 
is to be reported}
\item{include.means}{weather mean values of predictors are to 
be reported}
\item{singular}{what to do with linear dependencies in X (now only 
'ignore' implemented)}
\item{with.lm}{whether LM should be run along (only test purposes; always 
falls back to 'old' R-only implementation)}
\item{old}{if TRUE, old R-only code implementation is running (testing purposes, 
slow)}}
\examples{library(MASS)
library(mvtnorm)
library(GenABEL)
data(ge03d2.clean)
NIDS <- 100
df <- ge03d2.clean[1:NIDS,autosomal(ge03d2.clean)]
s <- summary(df@gtdata)
maf <- pmin(s$Q.2,1-s$Q.2)
df <- df[,(maf>0.05)]
gkin <- ibs(df[,autosomal(df)],w="freq")

modelh2 <- 0.8
covars <- c("sex","age")
s2 <- 12^2
betas_cov <- c(170,12,0.01)

X <- as.matrix(cbind(rep(1,NIDS),phdata(ge03d2.clean)[1:NIDS,covars]))
grel <- gkin
grel[upper.tri(grel)] <- t(grel)[upper.tri(grel)]
grel <- 2*grel
Y <- as.vector(rmvnorm(1,mean=(X \%*\% betas_cov),sigma=s2*(modelh2*grel+diag(dim(grel)[1])*(1-modelh2))))
length(Y)
Y[2] <- NA
df@phdata$Y <- Y

mdl <- Y~age+sex
h2 <- polygenic(mdl,data=df,kin=gkin)
h2$h2an
mdl <- Y~age+sex+SNP

gtOld <- df@gtdata[,1:10]
gtReal <- as.double(gtOld)
gtNew <- as(gtReal,"databel")
aIR <- GWFGLS(mdl,data=phdata(df),genodata = gtReal,verbosity=0) #,model="genotypic")
aIO <- GWFGLS(mdl,data=phdata(df),genodata = gtOld, verbosity=0) #,model="genotypic")
aIN <- GWFGLS(mdl,data=phdata(df),genodata = gtNew, verbosity=0) #,model="genotypic")
aWR <- GWFGLS(mdl,data=phdata(df),W=h2$InvSigma,genodata = gtReal,verbosity=0) #,model="genotypic")
aWO <- GWFGLS(mdl,data=phdata(df),W=h2$InvSigma,genodata = gtOld, verbosity=0) #,model="genotypic")
aWN <- GWFGLS(mdl,data=phdata(df),W=h2$InvSigma,genodata = gtNew, verbosity=0) #,model="genotypic")
aIN
aWN
complex_model <- GWFGLS(Y~age+sex*SNP,data=phdata(df),W=h2$InvSigma,
genodata = gtReal,verbosity=2,varcov=TRUE,include.means=TRUE,model="genotypic")
complex_model}
