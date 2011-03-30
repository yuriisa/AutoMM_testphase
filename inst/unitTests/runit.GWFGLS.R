### --- Test setup ---

if(FALSE) {
  ## Not really needed, but can be handy when writing tests
  library("RUnit")
  library("GenABEL")
  library("DatABEL")
  library("MixABEL")
  library("MASS")
  library("mvtnorm")
}

### do not run
# stop("SKIP THIS TEST")
###

### ---- common functions and data -----

source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.GWFGLS_all_data_types <- function()
{
  
  unlink("tmp*")
  
  
  data(ge03d2.clean)
  
  idprop <- runif(1,min=0.2,max=0.5)
  IDS <- sort(sample(1:nids(ge03d2.clean),floor(nids(ge03d2.clean)*idprop)))
  dfd <- ge03d2.clean[IDS,autosomal(ge03d2.clean)]
  gkin <- ibs(dfd,w="freq")
  
  snpprop <- runif(1,min=0.005,max=0.01)
  SNPS <- sort(sample(1:nsnps(dfd),floor(nsnps(dfd)*snpprop)))
  dfd <- dfd[,SNPS]
  s <- summary(dfd@gtdata)
  maf <- pmin(s$Q.2,1-s$Q.2)
  dfd <- dfd[,(maf>0.05)]
  
  modelh2 <- 0.8
  covars <- c("sex","age")
  s2 <- 12^2
  betas_cov <- c(170,12,0.01)
  
  X <- as.matrix(cbind(rep(1,length(IDS)),phdata(dfd)[,covars]))
  grel <- gkin
  grel[upper.tri(grel)] <- t(grel)[upper.tri(grel)]
  grel <- 2*grel
  Y <- as.vector(rmvnorm(1,mean=(X %*% betas_cov),sigma=s2*(modelh2*grel+diag(dim(grel)[1])*(1-modelh2)))) 
  #+ 10*as.numeric(dfd@gtdata[,10])
  length(Y)
  Y[2] <- NA
  dfd@phdata$Y <- Y
  
  mdl <- Y~age+sex
  h2 <- polygenic(mdl,data=dfd,kin=gkin,gradtol=1e-5,steptol=1e-5,maxdiffgls=1e-2)
  h2$h2an
  mdl <- Y~age+sex+SNP
  
  # run tests with different types of genodata
  gtOld <- dfd@gtdata
  gtReal <- as.double(gtOld)
  gtNew <- as(gtReal,"databel")
  for (varcov in c(TRUE,FALSE))
    for (include.means in c(TRUE,FALSE))
      for (verbosity in c(0,1,2))
        for (residuals in c(FALSE,TRUE))
          for (test in c("wald","score","robust"))
            for (model in c("additive","dominantB","recessiveB","overdominant","genotypic")) {
              print(c(varcov,verbosity,residuals,test,model))
              aIR <- GWFGLS(mdl,data=phdata(dfd),genodata = gtReal,
                  model.SNP=model, test = test, residuals = residuals, 
                  verbosity = verbosity, varcov = varcov, 
                  include.means = include.means)
              aIO <- GWFGLS(mdl,data=phdata(dfd),genodata = gtOld, 
                  model.SNP=model, test = test, residuals = residuals, 
                  verbosity = verbosity, varcov = varcov, 
                  include.means = include.means)
              aIO1<- GWFGLS(mdl,data=dfd, 
                  model.SNP=model, test = test, residuals = residuals, 
                  verbosity = verbosity, varcov = varcov, 
                  include.means = include.means)
              aIN <- GWFGLS(mdl,data=phdata(dfd),genodata = gtNew, 
                  model.SNP=model, test = test, residuals = residuals, 
                  verbosity = verbosity, varcov = varcov, 
                  include.means = include.means)
              checkEquals(aIR,aIO,tolerance=5*.Machine$double.eps^0.5)
              checkEquals(aIR,aIO1,tolerance=5*.Machine$double.eps^0.5)
              checkEquals(aIR,aIN,tolerance=5*.Machine$double.eps^0.5)
              aWR <- GWFGLS(mdl,data=phdata(dfd),W=h2$InvSigma,genodata = gtReal,
                  model.SNP=model, test = test, residuals = residuals, 
                  verbosity = verbosity, varcov = varcov, 
                  include.means = include.means)
              aWO <- GWFGLS(mdl,data=phdata(dfd),W=h2$InvSigma,genodata = gtOld,
                  model.SNP=model, test = test, residuals = residuals, 
                  verbosity = verbosity, varcov = varcov, 
                  include.means = include.means)
              aWO1<- GWFGLS(mdl,data=dfd,W=h2$InvSigma,
                  model.SNP=model, test = test, residuals = residuals, 
                  verbosity = verbosity, varcov = varcov, 
                  include.means = include.means)
              aWN <- GWFGLS(mdl,data=phdata(dfd),W=h2$InvSigma,genodata = gtNew,
                  model.SNP=model, test = test, residuals = residuals, 
                  verbosity = verbosity, varcov = varcov, 
                  include.means = include.means)
              checkEquals(aWR,aWO,tolerance=5*.Machine$double.eps^0.5)
              checkEquals(aWR,aWO1,tolerance=5*.Machine$double.eps^0.5)
              checkEquals(aWR,aWN,tolerance=5*.Machine$double.eps^0.5)
            }
  
  
  rm(list=ls());gc()
  
  unlink("tmp*")
}


test.GWFGLS_equal_to_LM <- function()
{
  
  unlink("tmp*")
  
  
  data(ge03d2.clean)
  
  idprop <- runif(1,min=0.25,max=0.5)
  IDS <- sort(sample(1:nids(ge03d2.clean),floor(nids(ge03d2.clean)*idprop)))
  dfd <- ge03d2.clean[IDS,autosomal(ge03d2.clean)]
  gkin <- ibs(dfd,w="freq")
  
  snpprop <- runif(1,min=0.005,max=0.01)
  SNPS <- sort(sample(1:nsnps(dfd),floor(nsnps(dfd)*snpprop)))
  dfd <- dfd[,SNPS]
  smr <- summary(dfd@gtdata)
  maf <- pmin(smr$Q.2,1-smr$Q.2)
  dfd <- dfd[,(maf>0.01 & smr$P.11>2 & smr$P.12>2 & smr$P.22>2)]
  
  modelh2 <- 0.8
  covars <- c("sex","age")
  s2 <- 12^2
  betas_cov <- c(170,12,0.01)
  
  X <- as.matrix(cbind(rep(1,length(IDS)),phdata(dfd)[,covars]))
  grel <- gkin
  grel[upper.tri(grel)] <- t(grel)[upper.tri(grel)]
  grel <- 2*grel
  Y <- as.vector(rmvnorm(1,mean=(X %*% betas_cov),sigma=s2*(modelh2*grel+diag(dim(grel)[1])*(1-modelh2)))) 
  #+ 10*as.numeric(dfd@gtdata[,10])
  length(Y)
  Y[2] <- NA
  dfd@phdata$Y <- Y
  
  mdl <- Y~age+sex
  h2 <- polygenic(mdl,data=dfd,kin=gkin,gradtol=1e-5,steptol=1e-5,maxdiffgls=1e-2)
  h2$h2an
  mdl <- Y~age+sex+SNP
  
  # check equivalence to LM
  for (model in c("additive","dominantB","recessiveB","overdominant","genotypic")) {
    print(model)
    aIR <- GWFGLS(mdl,data=dfd,
        model.SNP=model, test = "wald", residuals = FALSE, 
        verbosity = 2, varcov = TRUE, with.lm = TRUE)
    passed <- rep(TRUE,dim(dfd@gtdata)[2])
    for (i in 1:dim(dfd@gtdata)[2]) {
      ggg <- get_snp_data(dfd@gtdata,i,model.SNP=model,gtcoding="typed")
      vr <- apply(ggg,FUN=var,MAR=2)
      if (any(vr<1e-8) || any(is.na(vr))) {passed[i] <- FALSE}
    }
    aIR$FGLS <- aIR$FGLS[passed,]
    aIR$LM <- aIR$LM[passed]
    checkEqualsNumeric(aIR$LM,aIR$FGLS[,"P-value"],tolerance=5*.Machine$double.eps^0.5)
  }
  
  rm(list=ls());gc()
  
  unlink("tmp*")
}

test.GWFGLS_Old_vs_New <- function()
{
  
#  library("RUnit")
#  library("GenABEL")
#  library("DatABEL")
#  library("MixABEL")
#  library("MASS")
#  library("mvtnorm")
  
  data(ge03d2.clean)
  propNids <- runif(1,min=0.25,max=0.5)
  IDS <- sort(sample(1:nids(ge03d2.clean),round(propNids*nids(ge03d2.clean))))
  NIDS <- length(IDS)
  df <- ge03d2.clean[IDS,autosomal(ge03d2.clean)]
  s <- summary(df@gtdata)
  maf <- pmin(s$Q.2,1-s$Q.2)
  df <- df[,(maf>0.05)]
  gkin <- ibs(df[,autosomal(df)],w="freq")
  
  h2mod <- 0.5
  covars <- c("sex","age")
  s2 <- 12^2
  betas_cov <- c(170,12,0.01)
  
  X <- as.matrix(cbind(rep(1,NIDS),phdata(ge03d2.clean)[1:NIDS,covars]))
  grel <- gkin
  grel[upper.tri(grel)] <- t(grel)[upper.tri(grel)]
  grel <- 2*grel
  Y <- as.vector(rmvnorm(1,mean=(X %*% betas_cov),sigma=s2*(h2mod*grel+diag(dim(grel)[1])*(1-h2mod))))
  Y[2] <- NA
  df@phdata$Y <- Y
  
  h2 <- polygenic(Y~sex+age,data=df,kin=gkin,gradtol=1e-5,steptol=1e-5,maxdiffgls=1e-2)
  h2$h2an
  
  propSnps <- runif(1,min=0.005,max=0.01)
  SNPS <- sort(sample(1:nsnps(df),round(propSnps*nsnps(df))))
  NSNPS <- length(SNPS)
  
  dfn <- df[,SNPS]
  smr <- summary(dfn@gtdata)
  ## actually should automatically do that!
  maf <- pmin(smr$Q.2,1.0-smr$Q.2)
  #daf <- pmin(smr$P.11/smr$NoMeasured,1.0-smr$P.11/smr$NoMeasured)
  #raf <- pmin(smr$P.22/smr$NoMeasured,1.0-smr$P.22/smr$NoMeasured)
  #oaf <- pmin(smr$P.12/smr$NoMeasured,1.0-smr$P.12/smr$NoMeasured)
  dfn <- dfn[,which(maf>0.01 & smr$P.11 > 7 & smr$P.12 > 7 & smr$P.22 > 7)]
  #dfn <- dfn[,which(maf<0.01 | smr$P.11 < 3 | smr$P.12 < 3 & smr$P.22 < 3)]
  nsnps(dfn)
  #summary(dfn@gtdata)
  
  # check that iterator produces NAs with mono-snps
  resIteratorNew <- GWFGLS(Y~age+sex*SNP,data=dfn@phdata,W=h2$InvSigma,old=FALSE,ver=2,varcov=TRUE,
      include.means=FALSE,genodata=matrix(rep(1,length(Y)*2),ncol=2))
  print(resIteratorNew)
  checkTrue(all(is.na(resIteratorNew[1,4:dim(resIteratorNew)[2]])))
  checkTrue(all(is.na(resIteratorNew[2,4:dim(resIteratorNew)[2]])))
  
  t0 <- proc.time()
  resIteratorOld <- GWFGLS(Y~age+sex*SNP,data=dfn,W=h2$InvSigma,old=TRUE)
  proc.time() - t0
  t0 <- proc.time()
  resIteratorNew <- GWFGLS(Y~age+sex*SNP,data=dfn,W=h2$InvSigma,old=FALSE)
  proc.time() - t0
  newNA <- (resIteratorNew[,"Chisq"]<0)
  print(table(newNA))
  resIteratorNew <- resIteratorNew[!newNA,]
  resIteratorOld <- resIteratorOld[!newNA,]
  checkEquals(resIteratorOld,resIteratorNew)
  
  
  for (varcov in c(TRUE,FALSE))
    for (include.means in c(TRUE,FALSE))
      for (verbosity in c(1,2))
        for (residuals in c(FALSE,TRUE))
          for (test in c("wald","score","robust"))
            for (model in c("additive","dominantB","recessiveB","overdominant","genotypic")) {
              print(c(varcov,verbosity,residuals,test,model))
              #for (i in 1:dim(dfn@gtdata)[2]) {
              resIteratorOld <- GWFGLS(Y~age+sex*SNP,data=dfn,W=h2$InvSigma,old=TRUE,
                  verbosity=verbosity,varcov=varcov,include.means=include.means,
                  test=test,residuals=residuals,model=model,
                  with.lm=TRUE)
              resIteratorNew <- GWFGLS(Y~age+sex*SNP,data=dfn,W=h2$InvSigma,old=FALSE,
                  verbosity=verbosity,varcov=varcov,include.means=include.means,
                  test=test,residuals=residuals,model=model)
              newNA <- (resIteratorNew[,"Chisq"]<0)
              print(table(newNA))
              #print(resIteratorOld)
              #print(resIteratorNew)
              #print(summary(lm(Y~age+sex+as.numeric(dfn@gtdata[,i]),data=dfn@phdata)))
              #if (any(newNA)) checkEqualsNumeric(0,mean(resIteratorOld[newNA,"Chisq"]))
              resIteratorNew <- resIteratorNew[!newNA,,drop=FALSE]
              resIteratorOld <- resIteratorOld$FGLS[!newNA,,drop=FALSE]
              checkEquals(resIteratorOld,resIteratorNew,tolerance=5*(.Machine$double.eps^0.5))
              #}
            }
  
}
