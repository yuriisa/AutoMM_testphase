  library("RUnit")
  library("GenABEL")
  library("DatABEL")
  library("MixABEL")
  library("MASS")
  library("mvtnorm")
  
  data(ge03d2.clean)
  propNids <- 0.2 #runif(1,min=0.5,max=0.75)
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
  
  h2 <- polygenic(Y~sex+age,data=df,kin=gkin,gradtol=1e-6,steptol=1e-6,maxdiff=1e-2)
  h2$h2an
  
  propSnps <- runif(1,min=0.005,max=0.01)
  SNPS <- sort(sample(1:nsnps(df),round(propSnps*nsnps(df))))
  NSNPS <- length(SNPS)
  
  dfn <- df[,SNPS]
  smr <- summary(dfn@gtdata)
  ## actually should automatically do that!
  maf <- pmin(smr$Q.2,1.0-smr$Q.2)
  daf <- pmin(smr$P.11/smr$NoMeasured,1.0-smr$P.11/smr$NoMeasured)
  raf <- pmin(smr$P.22/smr$NoMeasured,1.0-smr$P.22/smr$NoMeasured)
  oaf <- pmin(smr$P.12/smr$NoMeasured,1.0-smr$P.12/smr$NoMeasured)
  dfn <- dfn[,which(maf>0.01 & smr$P.11>5 & smr$P.12 > 5 & smr$P.22 > 5)]
  #dfn <- dfn[,which(maf<0.01 | daf<0.01 | raf < 0.01 | oaf < 0.01)]
  
  # check that iterator produces NAs with mono-snps
  resIteratorNew <- GWFGLS(Y~age+sex*SNP,data=dfn@phdata,W=h2$InvSigma,old=FALSE,ver=2,varcov=TRUE,
      include.means=FALSE,genodata=matrix(rep(1,length(Y)*2),ncol=2))
  
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
      for (verbosity in c(0,1,2))
        for (residuals in c(FALSE,TRUE))
          for (test in c("wald","score","robust"))
            for (model in c("additive","dominantB","recessiveB","overdominant","genotypic")) {
              print(c(varcov,verbosity,residuals,test,model))
              #for (i in 1:dim(dfn@gtdata)[2]) {
              resIteratorOld <- GWFGLS(Y~age+sex*SNP,data=dfn,W=h2$InvSigma,old=TRUE,
                  verbosity=verbosity,varcov=varcov,include.means=FALSE,
                  test=test,residuals=residuals,model=model,
                  with.lm=TRUE)
              resIteratorNew <- GWFGLS(Y~age+sex*SNP,data=dfn,W=h2$InvSigma,old=FALSE,
                  verbosity=verbosity,varcov=varcov,include.means=FALSE,
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
  

