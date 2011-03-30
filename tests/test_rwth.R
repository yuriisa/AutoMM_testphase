library(GenABEL)
library(MixABEL)
data(ge03d2.clean)
df <- ge03d2.clean[,autosomal(ge03d2.clean)]

load("saved_gkin_and_h2.RData")
#gkin <- ibs(df,w="freq")

# toy example to see if the data transfer worked
NIDS <- 7
NSNPS <- 5
phi <- matrix(c(1:(NIDS^2)),ncol=NIDS)
phi
XRs <- as.numeric(gtdata(df[1:NIDS,1:NSNPS]))
XRs[is.na(XRs)] <- 0.5
XRs
y <- phdata(df)$height[1:NIDS]
y[is.na(y)] <- mean(y,na.rm=T)
y
XL <- as.matrix(phdata(df)[1:NIDS,c("sex","age")])
XL
if (any(is.na(XL))) stop("opan'ki!")

colMeans(XL)
colMeans(XRs)
a <- rwth_caller(1,phi,0.3,12.1,y,XL,XRs,1)
a

# repeat toy with different XR
tmp <- matrix(NA,ncol=2*dim(XRs)[2],nrow=dim(XRs)[1])
for (i in 1:dim(XRs)[2]) {
	tmp[,i*2-1] <- 1*(XRs[,i]==1)
	tmp[,i*2] <- 1*(XRs[,i]==2)
}
XRs <- tmp
colMeans(XRs)
a <- rwth_caller(1,phi,0.3,12.1,y,XL,XRs,2)
a


# realistic data
# note that h2 is actually on lower end for this example
# (where ~ 0.11, more frequently you see 0.25-0.75)
# also phi is on 'sparse' end -- we may have even more sparse
# data, but our typical phi is somewaht denser 
# (see 2,400 matrix Paolo has)
# estimate model 
# h2ht <- polygenic(height~sex+age,data=df,kin=gkin)
# save(gkin,h2ht,file="saved_gkin_and_h2.RData")
# h2ht$h2an
phi <- gkin[h2ht$measuredIDs,h2ht$measuredIDs]
phi[upper.tri(phi)] <- t(phi)[upper.tri(phi)]
phi[1:5,1:5]
y <- h2ht$residualY[h2ht$measuredIDs]
any(is.na(y))
h2 <- h2ht$est
h2
Sigma2 <- var(y)
Sigma2
XL <- NULL
XRs <- as.numeric(gtdata(df[h2ht$measuredIDs,1:10]))
cmn <- colMeans(XRs,na.rm=T)
for (i in 1:ncol(XRs)) {
	XRs[is.na(XRs[,i]),i] <- cmn[i]
}
#XRs
a <- rwth_caller(0,phi,h2,Sigma2,y,XL,XRs,1)
cmn
a

old_res <- GWFGLS(y~SNP,data=phdata(df[h2ht$measuredIDs,]),
		W=h2ht$InvSigma,genodata = XRs,verbosity=2,varcov=TRUE)
# here beta_(Int.. and beta_SNP are elemnts of beta, and 
# se... make diagonal of 'e' and cov_.. is an off-diag elem 
# (so 'e' is 2x2 matrix but symmetric, so 3 numbers)
old_res		
