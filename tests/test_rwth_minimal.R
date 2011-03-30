library(GenABEL)
library(MixABEL)
data(ge03d2.clean)
df <- ge03d2.clean[,autosomal(ge03d2.clean)]
#set.seed(1)
#gkin <- ibs(df[,sort(sample(1:nsnps(df),500))],w="freq")
#h2ht <- polygenic_hglm(height~sex+age,data=df,kin=gkin,conv=1e-6, maxit=50)
#save(gkin,h2ht,file="test_rwth_minimal.RData")
load("test_rwth_minimal.RData")

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
time0 <- system.time(
    fgls_full_res <- GWFGLS(y~SNP,data=df[,1],
        W=h2ht$InvSigma,verbosity=2,varcov=TRUE)
)
time0
time_full <- system.time(
		fgls_full_res <- GWFGLS(y~SNP,data=df,
				W=h2ht$InvSigma,verbosity=2,varcov=TRUE)
)
time_full-time0
#
# here beta_(Int.. and beta_SNP are elemnts of beta, and 
# se... make diagonal of 'e' and cov_.. is an off-diag elem 
# (so 'e' is 2x2 matrix but symmetric, so 3 numbers)
time_min <- system.time(
 fgls_min_res <- GWFGLS(y~SNP,data=df,
  W=h2ht$InvSigma,verbosity=2,varcov=TRUE,minimal=TRUE)
)
time_min-time0
fgls_full_res[1:3,]		
fgls_min_res[1:3,]