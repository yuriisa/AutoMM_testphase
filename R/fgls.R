#' Feasible GLS
#' 
#' Feasible Generalised Least Squares
#' 
#' @param Y dependent variable
#' @param X design matrix (including intercept, if necessary)
#' @param test test to be applied, one of 'wald', 'score' or 'robust'
#' @param whichtest which independent variables to be tested (set to 'TRUE')
#' @param W for GLS, inverse variance-covariance matrix, as such returned by 
#' GenABEL's polygenic(...)$InvSigma, or NULL for LS
#' 
#' @return List with elements 'beta' -- estimates fo the regression 
#' coefficients; 'V' -- variance covariance matrix for parameters 
#' estimates; 'T2' -- test statistics (distributed as Chi-squared under 
#' the null) for the testing of whichtest parameters; 'df' -- the number of 
#' degrees of freedom of the T2 test; 'tested' -- which parameters were 
#' tested with T2; 'meanColX' -- mean value of variable in columns of X;
#' 'n' -- length of Y (== height of X)
#'
#' @author Yurii Aulchenko
#'
FGLS <- function(Y,X,W=NULL,test="wald",whichtest = c(FALSE,rep(TRUE,dim(X)[2]-1))) 
{
	# do checks
	if (!is(whichtest,"logical") || length(whichtest)!=dim(X)[2]) 
		stop("argument whichtest should be logical of length dim(X)[2]")
	if (dim(X)[1] != length(Y))
		stop("dimensions of X and Y do not match")
	if (any(is.na(Y)) || any(is.na(X))) {
		warning("missing data points in Y or X, dropping")
		IScomplete <- !is.na(Y) & complete.cases(X)
		Y <- Y[IScomplete]
		X <- X[IScomplete,]
		if (!is.null(W)) W <- W[IScomplete,IScomplete]
	}
	
	# precompute tXW
	if (is.null(W)) {
		tXW <- t(X)
	} else {
		tXW <- t(X) %*% W
	}
	#print(tXW)
	
	# estimate beta
	XpWXm1 <- ginv(tXW %*% X)
	XpWY <- tXW %*% Y
	betaW <- XpWXm1 %*% XpWY
	#print(XpWXm1)
	
	# estimate V
	if (test=="wald") {
		YmXB <- Y - X %*% betaW
		if (is.null(W)) {
			sigma2 <- as.vector(((t(YmXB) %*% YmXB)/(dim(X)[1] - dim(X)[2])))
		} else {
			sigma2 <- as.vector(((t(YmXB) %*% W %*% YmXB)/(dim(X)[1] - dim(X)[2])))
		}
		#print(sigma2)
		V <- sigma2*XpWXm1
	} else if (test=="score") {
		Xm <- X[,!whichtest,drop=FALSE]
		if (is.null(W)) {
			bt <- ginv(t(Xm) %*% Xm) %*% (t(Xm) %*% Y)
			YmXB <- Y - Xm %*% bt
			sigma2 <- as.vector(((t(YmXB) %*% YmXB)/(dim(Xm)[1] - dim(Xm)[2])))
		} else {
			bt <- ginv(t(Xm) %*% W %*% Xm) %*% (t(Xm) %*% W %*% Y)
			YmXB <- Y - Xm %*% bt
			sigma2 <- as.vector(((t(YmXB) %*% W %*% YmXB)/(dim(Xm)[1] - dim(Xm)[2])))
		}
		V <- sigma2*XpWXm1
	} else if (test=="robust") {
		
		YmXB <- Y - X %*% betaW
		if (is.null(W)) {
			sigma2vec <- as.vector(YmXB * YmXB)
		} else {
			sigma2vec <- as.vector(as.vector(t(YmXB) %*% W) * YmXB)
		}
#  V <- XpWXm1 %*% (t(X) %*% W %*% (diag(as(YmXB,"vector")^2) %*% ginv(W)) %*% W %*% X) %*% XpWXm1
#  V <- XpWXm1 %*% (t(X) %*% W %*% diag(as(YmXB,"vector")^2) %*% X) %*% XpWXm1
#  V <- XpWXm1 %*% (t(X) %*% diag(as(YmXB,"vector")^2) %*% X) %*% XpWXm1
#  V <- ginv(t(X)%*%X) %*% (t(X) %*% diag(as(YmXB,"vector")^2) %*% X) %*% ginv(t(X)%*%X)
		##	V <- XpWXm1 %*% (tXW %*% diag(as(YmXB,"vector")^2) %*% X) %*% XpWXm1
		V <- XpWXm1 %*% (tXW %*% diag(sigma2vec) %*% X) %*% XpWXm1
# print(sigma2vec)
#		V <- XpWXm1 %*% (tXW %*% diag(sigma2vec) %*% t(tXW)) %*% XpWXm1
	} else {
		stop("test not recognised")
	}
	
	# do the test
	if (sum(whichtest)>1) {
		Vinv <- ginv(V)
		#print("V")
		#print(V)
		#print("Vinv")
		#print(Vinv)
		#print("det(V)")
		#print(det(V))
		#print(log(abs(det(V))))
		Vbg <- Vinv[whichtest,whichtest]-Vinv[whichtest,!whichtest] %*% ginv(Vinv[!whichtest,!whichtest]) %*% Vinv[!whichtest,whichtest]
		#Vbg <- Vinv[whichtest,whichtest]-Vinv[whichtest,!whichtest] %*% V[!whichtest,!whichtest,drop=F] %*% Vinv[!whichtest,whichtest]
		T2 <- t(betaW[whichtest]) %*% Vbg %*% betaW[whichtest]
	} else if (sum(whichtest) == 1) {
		T2 <- betaW[whichtest]^2/diag(V)[whichtest]
	} else {
		stop("unreachable statement")
	}
	
	rownames(betaW) <- colnames(X)
	dimnames(V) <- list(colnames(X),colnames(X))
	out <- list(beta=betaW,V=V,T2=T2,df=sum(whichtest),tested=whichtest,
			meanColX = apply(X,FUN=mean,MAR=2), n = dim(X)[1])
	class(out) <- "FGLS"
	out
}
#' summary for FGLS
#' 
#' summary for Feasible GLS
#' 
#' @param x object returned by FGLS function
#' @param verbosity what to report; 0 -- only test stats, 
#' 1 -- test stats and estimates concerning 'whichtest' parameters, 
#' 2 -- test stats and estimates of all parameters
#' @param varcov wheather var-cov matrix for estimated parameters 
#' is to be reported
#' @param include.means weather mean values of predictors are to 
#' be reported
#' 
#' @return one-line summary of output from FGLS function
#' 
#' @author Yurii Aulchenko
#'
summaryFGLS <- function(x,verbosity=1,varcov=FALSE,include.means=FALSE) {
	out <- c(x$n,x$T2,x$df,pchisq(x$T2,df=x$df,low=FALSE))
	names(out) <- c("n","Chisq","df","P-value")
	if (verbosity == 0) {
		return(out)
	} else if (verbosity == 1) {
		toout <- x$tested
	} else if (verbosity == 2) {
		toout <- rep(TRUE,length(x$tested))
	}
	
	if (is.null(colnames(x$V))) {
		xnames <- paste("pred",1:dim(x$V)[2],sep="")
	} else {
		xnames <- colnames(x$V)[toout]
	}
	names.save <- names(out)
	out <- c(out,x$beta[toout],sqrt(diag(x$V))[toout])
	names(out) <- c(names.save,paste("beta_",xnames,sep=""),
			paste("se_",xnames,sep=""))
	
	if (include.means) {
		names.save <- names(out)
		out <- c(out,x$meanColX[toout])
		names(out) <- c(names.save,paste("mean_",xnames,sep=""))
	}
	if (varcov) {
		Vout <- x$V[toout,toout]
		if (length(Vout) > 1) {
			names.Vout <- outer(rownames(Vout),colnames(Vout),FUN=paste,sep="_")
			Vout <- Vout[lower.tri(Vout)]
			names.Vout <- names.Vout[lower.tri(names.Vout)]
			names.save <- names(out)
			out <- c(out,Vout)
			names(out) <- c(names.save,paste("cov_",names.Vout,sep=""))
		}
	}
	return(out)
}
#' Genome-wide FGLS
#' 
#'  Genome-wide Feasible GLS
#' 
#' @param formula analysis formula; should contain special 'SNP' term
#' @param data phenotypic data frame, or 'gwaa.data-class' object
#' @param subset subset of data (individuals)
#' @param weights RESERVED FOR FURTURE USE
#' @param na.action RESERVED FOR FURTURE USE
#' @param contrasts RESERVED FOR FURTURE USE
#' @param offset RESERVED FOR FURTURE USE
#' @param W for GLS, (inverse of) variance-covariance matrix, as such returned by 
#' GenABEL's polygenic(...)\$InvSigma, or NULL for LS
#' @param inverse wheather W is already inverted 
#' @param na.SNP how to deal with missing SNP data; 'impute' -- substitute 
#' with mean, 'drop' -- drop rows with missing genotypes
#' @param mincall minimall call rate for a SNP (if less, the SNP is dropped)
#' @param residuals use residuals for analysis? (faster, but less precise)
#' @param test test to be applied, one of 'wald', 'score' or 'robust'
#' @param model.SNP SNP model to apply, one of 
#' c("additive","dominantB","recessiveB","overdominant","genotypic")
#' @param genodata genotypic data; can be missing if data is of 'gwaa.data-class'.
#' Otherwise can be regular matrix or 'databel' matrix
#' @param gtcoding one of c("typed","dose","probability")
#' 'typed' -- coded with NA, 0, 1, 2
#' @param verbosity what to report; 0 -- only test stats, 
#' 1 -- test stats and estimates concerning 'whichtest' parameters, 
#' 2 -- test stats and estimates of all parameters
#' @param varcov wheather var-cov matrix for estimated parameters 
#' is to be reported
#' @param include.means weather mean values of predictors are to 
#' be reported
#' @param singular what to do with linear dependencies in X (now only 
#' 'ignore' implemented)
#' @param with.lm whether LM should be run along (only test purposes; always 
#' falls back to 'old' R-only implementation)
#' @param old if TRUE, old R-only code implementation is running (testing purposes, 
#' slow)
#' 
#' @examples 
#' library(MASS)
#' library(mvtnorm)
#' library(GenABEL)
#' data(ge03d2.clean)
#' NIDS <- 100
#' df <- ge03d2.clean[1:NIDS,autosomal(ge03d2.clean)]
#' s <- summary(df@@gtdata)
#' maf <- pmin(s$Q.2,1-s$Q.2)
#' df <- df[,(maf>0.05)]
#' gkin <- ibs(df[,autosomal(df)],w="freq")
#' 
#' modelh2 <- 0.8
#' covars <- c("sex","age")
#' s2 <- 12^2
#' betas_cov <- c(170,12,0.01)
#' 
#' X <- as.matrix(cbind(rep(1,NIDS),phdata(ge03d2.clean)[1:NIDS,covars]))
#' grel <- gkin
#' grel[upper.tri(grel)] <- t(grel)[upper.tri(grel)]
#' grel <- 2*grel
#' Y <- as.vector(rmvnorm(1,mean=(X \%*\% betas_cov),sigma=s2*(modelh2*grel+diag(dim(grel)[1])*(1-modelh2))))
#' length(Y)
#' Y[2] <- NA
#' df@@phdata$Y <- Y
#' 
#' mdl <- Y~age+sex
#' h2 <- polygenic(mdl,data=df,kin=gkin)
#' h2$h2an
#' mdl <- Y~age+sex+SNP
#' 
#' gtOld <- df@@gtdata[,1:10]
#' gtReal <- as.double(gtOld)
#' gtNew <- as(gtReal,"databel")
#' aIR <- GWFGLS(mdl,data=phdata(df),genodata = gtReal,verbosity=0) #,model="genotypic")
#' aIO <- GWFGLS(mdl,data=phdata(df),genodata = gtOld, verbosity=0) #,model="genotypic")
#' aIN <- GWFGLS(mdl,data=phdata(df),genodata = gtNew, verbosity=0) #,model="genotypic")
#' aWR <- GWFGLS(mdl,data=phdata(df),W=h2$InvSigma,genodata = gtReal,verbosity=0) #,model="genotypic")
#' aWO <- GWFGLS(mdl,data=phdata(df),W=h2$InvSigma,genodata = gtOld, verbosity=0) #,model="genotypic")
#' aWN <- GWFGLS(mdl,data=phdata(df),W=h2$InvSigma,genodata = gtNew, verbosity=0) #,model="genotypic")
#' aIN
#' aWN
#' complex_model <- GWFGLS(Y~age+sex*SNP,data=phdata(df),W=h2$InvSigma,
#' 	genodata = gtReal,verbosity=2,varcov=TRUE,include.means=TRUE,model="genotypic")
#' complex_model
#' 
#' @author Yurii Aulchenko
#' 
#'
GWFGLS <- 
		function (formula, data, 
				subset, weights, na.action, 
				contrasts = NULL, offset, 
				W = NULL, inverse = TRUE,
				na.SNP = "impute", mincall = 0.95, residuals = FALSE, test = "wald", 
				model.SNP = "additive",
				genodata,gtcoding = "typed",verbosity=1, varcov = FALSE, 
				include.means = TRUE, singular = "ignore", with.lm = FALSE,
				old=FALSE)
# test in c("wald","score","robust")
# model.SNP in c("additive","dominantB","recessiveB","overdominant","genotypic")
# verbosity 0 : only chi2 and df, 1: only SNP-related; 2 -- all betas and se of betas
# varcov: TRUE -- report var-cov matrix
# genodata in c("matrix","data.frame","gtdata","databel")
# gtcoding in c("typed","dose","probability")
# singular: what to do with linear dependencies in X? in c("ignore","drop")
# mincall: do not analyze the SNP if call is less (<) then mincall
{
	
	if (gtcoding == "dose" && model.SNP != "additive") stop("when gtcoding is 'dose', only 'additive' model allowed")
# how to check consistency between genodata & gtcoding???
	
	allstattests <- c("wald","score","robust")
	intstattests = pmatch(test, allstattests, nomatch = -1, duplicates.ok = FALSE)
	if (intstattests < 0 || intstattests > length(allstattests)) 
		stop(paste("stat test",test,"not in any of ",allstattests,"\n"))
	
	
# start with checking consistency between phen- and geno-data
	if (is(data,"gwaa.data")) {
		if (!missing(genodata)) warning("data has gwaa.data-class; dropping data provided by genodata argument")
		genodata <- gtdata(data)
		data <- phdata(data)
	} else {
		if (missing(genodata)) stop("genodata argument missing with no default")
	}
	if (any(dim(data)[1] != dim(genodata)[1])) stop("dimensions of data and genodata do not match")
	
	
# prepare pheno-data & filter geno-data
	mydata <- 
			deliver.data(formula = formula, data = data) #, subset = subset, weights = weights, na.action = na.action, 
	#contrasts = contrasts, offset = offset) 
	
	genodata <- genodata[mydata$keep_IDs,,drop=FALSE]
	
	if (residuals) {
		mydata$Y <- lm.fit(x=mydata$Xx[,mydata$elXx_NOTin_Xg,drop=FALSE],y=mydata$Y,offset=mydata$offset)$residuals
		mydata$Xx <- mydata$Xx[,mydata$elXx_in_Xg,drop=FALSE]
		mydata$Xx <- cbind(1,mydata$Xx)
		dimnames(mydata$Xx)[[2]][1] <- "(Intercept)"
	}
	
# figure out the genotypic data type 
	charcoding <- c("typed","dose","probability")
	if (length(gtcoding) != 1) {
		cat("gtcoding argument (now",gtcoding,") should be of length one, one of",charcoding,"\n")
		stop()
	}
	intcoding = pmatch(gtcoding, charcoding, nomatch = -1, duplicates.ok = FALSE)
	if (intcoding < 0 || intcoding > length(charcoding)) stop(paste("coding",gtcoding,"not in any of ",charcoding,"\n"))
	gtcoding <- charcoding[intcoding]
	step <- 1
	if (gtcoding == "probability") step <- 2
	
# figure out the model
	genomodel <- c("additive","dominantB","recessiveB","overdominant","genotypic")
	if (length(model.SNP) != 1) {
		cat("model.SNP argument (now",model.SNP,") should be of length one, one of",genomodel,"\n")
		stop()
	}
	intmodelcoding = pmatch(model.SNP, genomodel, nomatch = -1, duplicates.ok = FALSE)
	if (intmodelcoding < 0 || intmodelcoding > length(genomodel)) stop(paste("coding",model.SNP,"not in any of ",genomodel,"\n"))
	model.SNP <- genomodel[intmodelcoding]
	nColG <- 1
	if (model.SNP == "genotypic") nColG <- 2
	
	
# prepare genotypic data
#  if (!is.null(names(mydata$Y))) idnamesY <- names(mydata$Y)
#  else if (!is.null(dimnames(mydata$Y)[[1]])) idnamesY <- dimnames(mydata$Y)[[1]]
#  else idnamesY <- NULL
	
# Prepare final data
	finaldata <- list()
	finaldata$Y <- mydata$Y
	if (nColG==1) {
		finaldata$X <- cbind(mydata$Xx,mydata$Xg)
		finaldata$snpinvolved <- c(rep(0,dim(mydata$Xx)[2]),rep(1,dim(mydata$Xg)[2]))
	} else if (nColG == 2) {
		names.save <- colnames(mydata$Xg)
		colnames(mydata$Xg) <- paste(names.save,".het",sep="")
		finaldata$X <- cbind(mydata$Xx,mydata$Xg)
		colnames(mydata$Xg) <- paste(names.save,".hom",sep="")
		finaldata$X <- cbind(finaldata$X,mydata$Xg)
		finaldata$snpinvolved <- c(rep(0,dim(mydata$Xx)[2]),rep(1,dim(mydata$Xg)[2]),rep(2,dim(mydata$Xg)[2]))
	} else {
		stop("nColG != 1, != 2")
	}
	#print(finaldata)
	
# pre-compute different staff
	
	fixed_df <- sum(finaldata$snpinvolved>0)
	
	resid.score <- lm.fit(mydata$Xx,mydata$Y)$resid
	if (class(mydata$Y) == "matrix") {
		stop("matrix-Y not implemented")
		#sigma2.score <- apply(resid.score^2,FUN=sum,MAR=2)/(dim(mydata$Xx)[1] - dim(mydata$Xx)[2])
	} else {
		Xm <- mydata$Xx
		Y <- finaldata$Y
		if (is.null(W)) {
			bt <- ginv(t(Xm) %*% Xm) %*% (t(Xm) %*% finaldata$Y)
			YmXB <- Y - Xm %*% bt
			sigma2.score <- as.vector(((t(YmXB) %*% YmXB)/(dim(Xm)[1] - dim(Xm)[2])))
		} else {
			bt <- ginv(t(Xm) %*% W %*% Xm) %*% (t(Xm) %*% W %*% finaldata$Y)
			YmXB <- Y - Xm %*% bt
			sigma2.score <- as.vector(((t(YmXB) %*% W %*% YmXB)/(dim(Xm)[1] - dim(Xm)[2])))
		}
		rm(Xm,YmXB,bt,Y);gc()
	}
	#print(sigma2.score)
	
	if (!is.null(W)) {
		require(MASS)
		if (!inverse) {
			V <- W
			W <- ginv(W)
		} else {
			V <- ginv(W)
		}
		tXW.fixed <- t(finaldata$X) %*% W
	} else {
		tXW.fixed <- t(finaldata$X)
	}
	
# run 1 time to get output dimensions and names
	mgls <- FGLS(Y=finaldata$Y,X=finaldata$X,whichtest=(finaldata$snpinvolved>0))
	smr <- summaryFGLS(mgls,verbosity=verbosity,varcov=varcov,include.means=include.means)
	#print(smr)
	
# NOW ITERATE OVER SNPs
	
	
	out <- list()
	snpseq <- seq(1,dim(genodata)[2],step) 
	outsnpnames <- c(dimnames(genodata)[[2]][snpseq])
	if (is.null(outsnpnames)) outsnpnames <- paste("mrk",c(1:length(snpseq)),sep="")
	if (with.lm) {
		out$LM <- rep(NA,length(snpseq))
		warning("with.lm specified; falling to 'old' R-only implementation (slow)");
		old <- TRUE
	}
	
	if (old) {
		out$FGLS <- matrix(NA,ncol=length(smr),nrow=length(snpseq))
		#out$FGLS <- matrix(NA,ncol=2*dim(finaldata$X)[2]+2,nrow=length(snpseq))
		dimnames(out$FGLS) <- list(outsnpnames,names(smr))
		#c(paste("beta_",colnames(finaldata$X),sep=""),paste("se_",colnames(finaldata$X),sep=""),"Chi2","P-value")
		
		#print("AAA")
		
		k <- 1
		for (csnp in snpseq) {
			
			GGG <- get_snp_data(genodata,csnp,model.SNP,gtcoding)
			nocalls <- sum(is.na(GGG))
			callrate <- (length(GGG)-nocalls)/length(GGG)
			
			if ( callrate >= mincall) {
				
				if (nocalls > 0) {
					if (na.SNP=="impute") {
						mn <- apply(GGG,MAR=2,FUN=mean,na.rm=TRUE)
						for (iii in 1:dim(GGG)[2]) GGG[is.na(GGG[,iii]),iii] <- mn[iii]
					} else if (na.SNP == "drop") {
						#stop("na.SNP != 'impute' and nocalls > 0; 'drop' not yet implemented")
						nomiss <- complete.cases(GGG)
						finaldata$Y <- finaldata$Y[nomiss]
						finaldata$X <- finaldata$X[nomiss,]
						GGG <- GGG[nomiss,]
					} else {
						stop(paste("na.SNP argument not recognised:",na.SNP))
					}
				}
				XXX <- finaldata$X
				for (i in which(finaldata$snpinvolved==1)) XXX[,i] <- XXX[,i]*GGG[,1]
				for (i in which(finaldata$snpinvolved==2)) XXX[,i] <- XXX[,i]*GGG[,2]
				if (singular == "drop") {
					# identify collinear and mono columns
					# probably use gsl_linalg_QRPT_decomp? or other QR with pivoting?
				} else {
					current_df <- fixed_df
				} 
				#print(c("NOW SNP",csnp))
				#print(cbind(mydata$Y,XXX))
				mgls <- FGLS(Y=mydata$Y,X=XXX,W=W,test=test,whichtest=(finaldata$snpinvolved>0))
				out$FGLS[k,] <- summaryFGLS(mgls,verbosity=verbosity,varcov=varcov,include.means=include.means)
				#out$FGLS[k,] <- c(mgls$beta,sqrt(diag(mgls$V)),mgls$T2,pchisq(mgls$T2,df=current_df,low=F))
				#print(mgls$beta)
				#print(pchisq(mgls$T2,df=sum(finaldata$snpinvolved>0),low=F))
				#print(det(t(XXX)%*%XXX))
				#print(rankMatrix(XXX))
				#print(qr(XXX))
				#print(svd(t(XXX)%*%XXX))
				if (with.lm) {
					lm0 <- lm(mydata$Y ~ mydata$Xx)
					lm1 <- lm(mydata$Y ~ 0+XXX)
					#print(summary(lm1))
					#print(anova(lm0,lm1,test="Chisq"))
					out$LM[k] <- anova(lm0,lm1,test="Chisq")$"P(>|Chi|)"[2]
				}
			} else {
				out$FGLS[k,] <- rep(NA,dim(out$FGLS[2]))
				if (with.lm) out$LM[k] <- NA
				warning(paste("call < mincall; skipping SNP",csnp))
			}
			
			k <- k + 1
		}
	} else {
		# new C++ procedure
		if (class(genodata) == "snp.data") {
			gtNrow <- dim(genodata)[1]
			gtNcol <- dim(genodata)[2]
			genodata <- genodata@gtps
		} else if (class(genodata) == "matrix") {
			gtNrow <- dim(genodata)[1]
			gtNcol <- dim(genodata)[2]
			storage.mode(genodata) <- "double"
		} else if (class(genodata)=="databel") {
			gtNrow <- dim(genodata)[1]
			gtNcol <- dim(genodata)[2]
			genodata <- genodata@data
		} else {
			stop(paste("genodata class not recognised ('",class(genodata),"')",sep=""))
		}
		
		# check NaNs and create gsl_matrix data
		if (any(is.na(finaldata$Y))) stop("NANs not allowed in Y")
		if (any(is.na(finaldata$X))) stop("NANs not allowed in X")
		if (!is.null(tXW.fixed)) if (any(is.na(tXW.fixed))) stop("NANs not allowed in tXW.fixed")
		if (!is.null(W)) if (any(is.na(W))) stop("NANs not allowed in W")
		SEXP_ptr_to_gsl_Y <- .Call("gslVector_from_SEXP",finaldata$Y)
		SEXP_ptr_to_gsl_X <- .Call("gslMatrix_from_SEXP",finaldata$X)
		if (!is.null(tXW.fixed)) 
		{SEXP_ptr_to_gsl_tXW_fixed <- .Call("gslMatrix_from_SEXP",tXW.fixed);} 
		else {SEXP_ptr_to_gsl_tXW_fixed <- NULL}
		if (!is.null(W)) {SEXP_ptr_to_gsl_W <- .Call("gslMatrix_from_SEXP",W);} else {SEXP_ptr_to_gsl_W <- NULL;}
		
		#mgls <- FGLS(Y=finaldata$Y,X=finaldata$X,whichtest=(finaldata$snpinvolved>0))
		#print("before fake_iterator")
		#print(finaldata$X)
		resIteratorNew <- .Call("iterator", genodata, 
				as.integer(gtNrow), as.integer(gtNcol),
				as.character("fgls"),
				"R", as.integer(2), 
				as.integer(1), # STEP 
				as.integer(10), 
				SEXP_ptr_to_gsl_Y, SEXP_ptr_to_gsl_X,
				SEXP_ptr_to_gsl_tXW_fixed, SEXP_ptr_to_gsl_W,
				as.integer(intstattests),
				as.integer(finaldata$snpinvolved),
				as.integer(which(model.SNP==genomodel)),
				as.integer(which(gtcoding==charcoding)),
				as.double(mincall),
				as.double(sigma2.score))
		#if (!is(resIteratorNew,"matrix")) resIteratorNew <- matrix(resIteratorNew,nrow=1)
		#print("after fake_iterator")
		mgls <- FGLS(Y=finaldata$Y,X=finaldata$X,whichtest=(finaldata$snpinvolved>0))
		smr <- summaryFGLS(mgls,verbosity=2,varcov=TRUE,include.means=TRUE)	
		#print(smr)
		#print(resIteratorNew)
		outdf <- sum(finaldata$snpinvolved != 0)
		outN <- dim(finaldata$X)[1]
		nBeta <- dim(finaldata$X)[2]
		tmpMtr <- matrix(NA,ncol=nBeta,nrow=nBeta)
		tmpMtr[lower.tri(tmpMtr,diag=T)] <- 1:(nBeta*(nBeta+1)/2)
		varSeq <- diag(tmpMtr)
		rIsize <- dim(resIteratorNew)[2]
		#print(rIsize)
		vcoSeq <- as.vector(tmpMtr[lower.tri(tmpMtr)])
		#print("before cbind")
		out$FGLS <- cbind(outN,resIteratorNew[,1,drop=FALSE],
				outdf,pchisq(resIteratorNew[,1,drop=FALSE],df=outdf,low=FALSE),
				resIteratorNew[,c(2:(1+nBeta)),drop=FALSE],
				sqrt(resIteratorNew[,c((1+nBeta)+varSeq),drop=FALSE]),
				resIteratorNew[,c((rIsize-nBeta+1):rIsize),drop=FALSE],
				resIteratorNew[,c(1+nBeta+vcoSeq),drop=FALSE]
		)
		#if (!is(out$FGLS,"matrix")) out$FGLS <- matrix(out$FGLS,nrow=1)
		#print(dim(out$FGLS))
		#print(out$FGLS)
		#print(smr)
		dimnames(out$FGLS)[[2]] <- names(smr)
		smr <- summaryFGLS(mgls,verbosity=verbosity,varcov=varcov,include.means=include.means)
		out$FGLS <- out$FGLS[,names(smr),drop=FALSE]
		#print(out$FGLS)
		dimnames(out$FGLS)[[1]] <- outsnpnames
		out$FGLS[which(out$FGLS[,"Chisq"]<(-0.1)),"P-value"] <- NA
		#print(out$FGLS)
		#print("before gc")
		gc();
		#print("after gc")
	}
	
	if (with.lm) return(out)
	else return(out$FGLS)
} 
