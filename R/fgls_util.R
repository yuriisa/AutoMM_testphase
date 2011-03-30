get_snp_data <- function(genodata,csnp,model.SNP,gtcoding)
{
	if (gtcoding == "dose" && model.SNP != "additive") 
		stop("'dose' type only possible with 'additive' model")
	if (gtcoding == "probability") {
		if (model.SNP=="additive") return(matrix(as(genodata[,csnp],"double")*2+as(genodata[,csnp],"double")*1,ncol=1))
		GGG <- matrix(c(as(genodata[,csnp+1],"double"),as(genodata[,csnp],"double")),ncol=2)
	} else if (gtcoding == "typed") {
		if (model.SNP=="additive") return(matrix(as(genodata[,csnp],"double"),ncol=1))
		tmp <- as(genodata[,csnp],"double")
		tmp <- as.integer(tmp)
		GGG <- matrix(c(1*(tmp==1),1*(tmp==2)),ncol=2)
	} else if (gtcoding == "dose") {
		if (model.SNP=="additive") return(matrix(as(genodata[,csnp],"double"),ncol=1))
		else stop("uhaha")
	} else {
		stop("brrrr")
	}
	if (model.SNP == "genotypic") return(GGG)
	else if (model.SNP == "dominantB") return(GGG[,1,drop=FALSE]+GGG[,2,drop=FALSE])
	else if (model.SNP == "recessiveB") return(GGG[,2,drop=FALSE])
	else if (model.SNP == "overdominant") return(GGG[,1,drop=FALSE])
	else stop("avava")
}

deliver.data <- 
		function (formula, data, subset, weights, na.action, 
				contrasts = NULL, offset) 
{
	#ret.x <- x
	#ret.y <- y
	cl <- match.call()
	#print(cl)
	mf <- match.call() #expand.dots = FALSE)
	#print(mf)
	m <- match(c("formula", "data", "subset", "weights", "na.action", 
					"offset"), names(mf), 0L)
	#print(m)
	mf <- mf[c(1L, m)]
	#print(mf)
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	#print(names(mf))
	#mf$formula <- as.formula(mf$formula)
	mf.save <- mf
	#print(class(mf$formula))
	#print(mf$formula)
	
	mf$formula <- eval(mf$formula)
	vrs <- as.character(attr(terms(mf$formula),"variables"))[-1] #all.vars(mf$formula)
	#print(vrs)
	if (!any(vrs=="SNP")) stop("no SNP term")
	notSNP <- (vrs!="SNP")
	vrs <- vrs[notSNP]
	
	if (sum(!notSNP)) {
		mfI <- mf
		#print(mf)
		#print(mfI)
		fla <- paste(vrs[1],"~")
		if (length(vrs)==1) vrs[2] <- "1"
		fla <- paste(fla,vrs[2])
		if (length(vrs)>2) for (i in 3:length(vrs))
				fla <- paste(fla,vrs[i],sep="+")
		#print(fla)
		mfI$formula <- as.formula(fla)
		#print(mfI$formula)
		#print(mfI)
		mfI.save <- mfI
		mfI <- eval(mfI, parent.frame())
		mfI.save$na.action <- na.pass
		tmp <- eval(mfI.save)
		#print(mf.save)
		keep_IDs <- rownames(tmp) %in% rownames(mfI)
		mfI$SNP <- 1
		#print("MF0")
		#print(mfI)
		#print(class(mfI))
		mf.save$data <- mfI
		#print(mf.save)
		#print(mf.save$formula)
		# somewhat messy way to access data 
		#attach(mfI,warn.conflicts = TRUE)
		#mf$formula <- eval(mf$formula)
		mf <- eval(mf.save, envir=mfI)
		#detach(mfI)
		#print("MF")
		#print(mf)
	} else {
		mf.save <- mf
		mf <- eval(mf, parent.frame())
		mf.save$na.action <- na.pass
		tmp <- eval(mf.save, envir=mfI)
		keep_IDs <- rownames(tmp) %in% rownames(mf)
	}
	
	
	#print("NOW MF")
	#print(class(mf))
	#print(attributes(mf))
	#print(attr(attr(mf, "terms"),"factors"))
	#if (method == "model.frame") 
	#    return(mf)
	#else if (method != "qr") 
	#    warning(gettextf("method = '%s' is not supported. Using 'qr'", 
	#        method), domain = NA)
	mt <- attr(mf, "terms")
	#print("MT")
	#print(mt)
	y <- model.response(mf, "numeric")
	#SNP <- rep(1,length(y))
	#print(y)
	w <- as.vector(model.weights(mf))
	#print(w)
	if (!is.null(w) && !is.numeric(w)) 
		stop("'weights' must be a numeric vector")
	offset <- as.vector(model.offset(mf))
	#print(offset)
	if (!is.null(offset)) {
		if (length(offset) != NROW(y)) 
			stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
							length(offset), NROW(y)), domain = NA)
	}
	if (is.empty.model(mt)) {
		xx <- NULL
		xg <- NULL
		#z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 
		#    3) else numeric(0L), residuals = y, fitted.values = 0 * 
		#    y, weights = w, rank = 0L, df.residual = if (is.matrix(y)) nrow(y) else length(y))
		#if (!is.null(offset)) {
		#    z$fitted.values <- offset
		#    z$residuals <- y - offset
		#}
	}
	else {
		x <- model.matrix(mt, mf, contrasts)
		#print(x)
		#print(attr(mt,"factors"))
		#print("aaa")
		attSNP <- (attr(mt,"factors")["SNP",] != 0)
		#print(attr(mt,"intercept")==0)
		#print(attSNP)
		
		#print("SUBM")
		subm <- attr(mt,"factors")[-which(dimnames(attr(mt,"factors"))[[1]]=="SNP"),,drop=FALSE]
		#print(subm)
		
		inNam <- c()
		for (iii in which(attSNP) ) {
			cmpvec <- subm[,iii]
			for (jjj in which(!attSNP)) {
				#print(c(iii,jjj))
				#print(cmpvec)
				#print(subm[,jjj])
				if (all(subm[,jjj]==cmpvec)) {
					inNam <- c(inNam,dimnames(subm)[[2]][jjj])
					#print(c("AAAAAA",dimnames(subm)[[2]][jjj]))
				}
				#print(inNam)
			}
		}
		#print(inNam)
		
		if (attr(mt,"intercept")==0) { 
			xx <- x[,!attSNP,drop=FALSE]
			xg <- x[,attSNP,drop=FALSE]
		} else {
			xx <- x[,c(TRUE,!attSNP),drop=FALSE]
			xg <- x[,c(FALSE,attSNP),drop=FALSE]
		}
		
		if (length(which(!(dimnames(xx)[[2]] %in% inNam)))>0)
			elXx_NOTin_Xg <- which(!(dimnames(xx)[[2]] %in% inNam))
		else
			elXx_NOTin_Xg <- NULL
		
		if (length(which(dimnames(xx)[[2]] %in% inNam))>0)
			elXx_in_Xg <- which(dimnames(xx)[[2]] %in% inNam)
		else
			elXx_in_Xg <- NULL
		
		#z <- if (is.null(w)) 
		#    lm.fit(x, y, offset = offset, singular.ok = singular.ok, 
		#        ...)
		#else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, 
		#    ...)
	}
	#class(z) <- c(if (is.matrix(y)) "mlm", "lm")
	#z$na.action <- attr(mf, "na.action")
	#z$offset <- offset
	#z$contrasts <- attr(x, "contrasts")
	#z$xlevels <- .getXlevels(mt, mf)
	#z$call <- cl
	#z$terms <- mt
	#if (model) 
	#    z$model <- mf
	#if (ret.x) 
	#    z$x <- x
	#if (ret.y) 
	#    z$y <- y
	#if (!qr) 
	#    z$qr <- NULL
	#z
	out <- list(Y=y,Xx=xx,Xg=xg,offset=offset,elXx_NOTin_Xg=elXx_NOTin_Xg,elXx_in_Xg=elXx_in_Xg,keep_IDs = keep_IDs)
	out
}
