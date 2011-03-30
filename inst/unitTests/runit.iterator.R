### --- Test setup ---

if(FALSE) {
	## Not really needed, but can be handy when writing tests
	library("RUnit")
	library("GenABEL")
	library("DatABEL")
	library("MixABEL")
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

skip_test.iterator_qtscore <- function() {
	
	# Witout RUnit, test with:
	# library(GenABEL); library(MixABEL); library(RUnit); path=getwd()
	# source('runit.iterator.R'); source('qtscore.R'); load('mc.RData')
	# res.score <- test.iterator_qtscore()
	
	# Test pocedure according to srdta demo file
	data(srdta); attach(srdta@phdata)
	# qtscore functions:
	source(paste(path, "/qtscore.R", sep=""))
	# Following file contains output from (file: GenABEL/demo/srdta.R, line: 183): 
	# mc <- check.marker(data=srdta,call=0.94,maf=(5/srdta@gtdata@nids),fdr=0.05,hweids=(srdta@phdata$bt==0),ibs.exclude="both")
	load(paste(path, '/mc.RData', sep=""))
	mymrk <- mc$snpok
	
	# Calling the iterator version of the qtscore function
	res.score <- qtscore(bt,data=srdta,snps=mymrk,trait="binomial")
	
	# Comparing with GenABEL version of qtscore, now fails
	checkIdentical(res.score, qtscore(bt,data=srdta,snps=mymrk,trait="binomial", package="GenABEL"))
	
	rm(list=ls());gc()
}

test.iterator_sum <- function()
{
#    testmatr <- make_random_matrix()
	
#	library("RUnit")
#	library("GenABEL")
#	library("DatABEL")
#	library("MixABEL")
	
	unlink("tmp*")
	
	data(srdta)
	
	nSn <- (sample(round(nsnps(srdta)/10):(round(nsnps(srdta)/5)),1))
	if ((nSn %% 2)!=0) nSn <- nSn - 1
	nIn <- (sample(round(nids(srdta)/10):(round(nids(srdta)/5)),1))
	if ((nIn %% 2)!=0) nIn <- nIn - 1
	dataOld <- srdta@gtdata[1:nIn,1:nSn]
	dataReal <- as.double(dataOld)
	dataNew <- as(dataReal,"databel")
	
	checkIdentical(dataReal,as.double(dataOld))
	checkIdentical(dataReal,as(dataNew,"matrix"))
	
	in_data <- list(dataNew@data,dataOld@gtps,dataReal) 
	out_data <- list("R","file")

	na_rm <- c(TRUE,FALSE)
	step <- c(1,2)
	mar <- c(1,2)
	
	nRow <- dim(dataReal)[1]
	nCol <- dim(dataReal)[2]
	
	for (cMar in mar) {
		for (cStep in step) {
			for (cOutData in out_data) {
				for (cNaRm in na_rm) {
					for (cInData in in_data) {
						
						if (cOutData != "R") # Create a new tmp file for each run
							cOutData = get_temporary_file_name()
						
						print("---------------------")
						print(c(
										class(cInData), 
										as.integer(nRow), as.integer(nCol),
										as.character("sum"),
										cOutData, as.integer(cMar), as.integer(cStep), 
										as.integer(1), as.integer(cNaRm)
								))
						resTemplate <- apply(dataReal,FUN=sum,MAR=cMar,na.rm=cNaRm)
						if (cStep>1) {
							tmp <- c()
							k <- 1;
							for (i in seq(1,length(resTemplate),cStep)) {
								tmp[k] <- 0
								for (j in 0:(cStep-1)) {tmp[k] <- tmp[k] + resTemplate[i+j]}
								k <- k+1
							}
							resTemplate <- tmp
						}
						print("Data as they should be: ")
						print(resTemplate)
						
						resIterator <- .Call("iterator", cInData, 
								as.integer(nRow), as.integer(nCol),
								as.character("sum"),
								cOutData, as.integer(cMar), as.integer(cStep), 
								as.integer(1), as.integer(cNaRm))
						
						# Iterator returns empty object when storing result in a file, read
						# in the data and perform checkEqual test
						if (cOutData != "R") {
							resIterator <- databel(cOutData)
							resIterator = as(resIterator, "vector")
						}
						
						print("Data from iterator: ")
						print(resIterator)
						
						checkEqualsNumeric(resTemplate, resIterator)
					}
				}
			}
		}
	}
	
	# ADD Fail-tests checkException(), e.g. missing last argument for sum, or similar
	# resIterator <- .Call("iterator", dataReal, 
	#								as.integer(nRow), as.integer(nCol),
	#								as.character("sum"),
	#								cOutData, as.integer(cMar), as.integer(cStep), 
	#								as.integer(0))
	
	rm(list=ls());gc()
	
	unlink("tmp*")
}
