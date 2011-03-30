get_snp_data <- function(genodata,csnp,model.SNP,gtcoding)
{
	if (gtcoding == "dose" && model.SNP != "additive") stop("brrA")
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


make_random_matrix <- function(range_dim1 = c(50,500), range_dim2 = c(50,500), range_data = c(-1e16,1e16), type="double")
#make_random_matrix <- function(range_dim1 = c(200,1000), range_dim2 = c(200,1000), range_data = c(-1e16,1e16), type="double")
#make_random_matrix <- function(range_dim1 = c(500,1919), range_dim2 = c(1000,5000), range_data = c(-1e16,1e16), type="double")
{
	dim1 <- round(runif(1,range_dim1[1],range_dim1[2]))
	dim2 <- round(runif(1,range_dim2[1],range_dim2[2]))
	data <- runif(dim1*dim2,range_data[1],range_data[2])
	data <- as(data,type)
	data <- matrix(data,nrow=dim1,ncol=dim2)
	namesCol <- paste("col",c(1:dim2),sep="_")
	namesRow <- paste("row",c(1:dim1),sep="_")
	dimnames(data) <- list(namesRow,namesCol)
	return(data)
}

checkNumEq <- function(testmatr,test_fv,tolmult=5)
{
    print("CheckNumEq()");
    matr <- as(test_fv,"matrix")
    print("testmatr = ");
    show(testmatr);
    print("matr = ");
    show(matr)
    print("test_fv = ");    
    show(test_fv)
	checkEqualsNumeric(testmatr,matr,tolerance=tolmult*sqrt(.Machine$double.eps))
	# not 1x1 matrix ?
	if (!is.null(dim(testmatr))) {
	    checkIdentical(dim(testmatr),dim(test_fv))
	}
	dmn <- dimnames(test_fv)
	if (is.null(dmn)) dmn <- get_dimnames(test_fv)
	cat("dmn=")
	show(dmn)
	cat("dimnames(testmatr)=")
	show(dimnames(testmatr))
	checkIdentical(dimnames(testmatr),dmn)
}	
