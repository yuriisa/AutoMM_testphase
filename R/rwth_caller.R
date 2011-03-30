"rwth_caller" <- function(
		printall,
# parameters and data to compute M
		Phi, H2, Sigma2,
# vector of length length_y
		y,
# matrix of height length_y and width widthXL
		XL,
# matrix of height length_y and width nXRs*widthXRs
# it contains nXRs matrices XR
		XRs,
# width of single XR
		WidthXR
)
{
# figure out data dims and some sanity checks
	Length_y = length(y);
	NXRs = dim(XRs)[2]/WidthXR;
	if ((dim(XRs)[2] %% WidthXR) != 0) stop("(dim(XRs)[2] %% widthXR) != 0");
	if (dim(XRs)[1] != length(y)) stop("dim(XRs)[1] != length(y)");
	if (is.null(XL)) {
		XL <- matrix(1,ncol=1,nrow=length(y));
	} else {
		tmp <- matrix(1,ncol=1,nrow=length(y));
		XL <- cbind(tmp,XL)
	}
	if (dim(XL)[1] != length(y)) stop("dim(XL)[1] != length(y)");
	WidthXL = dim(XL)[2];
	
	ncol_beta <- WidthXL + WidthXR;
	nrow_out <- NXRs;
	ncol_vcbeta <- ncol_beta*(ncol_beta+1)/2
	
	out <- .C("rwth_example",
			as.integer(printall),
# dimensions of the data, explained later
			as.integer(Length_y), as.integer(NXRs),
			as.integer(WidthXR), as.integer(WidthXL),
			# parameters and data to compute M
			as.double(Phi), as.double(H2), as.double(Sigma2),
			# vector of length length_y
			as.double(y),
			# matrix of height length_y and width widthXL
			as.double(XL),
			# matrix of height length_y and width nXRs*widthXRs
			# it contains nXRs matrices XR
			as.double(XRs),
			# this is output
			# 'beta's length is (widthXL + widthXRs)
			# for time-being, just fill it with mean of
			# columns of X		
			beta = double( nrow_out * ncol_beta )
	# 'vcbeta' is square matrix with
	# width/height = (widthXL + widthXRs)
	# for now, will skip it
	# ...
	)$beta
	out <- matrix(out,ncol=ncol_beta)
	out
}