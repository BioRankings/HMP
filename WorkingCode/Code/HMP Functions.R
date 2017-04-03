

library(dirmult)	# example code functions and dirmult, equalTheta functions
library(MASS)		# ginv function for Xsc.statistics
library(doParallel) # parallelizing KL calc
library(gplots)		# KL heatmap
library(ggplot2)	# pretty plot
library(vegan) 		# shannon/simpson diversity


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### External
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### ~~~~~~~~~~~~~~~~~~~~~
### Generation functions
### ~~~~~~~~~~~~~~~~~~~~~
# reviewed 9/29/16
Multinomial <- function(Nrs, probs){
	if(missing(Nrs) || missing(probs))
		stop("Nrs and/or probs missing.")
	
	# Create the data from the rmultinom
	mData <- matrix(0, length(Nrs), length(probs))
	for(i in 1:length(Nrs))	
		mData[i,] <- stats::rmultinom(1, Nrs[i], probs)
	
	# Label the created data
	colnames(mData) <- paste("Taxa", 1:ncol(mData))
	rownames(mData) <- paste("Sample", 1:nrow(mData))
	
	return(mData)
}

# reviewed 9/29/16
Dirichlet.multinomial <- function(Nrs, shape){
	if(missing(Nrs) || missing(shape))
		stop("Nrs and/or shape missing.")
	
	# Create the data from the rmultinom
	dmData <- matrix(0, length(Nrs), length(shape))
	for(i in 1:length(Nrs))	
		dmData[i,] <- stats::rmultinom(1, Nrs[i], dirmult::rdirichlet(1, shape))
	
	# Label the created data
	colnames(dmData) <- paste("Taxa", 1:ncol(dmData))
	rownames(dmData) <- paste("Sample", 1:nrow(dmData))
	
	return(dmData)
}



### ~~~~~~~~~~~~~~~~~~~~~
### Other functions
### ~~~~~~~~~~~~~~~~~~~~~
# reviewed 9/29/16
C.alpha.multinomial <- function(data){
	if(missing(data))
		stop("data missing.")
	
	perNumReadsSubs <- rowSums(data)/sum(data)
	
	# Get T statistic
	Ts <- T.statistics(data)
	
	M.alpha <- diag(perNumReadsSubs)- as.matrix(perNumReadsSubs) %*% t(as.matrix(perNumReadsSubs))
	g <- sum(diag(M.alpha %*% M.alpha)) / sum(diag(M.alpha))
	df <- (ncol(data)-1)*((sum(diag(M.alpha)))^2) / (sum(diag(M.alpha %*% M.alpha)))
	
	# Get pvalue
	pval <- 1-pchisq(q=Ts/g, df=df, ncp=0, lower.tail=TRUE)
	
	GoF.test <- list("T statistics"=Ts, "p value"=pval)
	
	return(GoF.test)
}

# reviewed 9/30/16 
DM.MoM <- function(data){
	if(missing(data))
		stop("data missing.")
	
	
	pi.MoM <- colSums(data)/sum(data)
	theta.MoM <- weirMoM(data, pi.MoM)$theta
	gamma.MoM <- pi.MoM*((1-theta.MoM)/theta.MoM)
	
	# Set LL to Inf if we only have 1 sample
	if(nrow(data) == 1){
		loglikdm <- Inf
	}else{
		loglikdm <- loglikDM(data, gamma.MoM)
	}
	
	fit.MoM <- list(loglik=loglikdm, gamma=gamma.MoM, pi=pi.MoM, theta=theta.MoM)
	
	return(fit.MoM)
}

# reviewed 9/30/16
pioest <- function(group.data){
	if(missing(group.data))
		stop("data.groups missing.")
	
	# Check every data set has the same number of taxa
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa))
		stop("Every data set must have matching taxa, including pi0")
	
	numGroups <- length(group.data)
	
	# Pull out pi and calculate xsc
	pis <- matrix(0, numTaxa, numGroups)
	xscs <- rep(0, numGroups)
	for(i in 1:numGroups){
		tempData <- group.data[[i]]
		
		numReadsSubs <- rowSums(tempData)
		totalReads <- sum(tempData)
		pi <- colSums(tempData)/totalReads
		theta <- weirMoM(tempData, pi)$theta
		
		xscs[i] <- (theta * (sum(numReadsSubs^2)-totalReads) + totalReads) / totalReads^2
		pis[,i] <- pi
	}
	
	# Remove any 0 taxa
	pis <- pis[rowSums(pis)!=0,]
	
	# Calculate pi0
	pi0 <- rowSums(pis/xscs)/sum(1/xscs)
	
	names(pi0) <- colnames(group.data[[1]])
	
	return(pi0)
}

Kullback.Leibler <- function(group.data, plot=TRUE, parallel=FALSE, cores=3){
	if(missing(group.data))
		stop("data missing.")
	
	# Check the number of groups
	numGrps <- length(group.data)
	if(numGrps < 2)
		stop("At least 2 data sets are required.")
	
	# Check every data set has the same number of taxa
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa))
		stop("Every data set must have matching taxa.")
	
	# Make sure we have group names
	if(is.null(names(group.data))){
		grpNames <- paste("Data Set", 1:numGrps)
	}else{
		grpNames <- names(group.data)
	}
	
	group.data <- lapply(group.data, function(x) x+1)  # Add 1 so we don't ever get an all 0 comparison
	
	if(parallel){
		cl <- parallel::makeCluster(min(cores, numGrps)) 
		doParallel::registerDoParallel(cl)
		
		results <- foreach::foreach(i=1:numGrps, .combine=list, .multicombine=TRUE, .inorder=TRUE, .packages=c("dirmult")) %dopar%{
			mle.param <- dirmult::dirmult(group.data[[i]], trace=FALSE)
			return(mle.param)
		}
		parallel::stopCluster(cl)
	}else{
		results <- vector("list", numGrps)
		for(i in 1:numGrps)	
			results[[i]] <- dirmult::dirmult(group.data[[i]], trace=FALSE)
	}
	
	alpha <- lapply(results, function(x) x$gamma)
	names(alpha) <- grpNames
	LL.list <- mapply(function(x, y) loglikDM(x, y), x=group.data, y=alpha)
	
	KLmat <- matrix(0, numGrps, numGrps)
	for(i in 1:numGrps){
		for(j in 1:numGrps){
			ll <- loglikDM(group.data[[i]], alpha[[j]])
			KLmat[i, j] <- LL.list[i]- ll
		}
	} 
	colnames(KLmat) <- grpNames
	rownames(KLmat) <- grpNames
	
	if(plot){
		KLdist <- dist(KLmat, method="euclidean")
		
		gplots::heatmap.2(as.matrix(KLdist), dendrogram="both", Rowv=TRUE, Colv=TRUE, trace="none", symm=FALSE,
				main="Kullback-Leibler Divergences", margins=c(12,9), density.info="none")
	}
	
	return(KLmat)
}

Get.Diversities <- function(data){
	if(missing(data))
		stop("data is missing.")
	
	# Get base diversities
	numTaxa <- ncol(data)
	numSamps <- nrow(data)
	
	shan <- apply(data, 1, function(x) vegan::diversity(x, "shannon"))
	shanSD <- stats::sd(shan)
	shanDI <- mean(shan)
	
	simp <- apply(data, 1, function(x) vegan::diversity(x, "simpson"))
	simpSD <- stats::sd(simp)
	simpDI <- mean(simp)
	
	# Get confidence on those diversities
	CI <- 1.96*(shanSD/sqrt(numSamps))
	shanCI <- c(shanDI - CI, shanDI + CI)
	CI <- 1.96*(simpSD/sqrt(numSamps))
	simpCI <- c(simpDI - CI, simpDI + CI)
	
	# Get effective diversities
	shanEff <- exp(shan)
	shanSD <- stats::sd(shanEff)
	shanEffDI <- mean(shanEff)
	
	simpEff <- 1/(1-simp)
	simpSD <- stats::sd(simpEff)
	simpEffDI <- mean(simpEff)
	
	# Get confidence on effective diversities
	CI <- 1.96*(shanSD/sqrt(numSamps))
	shanEffCI <- c(shanEffDI - CI, shanEffDI + CI)
	CI <- 1.96*(simpSD/sqrt(numSamps))
	simpEffCI <- c(simpEffDI - CI, simpEffDI + CI)
	
	# Build table to return
	res <- data.frame(matrix(0, 3, 6))
	colnames(res) <- c("Value", "Lower CI", "Upper CI", "Eff Value", "Lower Eff CI", "Upper Eff CI")
	rownames(res) <- c("Richness", "Shannon", "Simpson")
	res[,1] <- c(numTaxa, shanDI, simpDI)
	res[,2] <- c(numTaxa, shanCI[1], simpCI[1])
	res[,3] <- c(numTaxa, shanCI[2], simpCI[2])
	res[,4] <- c(numTaxa, shanEffDI, simpEffDI)
	res[,5] <- c(numTaxa, shanEffCI[1], simpEffCI[1])
	res[,6] <- c(numTaxa, shanEffCI[2], simpEffCI[2])
	
	return(res)
}

# reviewed 9/30/16 & 10/6/16
# changed mod cramer 
Xmcupo.effectsize <- function(group.data){
	if(missing(group.data))
		stop("group.data missing.")
	
	# Check every data set has the same number of taxa
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa))
		stop("Every data set must have matching taxa")
	
	numGroups <- length(group.data)
	totalReads <- sum(sapply(group.data, sum))
	
	if(numTaxa < numGroups)
		stop("The number of taxa must be greater than the number of groups.")
	
	# Get the parameters for every group
	groupParameter <- lapply(group.data, function(x){
				# Calc pi, theta and the number of reads
				numReadsSubs <- rowSums(x)
				pi.MoM <- colSums(x)/sum(x)
				theta.MoM <- weirMoM(x, pi.MoM)$theta
				
				return(list(pi=pi.MoM, theta=theta.MoM, nrs=numReadsSubs))
			})
	
	# Calculate Xmcupo stats for base data
	Xmcupo <- Xmcupo.statistics(groupParameter)
	
	# Edit parameters to use the biggest difference between pis
	groupParameterMax <- groupParameter
	for(i in 1:numGroups){
		newPi <- rep(0, numTaxa)
		newPi[i] <- 1
		groupParameterMax[[i]]$pi <- newPi
	}
	
	# Calculate Xmcupo stats for biggest difference
	XmcupoMax <- Xmcupo.statistics(groupParameterMax)
	
	# Calculate Cramers
	CramerV	<- sqrt(Xmcupo/(totalReads*min(numTaxa-1, numGroups-1)))
	Mod.CramerV <- sqrt(Xmcupo/XmcupoMax)
	
	# Calculate pvalue
	pval <- 1-pchisq(q=Xmcupo, df=(numGroups-1)*(numTaxa-1), ncp=0, lower.tail=TRUE)
	
	result 	<- c("Chi-Squared"=Xmcupo, "Cramer Phi"=CramerV, "Modified-Cramer Phi"=Mod.CramerV, "P value"=pval)
	return(result)
}



### ~~~~~~~~~~~~~~~~~~~~~
### Plot functions
### ~~~~~~~~~~~~~~~~~~~~~
# reviewed 9/30/16
Barchart.data <- function(data, title="Taxa Proportions"){
	if(missing(data))
		stop("data missing.")
	
	dataProp <- apply(data, 1, function(x){x/sum(x)})
	
	barplot(dataProp, col=rainbow(ncol(data)), horiz=TRUE, 
			main=title, axisnames=FALSE, font.main=20, font.sub=16)
}

Est.PI <- function(group.data, useMLE=TRUE, main=NULL, plot=TRUE){
	if(missing(group.data))
		stop("group.data is missing.")
	
	if(is.null(main)){
		if(useMLE){
			main <- "PI Vectors: Using MLE"
		}else{
			main <- "PI Vectors: Using MoM"
		}
	}
	
	# Check the number of groups
	numGroups <- length(group.data)
	
	# Check every data set has the same number of taxa
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa))
		stop("Every data set must have matching taxa.")
	
	# Make sure we have group names
	if(is.null(names(group.data))){
		grpNames <- paste("Data Set", 1:numGroups)
	}else{
		grpNames <- names(group.data)
	}
	
	# Calculate the pi and error bars for each group
	allParams <- data.frame(matrix(0, 0, 6))		
	for(i in 1:numGroups){
		# Pull out a single group and check for taxa with 0 column sums (add 1 to subject 1)
		tempData <- group.data[[i]]
		
		badTaxa <- which(colSums(tempData) == 0)
		if(length(badTaxa) != 0)
			tempData[1, badTaxa] <- tempData[1, badTaxa] + 1
		
		tempParam <- data.frame(matrix(0, ncol(tempData), 6))
		
		# Get the MoM for every taxa
		fsum <- dirmult::dirmult.summary(tempData, dirmult::dirmult(tempData, trace=FALSE))
		fsum <- fsum[-nrow(fsum),]
		
		# Turn the summary into a data frame we can plot from
		tempParam[,1] <- rownames(fsum)
		tempParam[,2] <- grpNames[i]
		
		if(useMLE){
			tempParam[,3] <- fsum$MLE
			tempParam[,4] <- fsum$se.MLE
			tempParam[,5] <- fsum$MLE + 1.96*fsum$se.MLE
			tempParam[,6] <- fsum$MLE - 1.96*fsum$se.MLE
		}else{
			tempParam[,3] <- fsum$MoM
			tempParam[,4] <- fsum$se.MOM
			tempParam[,5] <- fsum$MoM + 1.96*fsum$se.MOM
			tempParam[,6] <- fsum$MoM - 1.96*fsum$se.MOM
		}
		
		allParams <- rbind(allParams, tempParam)
	}
	colnames(allParams) <- c("Taxa", "Group", "PI", "SE", "Upper", "Lower")
	
	# Make sure none of our error bars go over 100 or below 0
	allParams$Upper <- ifelse(allParams$Upper > 1, 1, allParams$Upper)
	allParams$Lower <- ifelse(allParams$Lower < 0, 0, allParams$Lower)
	
	# Factor the data so it stays in the right order
	allParams$Group <- factor(allParams$Group, levels=grpNames)
	allParams$Taxa <- factor(allParams$Taxa, levels=unique(colnames(group.data[[1]])))
	
	if(plot){
		print(ggplot2::ggplot(allParams, aes_string(y="PI", x="Taxa", colour="Group")) +
						geom_point() + theme(legend.position = "top") +
						geom_errorbar(aes_string(ymax="Upper", ymin="Lower")) +
						labs(title=main, y="PI Vector", x="") +
						theme(axis.text.x=element_text(hjust=1, angle=45, size=6)
						)
		)
	}
	
	return(allParams)
}

Plot.Abundance <- function(group.data, main="Group Abundance"){
	if(missing(group.data))
		stop("group.data is missing.")
	
	# Check the number of groups
	numGroups <- length(group.data)
	
	# Check every data set has the same number of taxa
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa))
		stop("Every data set must have matching taxa.")
	
	# Make sure we have group names
	if(is.null(names(group.data))){
		grpNames <- paste("Data Set", 1:numGroups)
	}else{
		grpNames <- names(group.data)
	}
	
	# Calculate the percent abundance and error bars for each group
	allParams <- data.frame(matrix(0, 0, 5))		
	for(i in 1:numGroups){
		# Pull out a single group and check for taxa with 0 column sums (add 1 to subject 1)
		tempData <- group.data[[i]]
		
		badTaxa <- which(colSums(tempData) == 0)
		if(length(badTaxa) != 0)
			tempData[1, badTaxa] <- tempData[1, badTaxa] + 1
		
		tempParam <- data.frame(matrix(0, ncol(tempData), 5))
		
		# Get percent abundance for every taxa
		tempData <- tempData/rowSums(tempData) * 100
		
		tempParam[,1] <- colnames(tempData)
		tempParam[,2] <- grpNames[i]
		tempParam[,3] <- apply(tempData, 2, mean)
		tempParam[,4] <- apply(tempData, 2, max)
		tempParam[,5] <- apply(tempData, 2, min)
		
		allParams <- rbind(allParams, tempParam)
	}
	colnames(allParams) <- c("Taxa", "Group", "Avg", "Max", "Min")
	
	# Factor the data so it stays in the right order
	allParams$Group <- factor(allParams$Group, levels=grpNames)
	allParams$Taxa <- factor(allParams$Taxa, levels=unique(colnames(group.data[[1]])))
	
	print(ggplot2::ggplot(allParams, aes_string(y="Avg", x="Taxa", colour="Group")) +
					geom_point() + theme(legend.position = "top") +
					geom_errorbar(aes_string(ymax="Max", ymin="Min")) +
					labs(title=main, y="Percent Abundance", x="") +
					theme(axis.text.x=element_text(hjust=1, angle=45, size=6)
					)
	)
}

Plot.MDS <- function(group.data, main="Group MDS", retCords=FALSE){
	if(missing(group.data))
		stop("group.data is missing.")
	
	numGroups <- length(group.data)
	
	# Check every data set has the same number of taxa
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa))
		stop("Every data set must have matching taxa.")
	
	# Make sure we have group names
	if(is.null(names(group.data))){
		grpNames <- paste("Data Set", 1:numGroups)
	}else{
		grpNames <- names(group.data)
	}
	
	# Merge all the data sets together
	mData <- do.call("rbind", group.data)
	
	# Get their mds location 
	loc <- getBC(mData)
	
	# Set color
	availCols <- rainbow(numGroups)
	cols <- NULL
	for(i in 1:numGroups)
		cols <- c(cols, rep(availCols[i], nrow(group.data[[i]])))
	
	# Plot MDS
	graphics::plot(loc, pch=16, ylab="MDS 2", xlab="MDS 1", col=cols)
	graphics::legend("topright", legend=grpNames, pch=16, col=availCols)
	
	if(retCords)
		return(loc)
}


### ~~~~~~~~~~~~~~~~~~~~~
### Filter functions
### ~~~~~~~~~~~~~~~~~~~~~
formatDataSets <- function(group.data, data=NULL){
	if(missing(group.data) && is.null(data))
		stop("group.data missing.")
	
	# Check is data is still being used
	if(!is.null(data) && missing(group.data))
		group.data <- data
	
	numGroups <- length(group.data)
	if(numGroups < 2)
		stop("At least 2 data sets are required.")
	
	newData <- NULL
	for(i in 1:length(group.data)){		
		newData <- merge(newData, t(group.data[[i]]), by=0, all=TRUE)
		rownames(newData) <- newData[,1]
		newData <- newData[,-1]
	}
	newData[is.na(newData)] <- 0
	newData <- t(newData)
	newData <- newData[,colSums(newData) != 0]
	newData <- newData[rowSums(newData) != 0,]
	newData <- newData[,order(colSums(newData), decreasing=TRUE)]
	
	retData <- vector("list", numGroups)
	base <- 0
	for(i in 1:numGroups){
		retData[[i]] <- newData[(base+1):(nrow(group.data[[i]])+ base),]
		base <- base + nrow(group.data[[i]])
	}
	
	names(retData) <- names(group.data)
	return(retData)
}

Data.filter <- function(data, order.type="data", minReads=0, numTaxa=NULL, perTaxa=NULL, K=NULL, reads.crit=NULL){
	if(missing(data))
		stop("data is missing.")
	if(tolower(order.type) != "data" && tolower(order.type) != "sample")
		stop(sprintf("'%s' not recognized, order.type must be 'data' or 'sample'", as.character(order.type)))
	
	# Check if K is still being used
	if(is.null(numTaxa) && !is.null(K))
		numTaxa <- K
	# Check if reads.crit is still being used
	if(!is.null(reads.crit))
		minReads <- reads.crit
	
	# Check if numTaxa or perTaxa is being used
	if(!is.null(numTaxa) && !is.null(perTaxa))
		stop("numTaxa and perTaxa cannot be used at the same time")
	if(!is.null(numTaxa)){
		if(numTaxa >= ncol(data) || numTaxa <= 0)
			stop(sprintf("numTaxa must be between 0 and %i.", ncol(data)-1))
	}
	if(!is.null(perTaxa)){
		if(perTaxa >= 1 || perTaxa <= 0)
			stop("perTaxa must be between 0 and 1.")
	}
	
	taxaNames <- colnames(data)
	
	# Drop all subjects that don't have enough reads
	data <- data[rowSums(data)>minReads,, drop=FALSE]
	if(nrow(data) < 2)
		stop("minReads is too large and is excluding too many samples.  Please try lowering its value.")
	
	# Drop all 0 taxa
	data <- data[,colSums(data)>0, drop=FALSE]
	
	# Order the data based on order.type
	if(tolower(order.type) == "sample"){
		data <- t(apply(data, 1, function(x){x[order(x, decreasing=TRUE)]}))
	}else{
		data <- data[,order(colSums(data), decreasing=TRUE)]
	}
	
	# Use a percentage based approach to find the number of taxa to collapse
	if(!is.null(perTaxa)){
		perNumReadsTaxa <- colSums(data)/sum(data)
		cumSumReads <- cumsum(perNumReadsTaxa)
		numTaxa <- max(which(cumSumReads <= perTaxa))
	}
	
	# Pull out the taxa we want to collapse
	otherData <- data[,-c(1:numTaxa), drop=FALSE]
	
	# Put the data back together and relabel
	retData <- cbind(data[,1:numTaxa], Other=rowSums(otherData))
	
	return(retData)
}


### ~~~~~~~~~~~~~~~~~~~~~
### MC functions
### ~~~~~~~~~~~~~~~~~~~~~
MC.ZT.statistics <- function(Nrs, numMC=10, fit, type="ha", siglev=0.05, MC=NULL) {
	if(missing(Nrs) || missing(fit))
		stop("Nrs and/or fit missing.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	
	# Check if someone is still using MC
	if(!is.null(MC))
		numMC <- MC
	
	# Get all the ZT values
	ZTstatMatrix <- matrix(0, numMC, 2)
	for(i in 1:numMC)
		ZTstatMatrix[i,] <- ZT.statistics.Hnull.Ha(Nrs, fit, type)
	
	# Pull out z and t and remove NAs
	z <- ZTstatMatrix[,1]
	z <- z[!is.na(z)]
	t <- ZTstatMatrix[,2]
	t <- t[!is.na(t)]
	
	# Get a reference value from the real data
	qAlpha <- qchisq(p=(1-siglev), df=length(fit$pi)-1, ncp=0, lower.tail=TRUE)
	
	# Calculate our pvalues for z and t
	zpval <- sum(z > qAlpha)/length(z)
	tpval <- sum(t > qAlpha)/length(t)
	
	return(cbind(zpval, tpval))
}

MC.Xsc.statistics <- function(Nrs, numMC=10, fit, pi0=NULL, type="ha", siglev=0.05, MC=NULL) {
	if(missing(Nrs) || missing(fit))
		stop("Nrs and/or fit missing.")
	if(is.null(pi0) && tolower(type) == "ha")
		stop("pi0 cannot be null with type 'ha'.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	
	# Check if someone is still using MC
	if(!is.null(MC))
		numMC <- MC
	
	# Get all the XSC values
	XscStatVector <- rep(0, numMC)
	for(i in 1:numMC)
		XscStatVector[i] <- Xsc.statistics.Hnull.Ha(Nrs, fit, type, pi0)
	
	# Remove NAs
	XscStatVector <- XscStatVector[!is.na(XscStatVector)]
	
	# Get a reference value from the real data
	qAlpha <- qchisq(p=(1-siglev), df=length(fit$pi)-1, ncp=0, lower.tail=TRUE)
	
	# Calculate pvalues
	pval <- sum(XscStatVector > qAlpha)/length(XscStatVector)
	
	return(pval)
}

MC.Xmc.statistics <- function(group.Nrs, numMC=10, pi0, group.pi, group.theta, type="ha", siglev=0.05, MC=NULL, Nrs=NULL) {
	if(missing(group.theta) || missing(pi0) || (missing(group.Nrs) && is.null(Nrs)))
		stop("group.Nrs, pi0 and/or group.theta missing.")
	if(missing(group.pi) && tolower(type) == "ha")
		stop("group.pi missing.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	
	# Check if someone is still using MC
	if(!is.null(MC))
		numMC <- MC
	
	# Check if someone is still using Nrs
	if(!is.null(Nrs))
		group.Nrs <- Nrs
	
	numGroups <- length(group.Nrs)
	numTaxa <- length(pi0)
	
	# If the type is ha this will change in the for loop
	tempPi <- pi0
	
	# Create the parameters for every group
	groupParameter <- vector("list", numGroups)
	for (i in 1:numGroups){
		if(tolower(type) == "ha")
			tempPi <- group.pi[i,]
		groupParameter[[i]] <- list(pi=tempPi, theta=group.theta[i], nrs=group.Nrs[[i]])
	}
	
	# Get all the Xmc values
	XmcStatVector <- rep(0, numMC)
	for(i in 1:numMC)
		XmcStatVector[i] <- Xmc.statistics.Hnull.Ha(groupParameter, pi0)
	
	# Get a reference value from the real data
	qAlpha <- qchisq(p=(1-siglev), df=length(group.theta)*(numTaxa-1), ncp=0, lower.tail=TRUE)
	
	# Calculate pvalues
	pval <- sum(XmcStatVector > qAlpha)/length(XmcStatVector)
	
	return(pval)
}

MC.Xmcupo.statistics <- function(group.Nrs, numMC=10, pi0, group.pi, group.theta, type="ha", siglev=0.05, MC=NULL, Nrs=NULL) {
	if(missing(group.theta) || (missing(group.Nrs) && is.null(Nrs)))
		stop("group.Nrs and/or group.theta missing.")
	if(missing(group.pi) && tolower(type) == "ha")
		stop("group.pi missing.")
	if(missing(pi0) && tolower(type) == "hnull")
		stop("pi0 missing.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	
	# Check if someone is still using MC
	if(!is.null(MC))
		numMC <- MC
	
	# Check if someone is still using Nrs
	if(!is.null(Nrs))
		group.Nrs <- Nrs
	
	numGroups <- length(group.Nrs)
	
	# Create the parameters for every group
	groupParameter <- vector("list", numGroups)
	for (i in 1:numGroups){
		if(tolower(type) == "ha"){
			numTaxa <- ncol(group.pi)
			tempPi <- group.pi[i,]
		}else{
			numTaxa <- length(pi0)
			tempPi <- pi0
		}
		groupParameter[[i]] <- list(pi=tempPi, theta=group.theta[i], nrs=group.Nrs[[i]])
	}
	
	# Get all the Xmcupo values
	XmcupoStatVector <- rep(0, numMC)
	for(i in 1:numMC)
		XmcupoStatVector[i] <- Xmcupo.statistics.Hnull.Ha(groupParameter)
	
	# Get a reference value from the real data
	qAlpha <- qchisq(p=(1-siglev), df=length(group.theta)*(numTaxa-1), ncp=0, lower.tail=TRUE)
	
	# Calculate pvalues
	pval <- sum(XmcupoStatVector > qAlpha)/length(XmcupoStatVector)
	
	return(pval)
}

MC.Xdc.statistics <- function(group.Nrs, numMC=10, alphap, type="ha", siglev=0.05, est="mom", MC=NULL, Nrs=NULL) {
	if(missing(alphap) || (missing(group.Nrs) && is.null(Nrs)))
		stop("group.Nrs and/or alphap  missing.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	if(tolower(est) != "mom" && tolower(est) != "mle")
		stop(sprintf("Est '%s' not found. Est must be 'mle' or 'mom'.", as.character(est)))
	
	# Check if someone is still using MC
	if(!is.null(MC))
		numMC <- MC
	
	# Check if someone is still using Nrs
	if(!is.null(Nrs))
		group.Nrs <- Nrs
	
	numGroups <- length(group.Nrs)	
	
	if(tolower(type) == "hnull"){
		numTaxa <- length(alphap)
	}else{
		numTaxa <- ncol(alphap)
	}
	
	# Get all the Xdc values
	XdcStatVector <- rep(0, numMC)
	for(i in 1:numMC)
		XdcStatVector[i] <- Xdc.statistics.Hnull.Ha(alphap, group.Nrs, type, est)
	
	# Get a reference value from the real data
	qAlpha <- qchisq(p=(1-siglev), df=(numGroups-1)*numTaxa, ncp=0, lower.tail=TRUE)
	
	# Calculate pvalues
	pval <- sum(XdcStatVector > qAlpha)/length(XdcStatVector)
	
	return(pval)
}

MC.Xoc.statistics <- function(group.Nrs, numMC=10, group.alphap, type="ha", siglev=0.05, MC=NULL, Nrs=NULL) {
	if(missing(group.alphap) || (missing(group.Nrs) && is.null(Nrs)))
		stop("group.Nrs and/or group.alphap missing.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	
	# Check if someone is still using MC
	if(!is.null(MC))
		numMC <- MC
	
	# Check if someone is still using Nrs
	if(!is.null(Nrs))
		group.Nrs <- Nrs
	
	numGroups <- length(group.Nrs)	
	
	# Get all the Xoc values
	XocStatVector <- rep(0, numMC)
	for(i in 1:numMC)
		XocStatVector[i] <- Xoc.statistics.Hnull.Ha(group.Nrs, group.alphap, type)
	
	# Get a reference value from the real data
	qAlpha <- qchisq(p=(1-siglev), df=(numGroups-1), ncp=0, lower.tail=TRUE)
	
	# Calculate pvalues
	pval <- sum(XocStatVector > qAlpha)/length(XocStatVector)
	
	return(pval)
}



### ~~~~~~~~~~~~~~~~~~~~~
### Sample functions
### ~~~~~~~~~~~~~~~~~~~~~
Xsc.onesample <- function(data, pi0){
	if(missing(data) || missing(pi0))
		stop("data and/or pi0 missing.")
	
	numReadsSubs <- rowSums(data)
	numTaxa	<- length(pi0)
	
	# Check the data set has the same number of taxa
	numTaxa	<- length(pi0)
	if(ncol(data) != numTaxa)
		stop("Every data set must have the same length as pi0")
	
	# Get parameters
	fit.MoM <- DM.MoM(data)
	
	# Get Xsc and calculate pvalue
	Xsc <- Xsc.statistics(fit.MoM$pi, fit.MoM$theta, numReadsSubs, pi0)
	pval <- 1-pchisq(q=Xsc, df=numTaxa-1, ncp=0, lower.tail=TRUE)
	
	RAD.mean.test <- list("Xsc statistics"=Xsc, "p value"=pval)
	
	return(RAD.mean.test)			
}

Xmc.sevsample <- function(group.data, pi0){	
	if(missing(group.data) || missing(pi0))
		stop("group.data and/or pi0 missing.")
	
	# Check every data set has the same number of taxa
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- length(pi0)
	if(any(taxaCounts != numTaxa))
		stop("Every data set must have matching taxa, including pi0")
	
	numGroups <- length(group.data)
	
	# Get the parameters for every group
	groupParameter <- lapply(group.data, function(x){
				# Calc pi, theta and the number of reads
				numReadsSubs <- rowSums(x)
				pi.MoM <- colSums(x)/sum(x)
				theta.MoM <- weirMoM(x, pi.MoM)$theta
				
				return(list(pi=pi.MoM, theta=theta.MoM, nrs=numReadsSubs))
			})
	
	# Get Xmc and calculate pvalue
	Xmc <- Xmc.statistics(groupParameter, pi0)
	pval <- 1-pchisq(q=Xmc, df=numGroups*(numTaxa-1), ncp=0, lower.tail=TRUE)
	
	sevRAD.mean.test <- list("Xmc statistics"=Xmc, "p value"=pval)
	
	return(sevRAD.mean.test)					
}

Xmcupo.sevsample <- function(group.data){
	if(missing(group.data))
		stop("group.data is missing.")
	
	# Check every data set has the same number of taxa
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa))
		stop("Every data set must have matching taxa")
	
	numGroups <- length(group.data)
	
	# Get the parameters for every group
	groupParameter <- lapply(group.data, function(x){
				# Calc pi, theta and the number of reads
				numReadsSubs <- rowSums(x)
				pi.MoM <- colSums(x)/sum(x)
				theta.MoM <- weirMoM(x, pi.MoM)$theta
				
				return(list(pi=pi.MoM, theta=theta.MoM, nrs=numReadsSubs))
			})
	
	# Get Xmcupo and calculate pvalue
	Xmcupo <- Xmcupo.statistics(groupParameter)
	pval <- 1-pchisq(q=Xmcupo, df=(numGroups-1)*(numTaxa-1), ncp=0, lower.tail=TRUE)
	
	ret <- list("Xmcupo statistics"=Xmcupo, "p value"=pval)
	
	return(ret)	
}

# reviewed 9/29/16
Xoc.sevsample <- function(group.data, epsilon=10^(-4)){
	if(missing(group.data))
		stop("group.data missing.")
	
	# Check every data set has the same number of taxa
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa))
		stop("Every data set must have matching taxa")
	
	numGroups <- length(group.data)
	
	# Get Xoc and calculate pvalue
	Xoc <- Xoc.statistics(group.data, epsilon)
	pval <- 1-pchisq(q=Xoc, df=numGroups-1, ncp=0, lower.tail=TRUE)	
	
	sev.overd.test <- list("Xoc statistics"=Xoc, "p value"=pval)
	
	return(sev.overd.test)			
}

# reviewed 9/29/16
Xdc.sevsample <- function(group.data, epsilon=10^(-4), est="mom"){
	if(missing(group.data))
		stop("group.data missing.")
	if(tolower(est) != "mle" && tolower(est) != "mom")
		stop(sprintf("Est '%s' not found. Est must be 'mle' or 'mom'.", as.character(est)))
	
	
	# Check every data set has the same number of taxa
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa))
		stop("Every data set must have matching taxa")
	
	numGroups <- length(group.data)
	
	# Get Xdc and calculate pvalue
	if(tolower(est) == "mle"){
		Xdc <- Xdc.statistics(group.data, epsilon)
	}else{
		Xdc <- Xdc.statistics.MoM(group.data)
	}
	pval <- 1-pchisq(q=Xdc, df=(numGroups-1)*numTaxa, ncp=0, lower.tail=TRUE)		
	
	xdc.sevsamp.test <- list("Xdc statistics"=Xdc, "p value"=pval)
	
	return(xdc.sevsamp.test)	
}





### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Internal
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### ~~~~~~~~~~~~~~~~~~~~~
### Other functions
### ~~~~~~~~~~~~~~~~~~~~~
# reviewed 9/30/16
loglikDM <- function(data, alphap){	
	data <- data[,colSums(data)!=0, drop=FALSE]
	alphap <- alphap[alphap!=0]
	
	ll <- sum(lgamma(rowSums(data)+1) + lgamma(sum(alphap)) - lgamma(sum(alphap)+rowSums(data))) + 
			sum(rowSums(lgamma(sweep(data, 2, alphap, "+")) - lgamma(data+1) - lgamma(t(replicate(nrow(data), alphap)))))
	
	return(ll)
}

# reviewed 9/30/16
weirMoM <- function(data, MoM, se=FALSE){
	numTaxa <- ncol(data)
	numSamp <- nrow(data)
	rowSumsData <- rowSums(data) + 0.000001
	colSumsData <- colSums(data)
	
	if(numSamp == 1)
		return(list(theta=0, se=0))
	
	MSP <- (numSamp-1)^(-1) * sum(rowSums((data/rowSumsData - matrix(rep(MoM, numSamp), numSamp, numTaxa, byrow=TRUE))^2) * rowSumsData)
	MSG <- (sum(colSumsData)-numSamp)^(-1) * sum(rowSums(data/rowSumsData * (1-data/rowSumsData)) * rowSumsData)
	nc <- 1/(numSamp-1) * (sum(rowSumsData)-sum(rowSumsData^2)/sum(rowSumsData))
	MoM.wh <- (MSP-MSG)/(MSP+(nc-1)*MSG)
	
	std.er <- NULL
	if(se)
		std.er <- sqrt(2 * (1-MoM.wh)^2/(numSamp-1) * ((1+(nc-1) * MoM.wh)/nc)^2)
	
	return(list(theta=MoM.wh, se=std.er))
}

kullbackLeiber <- function(data, plot=TRUE, parallel=FALSE, cores=3){
	warning("This function has been spellchecked.  Please use 'Kullback.Leibler' instead.")
	kl <- Kullback.Leibler(data, plot, parallel, cores)
	return(kl)
}

getBC <- function(data){
	dataPer <- data/rowSums(data)
	bcDist <- vegan::vegdist(dataPer, method="bray")
	nonMetricMDS <- MASS::isoMDS(bcDist, trace=FALSE)
	mdsPoints <- vegan::postMDS(nonMetricMDS$points, bcDist)
	mds <- vegan::scores(mdsPoints)
	
	return(mds[,1:2])
}



### ~~~~~~~~~~~~~~~~~~~~~
### Stat functions
### ~~~~~~~~~~~~~~~~~~~~~
Xmcupo.statistics <- function(groupParameter){	
	numGroups <- length(groupParameter)	
	numTaxa <- length(groupParameter[[1]]$pi)
	
	# Pull out pi and calculate xsc
	pis <- matrix(0, numTaxa, numGroups)
	xscs <- rep(0, numGroups)
	for(i in 1:numGroups){
		theta <- groupParameter[[i]]$theta
		numReads <- groupParameter[[i]]$nrs
		totalReads <- sum(numReads)
		
		# Calculate the Xsc for each group
		xscs[i] <- (theta * (sum(numReads^2)-totalReads) + totalReads) / totalReads^2
		pis[,i] <- groupParameter[[i]]$pi
	}
	
	# Remove any 0 taxa
	pis <- pis[rowSums(pis)!=0,]
	
	# Calculate pi0
	pi0 <- rowSums(pis/xscs)/sum(1/xscs)
	
	# Calculate Xmcupo
	Xmcupo <- sum(colSums((pis-pi0)^2/pi0)/xscs)		
	
	return(Xmcupo)
}

Z.statistics <- function(data){
	numTaxa <- ncol(data)
	numReadsTaxa <- colSums(data)
	numReadsSubs <- rowSums(data)
	totalReads <- sum(data)
	
	taxaSqSum <- sum(apply(data, 2, function(x){sum((x-1)*x)})/numReadsTaxa) 
	subSqSum <- sum(numReadsSubs*(numReadsSubs-1))
	
	denom <- sqrt(2*(numTaxa-1) * subSqSum)
	
	Zs <- (totalReads*taxaSqSum-subSqSum)/denom
	
	return(Zs)
}

T.statistics <- function(data){	
	numReadsTaxa <- colSums(data)
	numReadsSubs <- rowSums(data)
	totalReads <- sum(data)
	
	Ts <- sum(colSums((data - (numReadsSubs%*%t(numReadsTaxa))/totalReads)^2) / numReadsTaxa)
	
	return(Ts)
}

Xmc.statistics <- function(groupParameter, pi0){
	numGroups <- length(groupParameter)	
	numTaxa <- length(pi0)
	
	xsc <- rep(0, numGroups)
	for(i in 1:numGroups){
		pi <- groupParameter[[i]]$pi
		theta <- groupParameter[[i]]$theta
		numReads <- groupParameter[[i]]$nrs
		
		# Get Xsc values
		xsc[i] <- Xsc.statistics(pi, theta, numReads, pi0)
	}
	xmc <- sum(xsc)
	
	return(xmc)
}

Xsc.statistics <- function(pi1, theta, numReads, pi0){
	totalReads <- sum(numReads)
	
	# Get Xsc value
	tempVal <- ((theta*(sum(numReads^2)-totalReads) + totalReads) / totalReads^2) * (diag(pi0)-pi0 %*% t(pi0))
	xsc <- as.vector(t(pi1-pi0) %*% MASS::ginv(tempVal) %*% (pi1-pi0))	
	
	return(xsc)
}

Xoc.statistics <- function(group.data, epsilon=10^(-4)){
	numGroups <- length(group.data)
	
	# Get the fit for every data set
	thetas <- rep(0, numGroups)
	logliks <- rep(0, numGroups)
	pis <- vector("list", numGroups)
	for(i in 1:numGroups){
		fit <- dirmult::dirmult(group.data[[i]], initscalar=DM.MoM(group.data[[i]])$theta, epsilon=epsilon, trace=FALSE)
		
		thetas[i] <- fit$theta
		logliks[i] <- fit$loglik
		pis[[i]] <- fit$pi
	}
	
	# Get the fit assuming equal thetas
	equalFit <- dirmult::equalTheta(group.data, mean(thetas), epsilon, FALSE, pis)
	
	# Calculate the xoc
	xoc <- as.vector(-2*(equalFit$loglik-sum(logliks)))
	
	return(xoc)
}

Xdc.statistics <- function(group.data, epsilon=10^(-4)){ 
	# Get the loglik from the fit from every data set
	logliks <- sapply(group.data, function(x, epsilon){
				dirmult::dirmult(x, initscalar=DM.MoM(x)$theta, epsilon=epsilon, trace=FALSE)$loglik
			}, epsilon=epsilon)
	
	# Get the fit assuming all in the same group
	groupDataC <- do.call(rbind, group.data)
	groupFit <- dirmult::dirmult(groupDataC, initscalar=DM.MoM(groupDataC)$theta, epsilon=epsilon, trace=FALSE)	
	
	# Calculate the xdc
	xdc <- -2*(groupFit$loglik-sum(logliks))
	
	return(xdc)
}

# reviewed 9/29/16
Xdc.statistics.MoM <- function(group.data){	
	# Get the loglik from the fit from every data set
	logliks <- sapply(group.data, function(x){DM.MoM(x)$loglik})
	
	# Get the fit assuming all in the same group
	groupDataC <- do.call(rbind, group.data)
	groupFit <- DM.MoM(groupDataC)	
	
	# Calculate the xdc
	xdc <- -2*(groupFit$loglik-sum(logliks))
	
	return(xdc)
}



### ~~~~~~~~~~~~~~~~~~~~~
### Hnull / Ha functions
### ~~~~~~~~~~~~~~~~~~~~~
Xmcupo.statistics.Hnull.Ha <- function(groupParameter){		
	numGroups <- length(groupParameter)
	numTaxa <- length(groupParameter[[1]]$pi)
	
	genGroupParameter <- vector("list", numGroups)
	for(i in 1:numGroups){
		pi <- groupParameter[[i]]$pi
		theta <- groupParameter[[i]]$theta
		numReads <- groupParameter[[i]]$nrs
		
		# Generate a new set of data
		genData <- Dirichlet.multinomial(numReads, pi*(1-theta)/theta)
		genTotalReads <- sum(genData)
		genPi <- colSums(genData)/genTotalReads
		
		# Replace any 0 pi values with a small number
		# This will subtract that value from the other data so a total value of 1 is maintained
		if(any(genPi==0)){
			numZero <- sum(genPi==0)
			numNonZero <- numTaxa - numZero
			genPi[which(genPi!=0)] <- genPi[which(genPi!=0)] - numZero/(numNonZero*2*(genTotalReads+1))
			genPi[which(genPi==0)] <- 1/(2*(genTotalReads+1))
		}
		
		genTheta <- weirMoM(genData, genPi)$theta
		genGroupParameter[[i]] <- list(pi=genPi, theta=genTheta, nrs=numReads)			
	}
	
	# Get the Xmcupo stats for the generated data
	Xmcupo <- Xmcupo.statistics(genGroupParameter)		
	
	return(Xmcupo)
}

ZT.statistics.Hnull.Ha <- function(Nrs, fit, type){
	if(tolower(type) == "hnull"){
		genData <- Multinomial(Nrs, fit$pi)
	}else{
		genData <- Dirichlet.multinomial(Nrs, fit$gamma)
	}
	
	ZT <- c(Z.statistics(genData), T.statistics(genData))
	
	return(ZT)
}

Xmc.statistics.Hnull.Ha <- function(groupParameter, pi0){	
	numGroups <- length(groupParameter)	
	numTaxa <- length(pi0)
	
	genGroupParameter <- vector("list", numGroups)
	for(i in 1:numGroups){
		pi <- groupParameter[[i]]$pi
		theta <- groupParameter[[i]]$theta
		numReads <- groupParameter[[i]]$nrs
		
		# Generate a new set of data
		genData <- Dirichlet.multinomial(numReads, pi*(1-theta)/theta)				
		genPi <- colSums(genData)/sum(genData)
		genTheta <- weirMoM(genData, genPi)$theta
		
		genGroupParameter[[i]] <- list(pi=genPi, theta=genTheta, nrs=numReads)		
	}
	
	# Get the Xmc stats for the generated data
	Xmc <- Xmc.statistics(genGroupParameter, pi0)	
	
	return(Xmc)
}

Xsc.statistics.Hnull.Ha <- function(Nrs, fit, type, pi0){
	# Generate a new set of data
	genData <- Dirichlet.multinomial(Nrs, fit$gamma)
	fit.gen <- DM.MoM(genData)			
	
	tempPi <- fit$pi
	if(tolower(type) == "ha")
		tempPi <- pi0
	
	# Calculate Xsc stat
	xsc <- Xsc.statistics(fit.gen$pi, fit.gen$theta, Nrs, tempPi)
	
	return(xsc)
}

Xoc.statistics.Hnull.Ha <- function(group.Nrs, group.alphap, type){
	numGroups <- length(group.Nrs)
	tempShape <- group.alphap
	
	# Generate a new set of data
	genGroupData <- vector("list", numGroups)
	for(i in 1:numGroups){
		if(tolower(type) == "ha")
			tempShape <- group.alphap[i,]
		
		genGroupData[[i]] <- Dirichlet.multinomial(group.Nrs[[i]], tempShape)
	}
	
	# Get the xoc stats for the generated data
	xoc <- Xoc.statistics(genGroupData)	
	
	return(xoc)
}

Xdc.statistics.Hnull.Ha <- function(alphap, group.Nrs, type, est){
	numGroups <- length(group.Nrs)
	tempShape <- alphap
	
	# Generate a new set of data
	genGroupData <- vector("list", numGroups)
	for(i in 1:numGroups){
		if(tolower(type) == "ha")
			tempShape <- alphap[i,]
		
		genGroupData[[i]] <- Dirichlet.multinomial(group.Nrs[[i]], tempShape)
	}
	
	# Get the xdc stats for the generated data
	if(tolower(est) == "mle"){
		xdc <- Xdc.statistics(genGroupData)
	}else{
		xdc <- Xdc.statistics.MoM(genGroupData)
	}
	
	return(xdc)
}





### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### OLD / Removable???
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Xoc.statistics.MoM <- function(group.data){
#	numGroups <- length(group.data)
#	
#	# Get the theta from the fit assuming all in the same group
#	groupDataC <- do.call(rbind, group.data)
#	groupFit <- DM.MoM(groupDataC)
#	groupTheta <- groupFit$theta
#	
#	# Get the loglik from the fit from every data set
#	# Also get the loglik if thetas were equal
#	logliks <- rep(0, numGroups)
#	equalThetaLogliks <- rep(0, numGroups)
#	for(i in 1:numGroups){
#		fit <- dirmult::dirmult(group.data[[i]], initscalar=DM.MoM(group.data[[i]])$theta, epsilon=epsilon, trace=FALSE)
#		logliks[i] <- fit$loglik
#		
#		equalThetaLogliks[i] <- loglikDM(group.data[[i]], fit$pi*(1-groupTheta)/groupTheta)	
#	}
#	
#	# Calculate the xoc
#	xoc <- -2*(sum(equalThetaLogliks)-sum(logliks))	
#	
#	return(xoc)
#}
#
#LLDM <- function(data, alpha){
#	numReadsSubs <- rowSums(data)
#	
#	# Calculate the lprob for every row
#	lprob <- rep(0, nrow(data))
#	for(i in 1:nrow(data)){
#		lprob[i] <- log(numReadsSubs[i]) + lbeta(sum(alpha), numReadsSubs[i]) -
#				sum(log(data[i,])) - sum(lbeta(alpha , as.numeric(data[i,])))
#	}
#	
#	# Sum the lprob to get the loglik
#	ll <- sum(lprob)
#	
#	return(ll)
#}
#
#Plot.MLE.old <- function(group.data, main="Group MLE Plot", returnData=FALSE){
#	if(missing(group.data))
#		stop("group.data is missing.")
#	
#	# Check the number of groups
#	numGroups <- length(group.data)
#	
#	# Check every data set has the same number of taxa
#	taxaCounts <- sapply(group.data, ncol)
#	numTaxa	<- taxaCounts[1]
#	if(any(taxaCounts != numTaxa))
#		stop("Every data set must have matching taxa.")
#	
#	# Make sure we have group names
#	if(is.null(names(group.data))){
#		grpNames <- paste("Data Set", 1:numGroups)
#	}else{
#		grpNames <- names(group.data)
#	}
#	
#	# Remove periods in taxa names and check for taxa with 0 column sums (add 1 to subject 1)
#	group.data <- lapply(group.data, function(x){
#				colnames(x) <- gsub(pattern="[.]", replacement="_", x=colnames(x))
#				badTaxa <- which(colSums(x) == 0)
#				if(length(badTaxa) != 0)
#					x[1, badTaxa] <- x[1, badTaxa] + 1
#				return(x)
#			})
#	
#	# Get the mle and mom for every group
#	fit <- lapply(group.data, function(x){dirmult(x, trace=FALSE)})
#	params <- mapply(function(x, y){dirmult.summary(x, y)}, x=group.data, y=fit, SIMPLIFY=FALSE)
#	
#	# Pull out just the mle
#	params <- lapply(params, function(x){x[,c("MLE", "se.MLE")]})
#	params <- lapply(params, function(x){x[-nrow(x),]})
#	
#	# Add taxa columns
#	params <- lapply(params, function(x){cbind(x, Taxa=rownames(x), stringsAsFactors=FALSE)})
#	
#	# Collapse into a single dataframe
#	allParams <- do.call(rbind, params)
#	colnames(allParams) <- c("pi", "se", "Taxa")
#	allParams$Group	<- rep(grpNames, each=numTaxa)
#	
#	#  Make Taxa an ordered factor
#	allParams$Group <- factor(allParams$Group, levels=unique(grpNames))
#	allParams$Taxa <- factor(allParams$Taxa, levels=unique(colnames(group.data[[1]])))
#	
#	print(ggplot(allParams, aes(y=pi, x=Taxa, colour=Group)) +
#					geom_point() + theme(legend.position = "top") +
#					geom_errorbar(aes(ymax = pi + 1.96*se, ymin = pi - 1.96*se)) +
#					labs(title=main, y=expression(paste("MLE ", pi)), x="") +
#					theme(axis.text.x=element_text(hjust=1, angle=45, size=6)
#					)
#	)
#	
#	if(returnData)
#		return(allParams)
#}











