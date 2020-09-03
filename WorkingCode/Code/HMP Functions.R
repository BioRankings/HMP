

library(dirmult)	# example code functions and dirmult, equalTheta functions
library(ggplot2)	# pretty plots for pi
library(doParallel) # parallelizing 
library(MASS)		# ginv function for Xsc.statistics
library(vegan) 		# ga distances
library(gplots)		# heatmap plots
library(rpart)		# base rpart
library(rpart.plot)	# rpart plotting
library(lattice)	# Repeated measures plotting


### Define a global environment to use with rpart
hmp.pkg.env <- new.env(parent=emptyenv())
hmp.pkg.env$EVAL_COUNT_RPART <- 1


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### External
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### ~~~~~~~~~~~~~~~~~~~~~
### Generation functions
### ~~~~~~~~~~~~~~~~~~~~~
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

Kullback.Leibler <- function(group.data, plot=TRUE, main="Kullback Leibler Divergences", parallel=FALSE, cores=3){
	if(missing(group.data))
		stop("data missing.")
	
	# Check the number of groups
	numGrps <- length(group.data)
	if(numGrps < 2)
		stop("At least 2 data sets are required.")
	
	# Make sure we have the same columns
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa)){
		warning("Group columns do not match, running formatDataSets.")
		group.data <- formatDataSets(group.data)
	}
	
	# Make sure we have group names
	if(is.null(names(group.data))){
		grpNames <- paste("Data Set", 1:numGrps)
	}else{
		grpNames <- names(group.data)
	}
	
	# Add 1 so we don't ever get an all 0 comparison
	group.data <- lapply(group.data, function(x) x+1)  
	
	# Run dirmult on every group
	if(parallel){
		cl <- parallel::makeCluster(min(cores, numGrps)) 
		doParallel::registerDoParallel(cl)
		
		tryCatch({
					results <- foreach::foreach(i=1:numGrps, .combine=list, .multicombine=TRUE, .inorder=TRUE, .packages=c("dirmult")) %dopar%{
						param <- DM.MoM(group.data[[i]])
						return(param)
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		results <- vector("list", numGrps)
		for(i in 1:numGrps)	
			results[[i]] <- DM.MoM(group.data[[i]])
	}
	
	# Get alpha for every group
	alpha <- lapply(results, function(x) x$gamma)
	names(alpha) <- grpNames
	
	# Get LL given alpha
	LL.vals <- sapply(results, function(x) x$loglik)
	
	# Get LL for every group using another alpha
	KLmat <- matrix(0, numGrps, numGrps)		
	for(i in 1:numGrps){
		for(j in i:numGrps){
			if(i == j)
				next
			KLval1 <- LL.vals[i] - loglikDM(group.data[[i]], alpha[[j]])
			KLval2 <- LL.vals[j] - loglikDM(group.data[[j]], alpha[[i]])
			KLmat[i, j] <- KLval1 + KLval2
			KLmat[j, i] <- KLval1 + KLval2
		}
	}
	colnames(KLmat) <- grpNames
	rownames(KLmat) <- grpNames
	
	if(plot){
		gplots::heatmap.2(as.matrix(KLmat), dendrogram="both", Rowv=TRUE, Colv=TRUE, 
				trace="none", symm=TRUE, margins=c(12, 9), density.info="none",
				main=main
		)
	}
	
	return(KLmat)
}

Xmcupo.effectsize <- function(group.data){
	if(missing(group.data))
		stop("group.data missing.")
	
	# Make sure we have the same columns
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa)){
		warning("Group columns do not match, running formatDataSets.")
		group.data <- formatDataSets(group.data)
		numTaxa <- ncol(group.data[[1]])
	}
	
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

Est.PI <- function(group.data, conf=.95){
	if(missing(group.data))
		stop("group.data is missing.")
	
	# Check the number of groups
	numGroups <- length(group.data)
	
	# Make sure we have the same columns
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa)){
		warning("Group columns do not match, running formatDataSets.")
		group.data <- formatDataSets(group.data)
	}
	
	# Make sure we have group names
	if(is.null(names(group.data))){
		grpNames <- paste("Data Set", 1:numGroups)
	}else{
		grpNames <- names(group.data)
	}
	
	# Calculate the pi and error bars for each group
	allParamsMLE <- data.frame(matrix(0, 0, 6))
	allParamsMOM <- data.frame(matrix(0, 0, 6))
	thetaMLE <- data.frame(matrix(0, numGroups, 3))
	thetaMOM <- data.frame(matrix(0, numGroups, 3))
	for(i in 1:numGroups){
		tempData <- group.data[[i]]
		
		# Check the data has samples
		numSub <- nrow(tempData)
		if(numSub < 1)
			stop("At least one data set in group.data is empty")
		
		tempParam1 <- data.frame(matrix(0, ncol(tempData), 6))
		tempParam2 <- data.frame(matrix(0, ncol(tempData), 6))
		
		tempParam2[,2] <- grpNames[i]
		tempParam1[,2] <- grpNames[i]
		
		# Check for taxa with 0 column sums (add 1 to everything if this happens)
		badTaxa <- which(colSums(tempData) == 0)
		if(length(badTaxa) != 0)
			tempData <- tempData + 1
		
		# Handle having 1 sample
		if(numSub == 1){
			tempParam1[,1] <- colnames(tempData)
			tempParam1[,3] <- unlist(tempData[1,]/sum(tempData))
			tempParam1[,4] <- NA
			tempParam1[,5] <- tempParam1[,3]
			tempParam1[,6] <- tempParam1[,3]
			tempParam1 <- tempParam1[order(tempParam1[,1]),]
			
			tempTheta1 <- c(0, NA)
			
			tempParam2 <- tempParam1
			tempTheta2 <- tempTheta1
		}else{	
			# Get the MoM and MLE for every taxa
			fsum <- dirmult::dirmult.summary(tempData, dirmult::dirmult(tempData, trace=FALSE))
			tempTheta <- fsum[nrow(fsum),]
			fsum <- fsum[-nrow(fsum),] 
			
			# Turn the summary into a data frame we can plot from
			tempParam1[,1] <- rownames(fsum)
			tempParam1[,3] <- fsum$MLE
			tempParam1[,4] <- fsum$se.MLE
			tempTheta1 <- tempTheta[,2:3]
			
			tempParam2[,1] <- rownames(fsum)
			tempParam2[,3] <- fsum$MoM
			tempParam2[,4] <- fsum$se.MOM
			tempTheta2 <- tempTheta[,4:5]
			
			# Calc Upper and Lower bounds for CI
			minSubj <- min(sapply(group.data, function(x) nrow(x)))
			if(minSubj < 30){
				val <- stats::qt(0.5 + conf *0.5, df=minSubj-1)
			}else{
				val <- stats::qnorm(0.5 + conf*0.5)
			}
			
			tempParam1[,5] <- tempParam1[,3] + val*tempParam1[,4]
			tempParam1[,6] <- tempParam1[,3] - val*tempParam1[,4]
			
			tempParam2[,5] <- tempParam2[,3] + val*tempParam2[,4]
			tempParam2[,6] <- tempParam2[,3] - val*tempParam2[,4]
		}
		
		# Save outside of loop
		allParamsMLE <- rbind(allParamsMLE, tempParam1)
		thetaMLE[i,] <- c(grpNames[i], tempTheta1)
		
		allParamsMOM <- rbind(allParamsMOM, tempParam2)
		thetaMOM[i,] <- c(grpNames[i], tempTheta2)
	}
	colnames(allParamsMLE) <- c("Taxa", "Group", "PI", "SE", "Upper", "Lower")
	colnames(thetaMLE) <- c("Group", colnames(tempTheta1))
	colnames(allParamsMOM) <- c("Taxa", "Group", "PI", "SE", "Upper", "Lower")
	colnames(thetaMOM) <- c("Group", colnames(tempTheta2))
	
	# Make sure none of our error bars go over 100 or below 0
	allParamsMLE$Upper <- ifelse(allParamsMLE$Upper > 1, 1, allParamsMLE$Upper)
	allParamsMLE$Lower <- ifelse(allParamsMLE$Lower < 0, 0, allParamsMLE$Lower)
	allParamsMOM$Upper <- ifelse(allParamsMOM$Upper > 1, 1, allParamsMOM$Upper)
	allParamsMOM$Lower <- ifelse(allParamsMOM$Lower < 0, 0, allParamsMOM$Lower)
	
	# Factor the data so it stays in the right order
	allParamsMLE$Group <- factor(allParamsMLE$Group, levels=grpNames)
	allParamsMLE$Taxa <- factor(allParamsMLE$Taxa, levels=unique(colnames(group.data[[1]])))
	allParamsMOM$Group <- factor(allParamsMOM$Group, levels=grpNames)
	allParamsMOM$Taxa <- factor(allParamsMOM$Taxa, levels=unique(colnames(group.data[[1]])))
	
	MLE <- list(params=allParamsMLE, theta=thetaMLE)
	MOM <- list(params=allParamsMOM, theta=thetaMOM)
	
	return(list(MLE=MLE, MOM=MOM))
}

Test.Paired <- function(group.data, numPerms=1000, parallel=FALSE, cores=3){
	if(missing(group.data))
		stop("group.data is missing.")
	
	if(length(group.data) != 2)
		stop("group.data must have exactly 2 data sets.")
	
	if(numPerms <= 0)
		stop("The number of permutations must be an integer greater than 0.")
	
	# Make sure we have the same columns
	if(ncol(group.data[[1]]) != ncol(group.data[[2]])){
		warning("Group columns do not match, running formatDataSets.")
		group.data <- formatDataSets(group.data)
	}
	
	# Check they have the same number of subjects
	numSub <- nrow(group.data[[1]])
	if(numSub != nrow(group.data[[2]]))
		stop("Groups must have the same number of subjects.")
	
	# Check row names match
	rNames1 <- rownames(group.data[[1]])
	rNames2 <- rownames(group.data[[2]])
	if(!all(rNames1 == rNames2)){ # check names in the same order
		if(all(rNames1 %in% rNames2)){ # check names match in wrong order
			group.data[[1]] <- group.data[[1]][order(rNames1),]
			group.data[[2]] <- group.data[[2]][order(rNames2),]
		}else{
			warning("Subject names do not match, assuming data is ordered correctly.")
		}
	}
	
	# Turn into abundances
	group.data[[1]] <- group.data[[1]]/rowSums(group.data[[1]])
	group.data[[2]] <- group.data[[2]]/rowSums(group.data[[2]])
	
	# Merge data1 and data2 together
	dataComb <- rbind(group.data[[1]], group.data[[2]])
	
	# Get the differences between the groups
	dataDiff <- group.data[[1]] - group.data[[2]]
	meanDiff <- apply(dataDiff, 2, mean)
	
	# Calculate the sum of squares
	obsDiff <- sum(meanDiff^2)
	
	# Permute the group membership
	if(parallel){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		
		tryCatch({ 
					permDiffs <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE) %dopar%{
						# Randomly swap group membership by reverseing difference sign
						swaps <- sample(c(1, -1), numSub, replace=TRUE)
						dataDiffTemp <- dataDiff * swaps
						meanDiffTemp <- apply(dataDiffTemp, 2, mean)
						
						# Calculate the sum of squares
						obsDiffTemp <- sum(meanDiffTemp^2)
						
						return(obsDiffTemp)
					}	
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		permDiffs <- rep(0, numPerms)
		for(i in 1:numPerms){ 	
			# Randomly swap group membership by reverseing difference sign
			swaps <- sample(c(1, -1), numSub, replace=TRUE)
			dataDiffTemp <- dataDiff * swaps
			meanDiffTemp <- apply(dataDiffTemp, 2, mean)
			
			# Calculate the sum of squares
			permDiffs[i] <- sum(meanDiffTemp^2)
		}
	}	
	
	# Calculate pvalue
	pval <- (sum(permDiffs >= obsDiff) + 1)/(numPerms + 1)
	
	return(pval)
}

DM.Rpart <- function(data, covars, plot=TRUE, minsplit=1, minbucket=1, cp=0, numCV=10, numCon=0, parallel=FALSE, cores=3, use1SE=FALSE, lowerSE=TRUE){
	if(missing(data) || missing(covars))
		stop("data and/or covars are missing.")
	
	if(numCV < 2){
		ret <- DM.Rpart.Base(data, covars, plot, minsplit, minbucket, cp)
	}else if(numCon < 2){
		ret <- DM.Rpart.CV(data, covars, plot, minsplit, minbucket, cp, numCV, parallel, cores, use1SE, lowerSE)
	}else{
		ret <- DM.Rpart.CV.Consensus(data, covars, plot, minsplit, minbucket, cp, numCV, numCon, parallel, cores, use1SE, lowerSE)
	}
	
	return(ret)
}

DM.Rpart.Base <- function(data, covars, plot=TRUE, minsplit=1, minbucket=1, cp=0){
	if(missing(data) || missing(covars))
		stop("data and/or covars are missing.")
	
	# Set the methods to use and call rpart
	methods <- list(init=rpartInit, eval=rpartEval, split=rpartSplit)
	rpartRes <- rpart::rpart(as.matrix(data) ~., data=covars, method=methods, minsplit=minsplit, minbucket=minbucket, cp=cp)
	
	cpInfo <- rpartRes$cptable
	size <- cpInfo[nrow(cpInfo), 2] + 1
	
	# Get split info from best tree
	splits <- NULL
	if(size > 1)
		splits <- rpartCS(rpartRes)
	
	# Plot the rpart results
	if(plot)
		suppressWarnings(rpart.plot::rpart.plot(rpartRes, type=2, extra=101, box.palette=NA, branch.lty=3, shadow.col="gray", nn=FALSE))
	
	return(list(cpTable=cpInfo, fullTree=rpartRes, bestTree=rpartRes, errorRate=NULL, size=size, splits=splits))
}

DM.Rpart.CV <- function(data, covars, plot=TRUE, minsplit=1, minbucket=1, cp=0, numCV=10, parallel=FALSE, cores=3, use1SE=FALSE, lowerSE=TRUE){
	if(missing(data) || missing(covars))
		stop("data and/or covars are missing.")
	if(numCV < 2)
		stop("numCV must be at least 2.")
	
	# Run initial Rpart
	rpartBase <- DM.Rpart.Base(data, covars, FALSE, minsplit, minbucket, cp)
	rpartRes <- rpartBase$fullTree
	
	# Check for a valid starting tree
	if(nrow(rpartRes$cptable) == 1){
		warning("No splits in the data.")
		return(rpartBase)
	}
	
	cvRes <- rpartCV(data, covars, rpartRes, minsplit, minbucket, numCV, parallel, cores)
	
	# Calculate the best tree
	ciInfo <- as.data.frame(cvRes$ciInfo)
	
	# Find the tree with the lowest MSE
	minMSE <- min(ciInfo$MSE)
	lowestMSELoc <- which(ciInfo$MSE == minMSE)[1]
	
	# Find which trees are within 1 SE of the lowest mse tree
	cutoffU <- ciInfo$MSE[lowestMSELoc] + ciInfo$SE[lowestMSELoc]
	cutoffL <- ciInfo$MSE[lowestMSELoc] - ciInfo$SE[lowestMSELoc]
	ciInfo$within1SE <- ifelse(ciInfo$MSE <= cutoffU & ciInfo$MSE >= cutoffL, 1, 0)   
	
	if(use1SE){	
		# Find the smallest/biggest tree within 1 SE
		within <- which(ciInfo$within1SE == 1)
		if(lowerSE){
			bestTreeLoc <- min(within)
		}else{
			bestTreeLoc <- max(within)
		}
	}else{
		bestTreeLoc <- lowestMSELoc
	}
	
	# Pull out the best tree
	size <- ciInfo[bestTreeLoc, 2] + 1
	best <- rpart::prune(rpartRes, cp=ciInfo[bestTreeLoc, 1])
	
	# Get split info from best tree
	splits <- NULL
	if(size > 1)
		splits <- rpartCS(best)
	
	if(plot)
		suppressWarnings(rpart.plot::rpart.plot(best, type=2, extra=101, box.palette=NA, branch.lty=3, shadow.col="gray", nn=FALSE))
	
	return(list(cpTable=ciInfo, fullTree=rpartRes, bestTree=best, errorRate=cvRes$errorRate, size=size, splits=splits))
}

DM.Rpart.CV.Consensus <- function(data, covars, plot=TRUE, minsplit=1, minbucket=1, cp=0, numCV=10, numCon=100, parallel=FALSE, cores=3, use1SE=FALSE, lowerSE=TRUE){
	if(missing(data) || missing(covars))
		stop("data and/or covars are missing.")
	if(numCV < 2)
		stop("numCV must be at least 2.")
	if(numCon < 2)
		stop("numCon must be at least 2.")
	
	if(parallel){
		cl <- parallel::makeCluster(min(cores, numCon)) 
		doParallel::registerDoParallel(cl)
		
		tryCatch({
					results <- foreach::foreach(i=1:numCon, .combine=append, .multicombine=FALSE, .inorder=FALSE, .errorhandling="pass", .packages=c("rpart", "HMP")) %dopar%{
						cvList <- DM.Rpart.CV(data, covars, FALSE, minsplit, minbucket, cp, numCV, FALSE, 1, use1SE, lowerSE)
						return(list(cvList))
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		results <- vector("list", numCon)
		for(i in 1:numCon)	
			results[[i]] <- DM.Rpart.CV(data, covars, FALSE, minsplit, minbucket, cp, numCV, FALSE, 1, use1SE, lowerSE)
	}
	
	# Combine cv results
	MSETab <- do.call("cbind", lapply(results, function(x){x$cpTable[,4]}))
	rankTab <- do.call("cbind", lapply(results, function(x){x$cpTable[,8]}))
	ciInfo <- cbind(
			results[[1]]$cpTable[,1:3],
			"MeanMSE"=rowMeans(MSETab), 
			"sdMSE"=apply(MSETab, 1, sd), 
			"MeanRank"=rowMeans(rankTab), 
			"sdRank"=apply(rankTab, 1, sd)
	)
	
	# Find the tree with the lowest MSE
	minMSE <- min(ciInfo$MeanMSE)
	bestTreeLoc <- which(ciInfo$MeanMSE == minMSE)[1]
	
	# Pull out the best tree
	size <- ciInfo[bestTreeLoc, 2] + 1
	best <- rpart::prune(results[[1]]$fullTree, cp=ciInfo[bestTreeLoc, 1])
	
	# Get split info from best tree
	splits <- NULL
	if(size > 1)
		splits <- rpartCS(best)
	
	if(plot)
		suppressWarnings(rpart.plot::rpart.plot(best, type=2, extra=101, box.palette=NA, branch.lty=3, shadow.col="gray", nn=FALSE))
	
	return(list(cpTable=ciInfo, fullTree=results[[1]]$fullTree, bestTree=best, errorRate=NULL, size=size, splits=splits))
}

Gen.Alg <- function(data, covars, iters=50, popSize=200, earlyStop=0, dataDist="euclidean", covarDist="gower", 
		verbose=FALSE, plot=TRUE, minSolLen=NULL, maxSolLen=NULL, custCovDist=NULL, penalty=0){
	if(missing(data) || (missing(covars) && is.null(custCovDist)))
		stop("data and/or covars are missing.")
	
	# Check for any bad numbers
	if(iters <= 0)
		stop("iters must be an integer greater than 0")
	if(popSize <= 0)
		stop("popSize must be an integer greater than 0")
	if(earlyStop < 0)
		stop("earlyStop must be an integer greater than or equal to 0")
	if(penalty < 0 || penalty > 1)
		stop("penalty must be between 0 and 1")
	
	# Check distances
	if(dataDist != "euclidean" && dataDist != "gower")
		stop("data.dist must be euclidean or gower.")
	if(covarDist != "euclidean" && covarDist != "gower")
		stop("covars.dist must be euclidean or gower.")
	
	# Define size
	size <- ncol(data)
	
	# Check stopping rules
	if(!is.null(minSolLen))
		if(minSolLen < 0 || minSolLen >= size)
			stop("minSolLen must be 0 or greater and less than the number of columns in data.")
	if(!is.null(maxSolLen))
		if(maxSolLen <= 0 || maxSolLen > size)
			stop("maxSolLen must be greater than 0 and less than or equal to the number columns in data.")
	if(!is.null(maxSolLen) && !is.null(minSolLen))
		if(maxSolLen < minSolLen)
			stop("maxSolLen must be bigger than minSolLen.")
	
	# Define some variables for use in the GA loop
	mutationChance <- 1/(size+1)
	elitism <- floor(popSize/5)
	evalSumm <- matrix(NA, iters, 6)
	newPopSize <- popSize - elitism
	newPopulation <- matrix(NA, newPopSize, size)
	parentProb <- stats::dnorm(1:popSize, mean=0, sd=(popSize/3))
	
	if(verbose){
		print("X. Current Step : Current Time Taken")
		runningTime <- proc.time()
		print(paste("1. Calculating Distances:", round((proc.time() - runningTime)[3], 3)))
	}
	
	# Set up our base distance matrix
	if(is.null(custCovDist)){
		covarDists <- vegan::vegdist(covars, covarDist)
	}else{
		covarDists <- custCovDist
	}
	
	# Get each columns distance contribution
	colDists <- vector("list", ncol(data))
	for(i in 1:ncol(data))
		colDists[[i]] <- vegan::vegdist(data[,i], dataDist)
	if(dataDist == "euclidean")
		colDists <- lapply(colDists, function(x) x^2)
	
	if(verbose)
		print(paste("2. Creating Starting Data:", round((proc.time() - runningTime)[3], 3)))
	
	# Create our starting data
	population <- gaCreation(data, popSize)
	
	if(verbose)
		print(paste("3. Scoring Starting Data:", round((proc.time() - runningTime)[3], 3)))
	
	# Score and sort
	evalVals <- rep(NA, popSize)
	for(e in 1:popSize)
		evalVals[e] <- gaScoring(population[e,], covarDists, colDists, dataDist, penalty, minSolLen, maxSolLen)
	population <- population[order(evalVals, decreasing=TRUE),]
	bestScoreValue <- max(evalVals)
	bestScoreCounter <- 0
	
	if(verbose)
		print(paste("4. Running Iterations:", round((proc.time() - runningTime)[3], 3)))
	
	# Run GA
	ptr <- proc.time()
	for(i in 1:iters){
		if(verbose){
			if(i %% round(iters/10) == 0)
				print(paste("Iteration - ", i, ": ", round((proc.time() - runningTime)[3], 3), sep=""))
		}
		# Cross over to fill rest of new population
		for(child in 1:newPopSize){
			parentIDs <- sample(1:popSize, 2, prob=parentProb)
			parents <- population[parentIDs,]
			crossOverPoint <- sample(0:size, 1)
			if(crossOverPoint == 0){
				newPopulation[child,] <- parents[2,]
			}else if(crossOverPoint == size){
				newPopulation[child,] <- parents[1,]
			}else{
				newPopulation[child,] <- c(parents[1,][1:crossOverPoint], parents[2,][(crossOverPoint+1):size])
			}
		}
		
		# Mutate all but elite
		if(mutationChance > 0){
			population[(elitism+1):popSize,] <- apply(newPopulation, 2, function(x){ifelse(stats::runif(newPopSize) < mutationChance, 1-x, x)})
		}else{
			population[(elitism+1):popSize,] <- newPopulation
		}
		
		# Score and sort our new solutions
		for(e in 1:popSize)
			evalVals[e] <- gaScoring(population[e,], covarDists, colDists, dataDist, penalty, minSolLen, maxSolLen)
		population <- population[order(evalVals, decreasing=TRUE),]
		evalSumm[i,] <- summary(evalVals)
		
		# Check if we want to stop early
		if(bestScoreValue == max(evalVals)){
			bestScoreCounter <- bestScoreCounter + 1
		}else{
			bestScoreCounter <- 0
			bestScoreValue <- max(evalVals)
		}
		
		if(bestScoreCounter == earlyStop && earlyStop != 0)
			break
	}	
	gaTime <- (proc.time() - ptr)[3]
	
	if(verbose)
		print(paste("5. Prettying Results", round((proc.time() - runningTime)[3], 3)))
	
	# Pretty up our data for returning
	rownames(population) <- paste("Solution", 1:nrow(population))
	colnames(population) <- colnames(data)
	rownames(evalSumm) <- paste("Iteration", 1:nrow(evalSumm))
	colnames(evalSumm) <- c("Worst", "25%ile", "Median", "Mean", "75%ile", "Best")
	
	evalVals <- matrix(evalVals[order(evalVals, decreasing=TRUE)], 1, length(evalVals))
	colnames(evalVals) <- paste("Solution ", 1:length(evalVals))
	rownames(evalVals) <- "Score"
	
	# Get selected columns using a consensus
	selIndex <- which(population[1,] == 1)
	sel <- colnames(data)[selIndex]
	
	# Get the nonselected columns
	nonSel <- colnames(data)[-selIndex]
	
	# Plot scoring summary
	if(plot)
		gaPlot(evalSumm)
	
	return(list(scoreSumm=evalSumm, solutions=population, scores=evalVals, time=gaTime, selected=sel, nonSelected=nonSel, selectedIndex=selIndex))
}

Gen.Alg.Consensus <- function(data, covars, consensus=.5, numRuns=10, parallel=FALSE, cores=3, ...){
	if(missing(data) || missing(covars))
		stop("data and/or covars are missing.")
	
	if(consensus <= 0 || consensus > 1)
		stop("consensus must be greater than 0 and equal or less than 1")
	
	# Run the GA X times
	if(parallel){
		cl <- parallel::makeCluster(min(cores, numRuns)) 
		doParallel::registerDoParallel(cl)
		
		tryCatch({
					gaRes <- foreach::foreach(i=1:numRuns, .combine=list, .multicombine=TRUE, .inorder=FALSE, .packages=c("vegan", "HMP")) %dopar%{
						tempResults <- Gen.Alg(data, covars, plot=FALSE, verbose=FALSE, ...)
						return(tempResults)
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		gaRes <- vector("list", numRuns)
		for(i in 1:numRuns)
			gaRes[[i]] <- Gen.Alg(data, covars, plot=FALSE, verbose=FALSE, ...)
	}
	
	# Get all the best solutions
	bestSols <- sapply(gaRes, function(x){x$solutions[1,]})
	
	# Get the consensus solution vector
	consSol <- (rowSums(bestSols) >= (numRuns * consensus)) * 1
	
	# Get the selected Index's
	selInd <- which(consSol == 1)
	
	return(list(solutions=bestSols, consSol=consSol, selectedIndex=selInd))
}



### ~~~~~~~~~~~~~~~~~~~~~
### Plot functions
### ~~~~~~~~~~~~~~~~~~~~~
Barchart.data <- function(data, title="Taxa Proportions"){
	if(missing(data))
		stop("data missing.")
	
	dataProp <- apply(data, 1, function(x){x/sum(x)})
	
	barplot(dataProp, col=rainbow(ncol(data)), horiz=TRUE, 
			main=title, axisnames=FALSE, font.main=20, font.sub=16)
}

Plot.PI <- function(estPi, errorBars=TRUE, logScale=FALSE, main="PI Vector", ylab="Fractional Abundance"){
	if(missing(estPi))
		stop("estPi is missing.")
	
	# Move title to the middle
	ggplot2::theme_update(plot.title=ggplot2::element_text(hjust=0.5))
	
	# Make the base plot
	piPlot <- ggplot2::ggplot(estPi$params, ggplot2::aes_string(y="PI", x="Taxa", colour="Group")) +
			ggplot2::geom_point() + 
			ggplot2::theme(legend.position = "top", text=ggplot2::element_text(size=15)) +
			ggplot2::labs(title=main, y=ylab, x="") +
			ggplot2::theme(axis.text.x=ggplot2::element_text(hjust=1, angle=45, size=10))
	
	# Add error bars
	if(errorBars){
		piPlot <- piPlot + ggplot2::geom_errorbar(ggplot2::aes_string(ymax="Upper", ymin="Lower"))
	}else{
		piPlot <- piPlot + ggplot2::geom_line(ggplot2::aes_string(group="Group"))
	}
	
	# Do log scaling
	if(logScale)
		piPlot <- piPlot + ggplot2::scale_y_log10()
	
	if(logScale)
		piPlot <- piPlot + ggplot2::labs(y=paste(ylab, "(Logged)"))
	
	print(piPlot)
}

Plot.MDS <- function(group.data, main="Group MDS", retCords=FALSE){
	if(missing(group.data))
		stop("group.data is missing.")
	
	numGroups <- length(group.data)
	
	# Make sure we have the same columns
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa)){
		warning("Group columns do not match, running formatDataSets.")
		group.data <- formatDataSets(group.data)
	}
	
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
	plot(loc, pch=16, ylab="MDS 2", xlab="MDS 1", col=cols, main=main)
	legend("topright", legend=grpNames, pch=15, col=availCols)
	
	if(retCords)
		return(loc)
}

Plot.RM.Barchart <- function(group.data, groups, times, plotByGrp=TRUE, col=NULL, conf=.95){
	if(missing(group.data) || missing(groups) || missing(times))
		stop("group.data, groups and/or times are missing.")
	
	numSamps <- length(group.data)
	
	### Get the pi params
	myEst <- Est.PI(group.data, conf)
	params <- myEst$MLE$params
	
	### Add the group and time information to the params
	myGroups <- NULL
	myTimes <- NULL
	for(i in 1:numSamps){
		myGroups <- c(myGroups, rep(groups[i], ncol(group.data[[1]])))
		myTimes <- c(myTimes, rep(times[i], ncol(group.data[[1]])))
	}
	params$Grp <- as.character(myGroups)
	params$Time <- as.character(myTimes)
	
	if(is.null(col))
		col <- rainbow(length(unique(params$Taxa)))
	
	if(plotByGrp){
		lattice::barchart(params$PI ~ params$Time | paste("Group", params$Grp), 
				ylab="Fractional Abundance", xlab="Time", 
				stack=TRUE, groups=params$Taxa, col=col,
				key=list(
						text=list(levels(params$Taxa)), 
						points=list(pch=19, col=col),
						columns=5
				)
		)
	}else{
		lattice::barchart(params$PI ~ params$Grp | paste("Time", params$Time),
				ylab="Fractional Abundance", xlab="Time", 
				stack=TRUE, groups=params$Taxa, col=col,
				key=list(
						text=list(levels(params$Taxa)), 
						points=list(pch=19, col=col),
						columns=5
				)
		)
	}
}

Plot.RM.Dotplot <- function(group.data, groups, times, errorBars=TRUE, col=NULL, conf=.95, alpha=1){
	if(missing(group.data) || missing(groups) || missing(times))
		stop("group.data, groups and/or times are missing.")
	
	numSamps <- length(group.data)
	numGrps <- length(unique(groups))
	
	### Get the pi params
	myEst <- Est.PI(group.data, conf)
	params <- myEst$MLE$params
	
	### Add the group and time information to the params
	myGroups <- NULL
	myTimes <- NULL
	for(i in 1:numSamps){
		myGroups <- c(myGroups, rep(groups[i], ncol(group.data[[1]])))
		myTimes <- c(myTimes, rep(times[i], ncol(group.data[[1]])))
	}
	params$Grp <- as.character(myGroups)
	params$Time <- as.character(myTimes)
	
	if(is.null(col))
		col <- rainbow(numGrps)
	### Add alpha to the colors
	col <- apply(sapply(col, grDevices::col2rgb)/255, 2, function(x){grDevices::rgb(x[1], x[2], x[3], alpha=alpha)})  
	
	if(errorBars){
		lattice::dotplot(params$Taxa ~ params$PI | paste("Time", params$Time), 
				pch=19, groups=params$Grp, col=col,
				ylab="Taxa", xlab="Fractional Abundance", 
				panel=lattice::panel.superpose, 
				panel.groups=function(x, y, subscripts, col, ...){
					lattice::panel.xyplot(x, y, ...)
					lattice::panel.segments(params$Lower[subscripts], y, params$Upper[subscripts], y, col=col)
				},
				key=list(
						text=list(as.character(unique(params$Grp))), 
						points=list(pch=19, col=col)
				)
		)
	}else{
		lattice::dotplot(params$Taxa ~ params$PI | paste("Time", params$Time), 
				pch=19, groups=params$Grp, col=col,
				ylab="Taxa", xlab="Fractional Abundance", 
				key=list(
						text=list(as.character(unique(params$Grp))), 
						points=list(pch=19, col=col)
				)
		)
	}
}

Plot.Theta <- function(estPi, main="Theta by Group"){
	if(missing(estPi))
		stop("estPi is missing.")
	
	# Create the theta table
	theta <- estPi$theta
	thetaci <- cbind(
			theta, 
			lci = theta[,2] - 1.96*theta[,3],
			lci = theta[,2] + 1.96*theta[,3]
	)
	thetaci <- thetaci[order(thetaci[2]),]
	
	xlim <- range(thetaci[,c(4, 5)])
	
	# Plot the tornado plot
	plot(thetaci[,2], 1:nrow(theta), pch=16, yaxt="n", xlim=xlim,
			main=main, ylab="", xlab="Theta +/- 95% CI")
	
	grid(ny=15, lwd=2)
	axis(2, at=1:nrow(theta), labels=thetaci[,1])
	
	for (i in 1:nrow(theta)) 
		lines(c(thetaci[i, 4], thetaci[i, 5]), c(i, i))
}

Plot.MDS.wPI <- function(group.data, pi, main="Group MDS", retCords=FALSE){
	if(missing(group.data) || missing(pi))
		stop("group.data and/or pi is missing.")
	
	numGroups <- length(group.data)
	
	# Make sure we have the same columns
	taxaCounts <- sapply(group.data, ncol)
	numTaxa <- taxaCounts[1]
	if(any(taxaCounts != numTaxa)){
		warning("Group columns do not match, running formatDataSets.")
		group.data <- formatDataSets(group.data)
	}
	
	# Make sure we have group names
	if(is.null(names(group.data))){
		grpNames <- paste("Data Set", 1:numGroups)
	}else{
		grpNames <- names(group.data)
	}
	
	# Merge all the data sets together
	mData <- do.call("rbind", group.data)
	numSubs <- nrow(mData)
	
	# Get the pi data, label, sort and combine
	numTaxa <- ncol(mData)
	piData <- matrix(pi$params$PI, length(group.data), numTaxa, byrow=TRUE)
	colnames(piData) <- pi$params$Taxa[1:numTaxa]
	rownames(piData) <- unique(pi$params$Group)
	piData <- piData[, colnames(mData)]
	mData <- rbind(mData, piData)
	
	# Get their mds location
	loc <- getBC(mData)
	
	# Set color
	availCols <- rainbow(numGroups)
	cols <- NULL
	for(i in 1:numGroups)
		cols <- c(cols, rep(availCols[i], nrow(group.data[[i]])))
	cols <- c(cols, availCols)
	
	# Set size and shape
	pch <- c(rep(16, numSubs), rep(17, numGroups))
	cex <- c(rep(1, numSubs), rep(2, numGroups))
	
	# Plot MDS
	plot(loc, pch=pch, ylab="MDS 2", xlab="MDS 1", col=cols, main=main, cex=cex)
	legend("topright", legend=grpNames, pch=15, col=availCols)
	
	if(retCords)
		return(loc)
}


### ~~~~~~~~~~~~~~~~~~~~~
### Filter functions
### ~~~~~~~~~~~~~~~~~~~~~
formatDataSets <- function(group.data){
	if(missing(group.data))
		stop("group.data missing.")
	
	# Make sure we have more than 1 data set
	numGroups <- length(group.data)
	if(numGroups < 2)
		stop("At least 2 data sets are required.")
	
	# Remove all 0 total read subjects from the data
	group.data <- lapply(group.data, function(x){x[rowSums(x) > 0,, drop=FALSE]})
	
	# Merge all the data together
	dataNames <- vector("list", numGroups)
	newData <- NULL
	for(i in 1:length(group.data)){		
		tempData <- group.data[[i]]
		
		# Save the current row names
		dataNames[[i]] <- rownames(tempData)
		
		newData <- merge(newData, t(group.data[[i]]), by=0, all=TRUE)
		rownames(newData) <- newData[,1]
		newData <- newData[,-1]
	}
	
	# Remove any nas
	newData[is.na(newData)] <- 0
	newData <- t(newData)
	
	# Remove any all 0 columns and sort them
	newData <- newData[,colSums(newData) != 0, drop=FALSE]
	newData <- newData[,order(colSums(newData), decreasing=TRUE)]
	
	# Turn the data back into a list
	retData <- vector("list", numGroups)
	base <- 0
	for(i in 1:numGroups){
		retData[[i]] <- newData[(base+1):(nrow(group.data[[i]])+ base),]
		rownames(retData[[i]]) <- dataNames[[i]]
		
		base <- base + nrow(group.data[[i]])
	}
	
	names(retData) <- names(group.data)
	return(retData)
}

Data.filter <- function(data, order.type="data", minReads=0, numTaxa=NULL, perTaxa=NULL){
	if(missing(data))
		stop("data is missing.")
	if(tolower(order.type) != "data" && tolower(order.type) != "sample")
		stop(sprintf("'%s' not recognized, order.type must be 'data' or 'sample'", as.character(order.type)))
	
	# Check if numTaxa or perTaxa is being used
	if(!is.null(numTaxa) && !is.null(perTaxa))
		stop("numTaxa and perTaxa cannot be used at the same time")
	if(!is.null(numTaxa)){
		if(numTaxa > ncol(data) || numTaxa <= 0)
			stop(sprintf("numTaxa must be between 0 and %i.", ncol(data)))
	}
	if(!is.null(perTaxa)){
		if(perTaxa >= 1 || perTaxa <= 0)
			stop("perTaxa must be between 0 and 1.")
	}
	if(is.null(numTaxa) && is.null(perTaxa))
		numTaxa <- ncol(data)
	
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
		taxaAboveThrs <- which(cumSumReads > perTaxa)
		if(length(taxaAboveThrs) == 0){
			numTaxa <- 1
		}else{
			numTaxa <- min(taxaAboveThrs)
		}
	}
	
	if(numTaxa >= ncol(data)){
		retData <- data
	}else{
		# Pull out the taxa we want to collapse
		otherData <- data[,-c(1:numTaxa), drop=FALSE]
		
		# Put the data back together and relabel
		retData <- cbind(data[,1:numTaxa], Other=rowSums(otherData))
	}
	
	return(retData)
}



### ~~~~~~~~~~~~~~~~~~~~~
### MC functions
### ~~~~~~~~~~~~~~~~~~~~~
MC.ZT.statistics <- function(Nrs, numMC=10, fit, type="ha", siglev=0.05) {
	if(missing(Nrs) || missing(fit))
		stop("Nrs and/or fit missing.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	
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
	zpval <- (sum(z > qAlpha) + 1)/(length(z) + 1)
	tpval <- (sum(t > qAlpha) + 1)/(length(t) + 1)
	
	return(cbind(zpval, tpval))
}

MC.Xsc.statistics <- function(Nrs, numMC=10, fit, pi0=NULL, type="ha", siglev=0.05) {
	if(missing(Nrs) || missing(fit))
		stop("Nrs and/or fit missing.")
	if(is.null(pi0) && tolower(type) == "ha")
		stop("pi0 cannot be null with type 'ha'.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	
	# Get all the XSC values
	XscStatVector <- rep(0, numMC)
	for(i in 1:numMC)
		XscStatVector[i] <- Xsc.statistics.Hnull.Ha(Nrs, fit, type, pi0)
	
	# Remove NAs
	XscStatVector <- XscStatVector[!is.na(XscStatVector)]
	
	# Get a reference value from the real data
	qAlpha <- qchisq(p=(1-siglev), df=length(fit$pi)-1, ncp=0, lower.tail=TRUE)
	
	# Calculate pvalues
	pval <- (sum(XscStatVector > qAlpha) + 1)/(length(XscStatVector) + 1)
	
	return(pval)
}

MC.Xmc.statistics <- function(group.Nrs, numMC=10, pi0, group.pi, group.theta, type="ha", siglev=0.05) {
	if(missing(group.theta) || missing(pi0) || missing(group.Nrs))
		stop("group.Nrs, pi0 and/or group.theta missing.")
	if(missing(group.pi) && tolower(type) == "ha")
		stop("group.pi missing.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	
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
	pval <- (sum(XmcStatVector > qAlpha) + 1)/(length(XmcStatVector) + 1)
	
	return(pval)
}

MC.Xmcupo.statistics <- function(group.Nrs, numMC=10, pi0, group.pi, group.theta, type="ha", siglev=0.05) {
	if(missing(group.theta) || missing(group.Nrs))
		stop("group.Nrs and/or group.theta missing.")
	if(missing(group.pi) && tolower(type) == "ha")
		stop("group.pi missing.")
	if(missing(pi0) && tolower(type) == "hnull")
		stop("pi0 missing.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	
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
	pval <- (sum(XmcupoStatVector > qAlpha) + 1)/(length(XmcupoStatVector) + 1)
	
	return(pval)
}

MC.Xdc.statistics <- function(group.Nrs, numMC=10, alphap, type="ha", siglev=0.05, est="mom") {
	if(missing(alphap) || missing(group.Nrs))
		stop("group.Nrs and/or alphap  missing.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	if(tolower(est) != "mom" && tolower(est) != "mle")
		stop(sprintf("Est '%s' not found. Est must be 'mle' or 'mom'.", as.character(est)))
	
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
	pval <- (sum(XdcStatVector > qAlpha) + 1)/(length(XdcStatVector) + 1)
	
	return(pval)
}

MC.Xoc.statistics <- function(group.Nrs, numMC=10, group.alphap, type="ha", siglev=0.05) {
	if(missing(group.alphap) || missing(group.Nrs))
		stop("group.Nrs and/or group.alphap missing.")
	if(tolower(type) != "ha" && tolower(type) != "hnull")
		stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
	
	numGroups <- length(group.Nrs)	
	
	# Get all the Xoc values
	XocStatVector <- rep(0, numMC)
	for(i in 1:numMC)
		XocStatVector[i] <- Xoc.statistics.Hnull.Ha(group.Nrs, group.alphap, type)
	
	# Get a reference value from the real data
	qAlpha <- qchisq(p=(1-siglev), df=(numGroups-1), ncp=0, lower.tail=TRUE)
	
	# Calculate pvalues
	pval <- (sum(XocStatVector > qAlpha) + 1)/(length(XocStatVector) + 1)
	
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
	
	# Make sure we have the same columns
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa)){
		warning("Group columns do not match, running formatDataSets.")
		group.data <- formatDataSets(group.data)
		numTaxa <- ncol(group.data[[1]])
	}
	
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

Xoc.sevsample <- function(group.data, epsilon=10^(-4)){
	if(missing(group.data))
		stop("group.data missing.")
	
	# Make sure we have the same columns
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa)){
		warning("Group columns do not match, running formatDataSets.")
		group.data <- formatDataSets(group.data)
	}
	
	numGroups <- length(group.data)
	
	# Get Xoc and calculate pvalue
	Xoc <- Xoc.statistics(group.data, epsilon)
	pval <- 1-pchisq(q=Xoc, df=numGroups-1, ncp=0, lower.tail=TRUE)	
	
	sev.overd.test <- list("Xoc statistics"=Xoc, "p value"=pval)
	
	return(sev.overd.test)			
}

Xdc.sevsample <- function(group.data, epsilon=10^(-4), est="mom"){
	if(missing(group.data))
		stop("group.data missing.")
	if(tolower(est) != "mle" && tolower(est) != "mom")
		stop(sprintf("Est '%s' not found. Est must be 'mle' or 'mom'.", as.character(est)))
	
	# Make sure we have the same columns
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa)){
		warning("Group columns do not match, running formatDataSets.")
		group.data <- formatDataSets(group.data)
		numTaxa <- ncol(group.data[[1]])
	}
	
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
kullbackLeiber <- function(data, plot=TRUE, parallel=FALSE, cores=3){
	warning("This function has been spellchecked.  Please use 'Kullback.Leibler' instead.")
	kl <- Kullback.Leibler(data, plot, parallel, cores)
	return(kl)
}

loglikDM <- function(data, alphap){	
	data <- data[,colSums(data) != 0, drop=FALSE]
	alphap <- alphap[alphap != 0]
	
	ll <- sum(lgamma(rowSums(data)+1) + lgamma(sum(alphap)) - lgamma(sum(alphap)+rowSums(data))) + 
			sum(rowSums(lgamma(sweep(data, 2, alphap, "+")) - lgamma(data+1) - lgamma(t(replicate(nrow(data), alphap)))))
	
	return(ll)
}

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

getBC <- function(data){
	dataPer <- data/rowSums(data)
	
	### Remove any points with exactly the same values
	dups <- duplicated(dataPer)
	if(any(dups))
		dataPer <- dataPer[!dups,]
	
	bcDist <- vegan::vegdist(dataPer, method="bray")
	nonMetricMDS <- MASS::isoMDS(bcDist, trace=FALSE)
	mdsPoints <- vegan::postMDS(nonMetricMDS$points, bcDist)
	mds <- vegan::scores(mdsPoints)
	
	return(mds[, 1:2])
}

pioest <- function(group.data){
	if(missing(group.data))
		stop("data.groups missing.")
	
	# Make sure we have the same columns
	taxaCounts <- sapply(group.data, ncol)
	numTaxa	<- taxaCounts[1]
	if(any(taxaCounts != numTaxa)){
		warning("Group columns do not match, running formatDataSets.")
		group.data <- formatDataSets(group.data)
		numTaxa <- ncol(group.data[[1]])
	}
	
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
	pis <- pis[rowSums(pis) != 0,]
	
	# Calculate pi0
	pi0 <- rowSums(pis/xscs)/sum(1/xscs)
	
	names(pi0) <- colnames(group.data[[1]])
	
	return(pi0)
}



### ~~~~~~~~~~~~~~~~~~~~~
### ga functions
### ~~~~~~~~~~~~~~~~~~~~~
gaScoring <- function(indices, covarDists, colDists, distType, lambda, minSolLen, maxSolLen) {
	BAD_RETURN <- -2 # Return worse than cor could do
	
	numSel <- sum(indices)
	
	# Check if nothing is selected
	if(numSel == 0) 
		return(BAD_RETURN) 
	# Check if we dont have enough selected
	if(!is.null(minSolLen)) 
		if(numSel < minSolLen)
			return(BAD_RETURN)
	# Check if we dont have too many selected
	if(!is.null(maxSolLen)) 
		if(numSel > maxSolLen)
			return(BAD_RETURN) 
	
	edges <- which(indices==1)
	
	# Combine the distance matrices based on distance type
	if(distType == "gower"){
		combinedSolDists <- Reduce("+", colDists[edges])/sum(indices)
	}else if(distType == "euclidean"){
		combinedSolDists <- sqrt(Reduce("+", colDists[edges]))
	}
	
	# Get the correlation and penalize it based on the number of columns selected
	mycor <- stats::cor(combinedSolDists, covarDists)
	mycor <- mycor - (lambda * (sum(indices)/length(indices)))
	
	return(mycor)
}

gaCreation <- function(data, popSize){
	ZERO_TO_ONE_RATIO <- 10 # Ratio of 0 to 1s for the random data
	SUGGESTION_COUNT <- 10 # Number starting points we should make from the data
	
	size <- ncol(data)
	population <- matrix(NA, popSize, size)
	
	# Make 10 starting points as long as our popSize is > 10
	if(popSize >= SUGGESTION_COUNT){
		# Get a rough starting point
		rstart <- apply(data, 2, mean)
		
		# Use the rough difference to make starting solutions
		breaks <- seq(.05, 1, 1/SUGGESTION_COUNT)
		breakVals <- stats::quantile(rstart, breaks)
		suggestions <- matrix(0, length(breaks), length(rstart))
		for(i in 1:length(breaks))
			suggestions[i,] <- ifelse(rstart >= breakVals[i], 1, 0)
		
		population[1:SUGGESTION_COUNT,] <- suggestions
		numCreated <- SUGGESTION_COUNT
	}else{
		numCreated <- 0
	}
	
	# Fill any remaining population spots with random solutions
	if(popSize != SUGGESTION_COUNT){
		for(child in (numCreated+1):popSize) 
			population[child,] <- sample(c(rep(0, ZERO_TO_ONE_RATIO), 1), size, replace=TRUE)
	}
	
	return(population)
}

gaPlot <- function(evalSumm){
	plot(evalSumm[,4], type="l", ylab="Score", ylim=c(0, 1), lwd=2, main="Eval Scores by Iteration", xlab="Iteration")
	lines(evalSumm[,6], col="red", lwd=2)
	lines(evalSumm[,1], col="blue", lwd=2)
	legend("topleft", colnames(evalSumm)[c(4, 6, 1)], pch=16, col=c("black", "red", "blue"))
}



### ~~~~~~~~~~~~~~~~~~~~~
### rpart functions
### ~~~~~~~~~~~~~~~~~~~~~
rpartInit <- function(y, offset, parms, wt){
	hmp.pkg.env$EVAL_COUNT_RPART <- 1	# reset eval counts
	sfun <- function(yval, dev, wt, ylevel, digits){
		paste(" mean=", round(mean(yval), 3), sep="")
	}
	environment(sfun) <- .GlobalEnv
	list(y=y, parms=NULL, numresp=1, numy=ncol(y), summary=sfun)
}

rpartEval <- function(y, wt, parms){
	# Set a unique label
	label <- hmp.pkg.env$EVAL_COUNT_RPART
	hmp.pkg.env$EVAL_COUNT_RPART <- hmp.pkg.env$EVAL_COUNT_RPART + 1
	
	dev <- DM.MoM(y)$loglik * -1
	
	# Skip any infinite LL comparisons (makes lrt 0)
	if(dev == Inf || dev == -Inf)
		dev <- 0
	
	list(label=label, deviance=dev)
}

rpartSplit <- function(y, wt, x, parms, continuous){
	# Get initial LL
	LL <- DM.MoM(y)$loglik
	
	uniqX <- sort(unique(x))
	numUni <- length(uniqX) - 1
	
	# Determine what we are comparing
	if(continuous){
		numTests <- length(x) - 1
		dir <- rep(-1, numTests)
	}else{
		numTests <- numUni
		dir <- uniqX
	}
	
	# Run through every comparison
	LRT <- rep(0, numTests)
	for(i in 1:numUni){
		if(continuous){
			id <- which(x <= uniqX[i])
			grp1 <- y[id,, drop=FALSE]
			grp2 <- y[-id,, drop=FALSE]
		}else{
			grp1 <- y[x == uniqX[i],, drop=FALSE]
			grp2 <- y[x != uniqX[i],, drop=FALSE]
		}
		# Skip any 1 subject groups
		if(nrow(grp1) == 1 || nrow(grp2) == 1)
			next
		
		LLgrp1 <- DM.MoM(grp1)$loglik
		LLgrp2 <- DM.MoM(grp2)$loglik
		
		# Skip any infinite LL comparisons (makes lrt 0)
		if(LLgrp1 == Inf || LLgrp2 == Inf)
			next
		
		if(continuous){
			LRT[id[length(id)]] <- -2*(LL-LLgrp1-LLgrp2)
		}else{
			LRT[i] <- -2*(LL-LLgrp1-LLgrp2)
		}
	}
	ret <- list(goodness=LRT, direction=dir)
	
	return(ret)
}

rpartCV <- function(data, covars, rpartRes, minsplit, minbucket, numCV, parallel, cores){
	# Pull out cp info
	cpTable <- rpartRes$cptable
	numCPLvls <- nrow(cpTable)
	
	# New cp for pruning
	if(numCPLvls > 2){
		oldCPs <- cpTable[,1]
		
		cpTable[1, 1] <- Inf
		cpTable[numCPLvls, 1] <- 0
		for(m in 2:(numCPLvls-1))
			cpTable[m, 1] <- sqrt(oldCPs[m] * oldCPs[m-1])
	}
	
	# Set up groups for dropping data
	numSub <- nrow(data)
	dropNums <- cut(1:numSub, numCV, FALSE)
	dropGrps <- sample(dropNums, numSub)
	
	errorRate <- vector("list", numCV)
	if(parallel){
		cl <- parallel::makeCluster(min(cores, numCV)) 
		doParallel::registerDoParallel(cl)
		tryCatch({
					cvRes <- foreach::foreach(k=1:numCV, .combine=append, .multicombine=FALSE, .inorder=FALSE, .errorhandling="pass", .packages=c("rpart", "HMP")) %dopar%{
						cvRes <- rpartCVSingle(data, covars, k, cpTable, dropGrps, numCPLvls, minsplit, minbucket)
						return(list(cvRes))
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
		errorRate <- lapply(cvRes, function(x)x[[1]])
		subTree <- lapply(cvRes, function(x)x[[2]])
	}else{
		errorRate <- vector("list", numCV)
		subTree <- vector("list", numCV)
		for(k in 1:numCV){
			cvRes <- rpartCVSingle(data, covars, k, cpTable, dropGrps, numCPLvls, minsplit, minbucket)
			errorRate[[k]] <- cvRes[[1]]
			subTree[[k]] <- cvRes[[2]]
		}
	}
	
	# Calculate the square root of the errors
	error <- sapply(errorRate, sqrt)
	
	# Calculate CI of MSE
	ciInfo <- matrix(NA, numCPLvls, 4)
	for(j in 1:numCPLvls){
		ciInfo[j, 1] <- mean(error[j,])
		ciInfo[j, 2:3] <- rpartCI(error[j,], 0.95)
		ciInfo[j, 4] <- sd(error[j,])/sqrt(ncol(error))
	}
	
	ciInfo <- cbind(ciInfo, rank(ciInfo[,1]))
	colnames(ciInfo) <- c("MSE", "Lower", "Upper", "SE", "Rank")
	
	# Add ci info back into cp table
	cpTable2 <- cbind(rpartRes$cptable, ciInfo)
	
	return(list(subTree=subTree, errorRate=do.call("cbind", errorRate), ciInfo=cpTable2))
}

rpartCVSingle <- function(data, covars, cvNum, cpTable, dropGrps, numCPLvls, minsplit, minbucket){
	# Get location of data to drop
	dropLoc <- which(dropGrps == cvNum)
	
	# Seperate dropped data
	subData <- data[-dropLoc,, drop=FALSE]
	dropData <- data[dropLoc,, drop=FALSE]
	
	# Seperate dropped covars
	subCovs <- covars[-dropLoc,, drop=FALSE]
	dropCovs <- covars[dropLoc,, drop=FALSE]
	
	# Run rpart on smaller data
	subRpartRes <- DM.Rpart.Base(subData, subCovs, FALSE, minsplit, minbucket)$fullTree
	
	# Need to be abandance for later code
	subData <- subData/(rowSums(subData))
	dropData <- dropData/(rowSums(dropData))
	
	# Calculate relative error
	MSEn <- rep(NA, numCPLvls)
	subTree <- vector("list", numCPLvls)
	for(i in 1:numCPLvls){
		subTree[[i]] <-  rpart::prune(subRpartRes, cp=cpTable[i, 1])
		pruneSubTree <- subTree[[i]]
		subPre <- predict(pruneSubTree, newdata=subCovs, type="vector")
		dropPre <- predict(pruneSubTree, newdata=dropCovs, type="vector")
		
		## 1.a. Distance: new subject to the mean Taxa in the signed Terminal node
		tempDist <- 0
		for(j in 1:length(dropPre))
			tempDist <- tempDist + (dist(rbind(dropData[j,], colMeans(subData[subPre == dropPre[j],]))))^2
		
		MSEn[i] <- tempDist/length(dropPre)
	}
	names(MSEn) <- cpTable[,2] + 1
	
	return(list(errorRate=MSEn, subTree=subTree))
}

rpartCI <- function(vector, interval) {
	vec_sd <- sd(vector)
	numSamp <- length(vector)
	vec_mean <- mean(vector)
	
	# Error according to t distribution
	error <- qt((interval + 1)/2, df = numSamp - 1) * vec_sd/sqrt(numSamp)
	
	# Confidence interval as a vector
	res <- c("Lower"=vec_mean - error, "Upper"=vec_mean + error)
	return(res)
}

rpartCS <- function(fit) {
	# Pull out split information
	splitNames <- rownames(fit$splits)
	allVars <- colnames(attributes(fit$terms)$factors)  
	
	# Rename splits to fit into data frame
	rownames(fit$splits) <- 1:nrow(fit$splits)
	splits <- data.frame(fit$splits)
	splits$var <- splitNames
	splits$type <- ""
	splits$primary <- ""
	
	# Get the frame
	frame <- as.data.frame(fit$frame)
	frame$var <- as.character(frame$var)
	primeFr <- frame[frame$var != "<leaf>",]
	
	# Go through every node and check competing and surrogaet splits
	index <- 0
	for(i in 1:nrow(primeFr)){
		spltPrimName <- paste("Split", primeFr$yval[i], primeFr$var[i])
		
		# Fill in primary info
		index <- index + 1
		splits$type[index] <- "primary"
		splits$primary[index] <- spltPrimName
		
		# Check for competing splits
		if(primeFr$ncompete[i] > 0){
			for(j in 1:primeFr$ncompete[i]){
				index <- index + 1
				splits$type[index] <- "competing"
				splits$primary[index] <- spltPrimName
			}
		}
		
		# Check for surrogate splits
		if(primeFr$nsurrogate[i] > 0){
			for(j in 1:primeFr$nsurrogate[i]){
				index <- index + 1
				splits$type[index] <- "surrogate"
				splits$primary[index] <- spltPrimName
			}
		}
	}
	
	return(splits)
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
	pi0 <- colSums(t(pis)/xscs)/sum(1/xscs)
	
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
		tempTheta <- DM.MoM(group.data[[i]])$theta
		fit <- dirmult::dirmult(group.data[[i]], initscalar=(1-tempTheta)/tempTheta, epsilon=epsilon, trace=FALSE)
		
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
				tempTheta <- DM.MoM(x)$theta
				dirmult::dirmult(x, initscalar=(1-tempTheta)/tempTheta, epsilon=epsilon, trace=FALSE)$loglik
			}, epsilon=epsilon)
	
	# Get the fit assuming all in the same group
	groupDataC <- do.call(rbind, group.data)
	tempTheta <- DM.MoM(groupDataC)$theta
	groupFit <- dirmult::dirmult(groupDataC, initscalar=(1-tempTheta)/tempTheta, epsilon=epsilon, trace=FALSE)	
	
	# Calculate the xdc
	xdc <- -2*(groupFit$loglik-sum(logliks))
	
	return(xdc)
}

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


### ~~~~~~~~~~~~~~~~~~~~~
### OLD
### ~~~~~~~~~~~~~~~~~~~~~
DM.Rpart.Base.Old <- function(data, covars, plot=TRUE, minsplit=1, minbucket=1, cp=0){
	if(missing(data) || missing(covars))
		stop("data and/or covars are missing.")
	
	# Set the methods to use and call rpart
	methods <- list(init=rpartInit, eval=rpartEval, split=rpartSplitOld)
	rpartRes <- rpart::rpart(as.matrix(data) ~., data=covars, method=methods, minsplit=minsplit, minbucket=minbucket, cp=cp)
	
	cpInfo <- rpartRes$cptable
	size <- cpInfo[nrow(cpInfo), 2] + 1
	
	# Get split info from best tree
	splits <- NULL
	if(size > 1)
		splits <- rpartCS(rpartRes)
	
	# Plot the rpart results
	if(plot)
		suppressWarnings(rpart.plot::rpart.plot(rpartRes, type=2, extra=101, box.palette=NA, branch.lty=3, shadow.col="gray", nn=FALSE))
	
	return(list(cpTable=cpInfo, fullTree=rpartRes, bestTree=rpartRes, subTree=NULL, errorRate=NULL, size=size, splits=splits))
}

rpartSplitOld <- function(y, wt, x, parms, continuous){
	# Get initial LL
	LL <- DM.MoM(y)$loglik
	
	# Determine what we are comparing
	if(continuous){
		numTests <- length(x) - 1
		dir <- rep(-1, numTests)
	}else{
		uniqX <- sort(unique(x))
		numTests <- length(uniqX) - 1
		dir <- uniqX
	}
	
	# Run through every comparison
	LRT <- rep(0, numTests)
	for(i in 1:numTests){
		if(continuous){
			grp1 <- y[1:i,, drop=FALSE]
			grp2 <- y[-c(1:i),, drop=FALSE]
		}else{
			grp1 <- y[x == uniqX[i],, drop=FALSE]
			grp2 <- y[x != uniqX[i],, drop=FALSE]
		}
		# Skip any 1 subject groups
		if(nrow(grp1) == 1 || nrow(grp2) == 1)
			next
		
		LLgrp1 <- DM.MoM(grp1)$loglik
		LLgrp2 <- DM.MoM(grp2)$loglik
		
		# Skip any infinite LL comparisons (makes lrt 0)
		if(LLgrp1 == Inf || LLgrp2 == Inf)
			next
		
		LRT[i] <- -2*(LL-LLgrp1-LLgrp2)
	}
	ret <- list(goodness=LRT, direction=dir)
	
	return(ret)
}












