

source("HMP Functions.R")


#data <- read.csv("MicroTest.csv", row.names=1)
#covars <- read.csv("CovTest.csv", row.names=1)

saliva <- read.csv("saliva.csv", row.names=1)
throat <- read.csv("throat.csv", row.names=1)
tonsils <- read.csv("tonsils.csv", row.names=1)

dmrp_data <- read.csv("dmrp_data.csv", row.names=1)
dmrp_covars <- read.csv("dmrp_covars.csv", row.names=1)


### ~~~~~~~~~~~~~~~~~~~~~
### MC functions
### ~~~~~~~~~~~~~~~~~~~~~
MC.Xdc.statistics.Test <- function(){
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Get a list of dirichlet-multinomial parameters for the data
	fit.saliva <- DM.MoM(saliva) 
	fit.throat <- DM.MoM(throat)
	fit.tonsils <- DM.MoM(tonsils)
	
	### Set up the number of Monte-Carlo experiments
	### We use 1 for speed, should be at least 1,000
	numMC <- 1
	
	### Generate the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	nrsGrp1 <- rep(12000, 9)
	nrsGrp2 <- rep(12000, 11)
	nrsGrp3 <- rep(12000, 12)
	group.Nrs <- list(nrsGrp1, nrsGrp2, nrsGrp3)
	
	### Computing size of the test statistics (Type I error)
	alphap <- fit.saliva$gamma
	pval1 <- MC.Xdc.statistics(group.Nrs, numMC, alphap, "hnull")
	pval1
	
	### Computing Power of the test statistics (Type II error)
	alphap <- rbind(fit.saliva$gamma, fit.throat$gamma, fit.tonsils$gamma)
	pval2 <- MC.Xdc.statistics(group.Nrs, numMC, alphap)
	pval2
}

MC.Xoc.statistics.Test <- function(){
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Get a list of dirichlet-multinomial parameters for the data
	fit.saliva <- DM.MoM(saliva) 
	fit.throat <- DM.MoM(throat)
	fit.tonsils <- DM.MoM(tonsils)
	
	### Set up the number of Monte-Carlo experiments
	### We use 1 for speed, should be at least 1,000
	numMC <- 1
	
	### Generate the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	nrsGrp1 <- rep(12000, 9)
	nrsGrp2 <- rep(12000, 11)
	nrsGrp3 <- rep(12000, 12)
	group.Nrs <- list(nrsGrp1, nrsGrp2, nrsGrp3)
	
	### Computing size of the test statistics (Type I error)
	alphap <- fit.tonsils$gamma
	pval1 <- MC.Xoc.statistics(group.Nrs, numMC, alphap, "hnull")
	pval1
	
	\dontrun{
		### Computing Power of the test statistics (Type II error)
		alphap <- rbind(fit.saliva$gamma, fit.throat$gamma, fit.tonsils$gamma)
		pval2 <- MC.Xoc.statistics(group.Nrs, numMC, alphap, "ha")
		pval2
	}
}

MC.Xmc.statistics.Test <- function(){
	data(saliva)
	data(throat) 
	data(tonsils)
	
	### Get a list of dirichlet-multinomial parameters for the data
	fit.saliva <- DM.MoM(saliva) 
	fit.throat <- DM.MoM(throat)
	fit.tonsils <- DM.MoM(tonsils) 
	
	### Set up the number of Monte-Carlo experiments
	### We use 1 for speed, should be at least 1,000
	numMC <- 1 
	
	### Generate the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	nrsGrp1 <- rep(12000, 9)
	nrsGrp2 <- rep(12000, 11)
	group.Nrs <- list(nrsGrp1, nrsGrp2)
	
	group.theta <- c(0.01, 0.05)
	pi0 <- fit.saliva$pi
	
	### Computing size of the test statistics (Type I error)
	pval1 <- MC.Xmc.statistics(group.Nrs, numMC, pi0, group.theta=group.theta, type="hnull")
	pval1
	
	### Computing Power of the test statistics (Type II error)
	group.pi <- rbind(fit.throat$pi, fit.tonsils$pi)
	pval2 <- MC.Xmc.statistics(group.Nrs, numMC, pi0, group.pi, group.theta)
	pval2
}

MC.Xmcupo.statistics.Test <- function(){
	data(saliva) 
	data(throat) 
	data(tonsils)
	
	### Get a list of dirichlet-multinomial parameters for the data
	fit.saliva <- DM.MoM(saliva) 
	fit.throat <- DM.MoM(throat)
	fit.tonsils <- DM.MoM(tonsils) 
	
	### Set up the number of Monte-Carlo experiments
	### We use 1 for speed, should be at least 1,000
	numMC <- 1
	
	### Generate the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	Nrs1 <- rep(12000, 10)
	Nrs2 <- rep(12000, 19)
	group.Nrs <- list(Nrs1, Nrs2)
	
	group.theta <- c(fit.throat$theta, fit.tonsils$theta)
	pi0 <- fit.saliva$pi
	
	### Computing size of the test statistics (Type I error)
	group.theta <- c(fit.throat$theta, fit.tonsils$theta)
	pval1 <- MC.Xmcupo.statistics(group.Nrs, numMC, pi0, group.theta=group.theta, type="hnull")
	pval1
	
	### Computing Power of the test statistics (Type II error)
	group.pi <- rbind(fit.throat$pi, fit.tonsils$pi)
	pval2 <- MC.Xmcupo.statistics(group.Nrs, numMC, group.pi=group.pi, group.theta=group.theta)
	pval2
}

MC.Xsc.statistics.Test <- function(){
	data(saliva)
	data(throat) 
	data(tonsils)
	
	### Get a list of dirichlet-multinomial parameters for the data
	fit.saliva <- DM.MoM(saliva) 
	fit.throat <- DM.MoM(throat)
	fit.tonsils <- DM.MoM(tonsils) 
	
	### Set up the number of Monte-Carlo experiments
	### We use 1 for speed, should be at least 1,000
	numMC <- 1
	
	### Generate the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	nrs <- rep(15000, 25)
	
	### Computing size of the test statistics (Type I error)
	pval1 <- MC.Xsc.statistics(nrs, numMC, fit.tonsils, fit.saliva$pi, "hnull")
	pval1
	
	### Computing Power of the test statistics (Type II error)
	pval2 <- MC.Xsc.statistics(nrs, numMC, fit.throat, fit.tonsils$pi)
	pval2
}

MC.ZT.statistics.Test <- function(){
	data(saliva) 
	
	### Get a list of dirichlet-multinomial parameters for the data
	fit.saliva <- DM.MoM(saliva) 
	
	### Set up the number of Monte-Carlo experiments
	### We use 1 for speed, should be at least 1,000
	numMC <- 1
	
	### Generate the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	nrs <- rep(15000, 25)
	
	### Computing size of the test statistics (Type I error)
	pval1 <- MC.ZT.statistics(nrs, numMC, fit.saliva, "hnull") 
	pval1
	
	### Computing Power of the test statistics (Type II error)
	pval2 <- MC.ZT.statistics(nrs, numMC, fit.saliva)
	pval2
}



### ~~~~~~~~~~~~~~~~~~~~~
### Plot functions
### ~~~~~~~~~~~~~~~~~~~~~
Barchart.data.Test <- function(){
	data(saliva)
	
	Barchart.data(saliva)
}

Plot.PI.Test <- function(){
	\dontrun{
		data(saliva)
		data(throat)
		data(tonsils)
		
		### Combine the data sets into a single list
		group.data <- list(saliva, throat, tonsils)
		
		### Get PI using MLE with CI
		mle <- Est.PI(group.data)$MLE
		
		### Plot with Error Bars
		Plot.PI(mle)
		
		### Plot without Error Bars
		Plot.PI(mle, FALSE)
		
		### Plot with Error Bars and scaling
		Plot.PI(mle, TRUE, TRUE)
	}
}

Plot.MDS.Test <- function(){
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Combine the data sets into a single list
	group.data <- list(saliva, throat, tonsils)
	
	Plot.MDS(group.data)
}

Plot.RM.Barchart <- function(){
	data(saliva)
	data(throat)
	
	### Reduce the size of the data
	saliva <- Data.filter(saliva, numTaxa=9)
	throat <- Data.filter(throat, numTaxa=9)
	
	### Get the gamma value for the data
	saliva.gamma <- DM.MoM(saliva)$gamma
	throat.gamma <- DM.MoM(throat)$gamma
	mid.gamma <- (saliva.gamma + throat.gamma)/2
	
	### Generate a the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	nrs <- rep(10000, 20)
	
	### Create data sets to be our time series in a list
	group.data <- list(
			Dirichlet.multinomial(nrs, saliva.gamma),
			Dirichlet.multinomial(nrs, saliva.gamma),
			Dirichlet.multinomial(nrs, saliva.gamma),
			Dirichlet.multinomial(nrs, saliva.gamma),
			Dirichlet.multinomial(nrs, mid.gamma),
			Dirichlet.multinomial(nrs, throat.gamma)
	)
	names(group.data) <- c(
			"Group 1, Time 1", "Group 2, Time 1",
			"Group 1, Time 2", "Group 2, Time 2",
			"Group 1, Time 3", "Group 2, Time 3"
	)
	
	### Set the group and time information for each list element
	groups <- c(1, 2, 1, 2, 1, 2)
	times <- c(1, 2, 3, 1, 2, 3)
	
	### Plot the data by Group
	Plot.RM.Barchart(group.data, groups, times)
	
	### Plot the data by Time
	Plot.RM.Barchart(group.data, groups, times, FALSE)
}

Plot.RM.Dotplot <- function(){
	data(saliva)
	data(throat)
	
	### Reduce the size of the data
	saliva <- Data.filter(saliva, numTaxa=9)
	throat <- Data.filter(throat, numTaxa=9)
	
	### Get the gamma value for the data
	saliva.gamma <- DM.MoM(saliva)$gamma
	throat.gamma <- DM.MoM(throat)$gamma
	mid.gamma <- (saliva.gamma + throat.gamma)/2
	
	### Generate a the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	nrs <- rep(10000, 20)
	
	### Create data sets to be our time series in a list
	group.data <- list(
			Dirichlet.multinomial(nrs, saliva.gamma),
			Dirichlet.multinomial(nrs, saliva.gamma),
			Dirichlet.multinomial(nrs, saliva.gamma),
			Dirichlet.multinomial(nrs, saliva.gamma),
			Dirichlet.multinomial(nrs, mid.gamma),
			Dirichlet.multinomial(nrs, throat.gamma)
	)
	names(group.data) <- c(
			"Group 1, Time 1", "Group 2, Time 1",
			"Group 1, Time 2", "Group 2, Time 2",
			"Group 1, Time 3", "Group 2, Time 3"
	)
	
	### Set the group and time information for each list element
	groups <- c(1, 2, 1, 2, 1, 2)
	times <- c(1, 2, 3, 1, 2, 3)
	
	### Plot the data with error bars
	Plot.RM.Dotplot(group.data, groups, times)
	
	### Plot the data without error bars
	Plot.RM.Dotplot(group.data, groups, times, FALSE)
}

Plot.Theta.Test <- function(){
	\dontrun{
		data(saliva)
		data(throat)
		data(tonsils)
		
		### Combine the data sets into a single list
		group.data <- list(saliva, throat, tonsils)
		
		### Get PI using MLE with CI
		mle <- Est.PI(group.data)$MLE
		
		### Plot theta
		Plot.Theta(mle)
	}
}



### ~~~~~~~~~~~~~~~~~~~~~
### Other functions
### ~~~~~~~~~~~~~~~~~~~~~
C.alpha.multinomial.Test <- function(){
	data(saliva)
	
	calpha <- C.alpha.multinomial(saliva)
	calpha
}

DM.MoM.Test <- function(){
	data(throat)
	
	fit.throat <- DM.MoM(throat)
	fit.throat
}

Kullback.Leibler.Test <- function(){
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Combine the data sets into a single list
	group.data <- list(saliva, throat, tonsils)
	
	\dontrun{
		kl <- Kullback.Leibler(group.data)
		kl
	}
}

Xmcupo.effectsize.Test <- function(){
	data(saliva)
	data(throat)
	
	### Combine the data sets into a single list
	group.data <- list(saliva, throat)
	
	effect <- Xmcupo.effectsize(group.data)
	effect
}

Est.PI.Test <- function(){
	\dontrun{
		data(saliva)
		data(throat)
		data(tonsils)
		
		### Combine the data sets into a single list
		group.data <- list(saliva, throat, tonsils)
		
		### Get PI using MLE and MOM with CI
		piEsts <- Est.PI(group.data)
		
		mle <- piEsts$MLE
		mom <- piEsts$MOM
	}
}

Test.Paired.Test <- function(){
	data(saliva)
	data(throat)
	
	
	### Since saliva and throat come from same subjects, the data is paired 
	saliva1 <- saliva[-24,] # Make saliva 23 subjects to match throat
	group.data <- list(throat, saliva1)
	
	### We use 1 for speed, should be at least 1,000
	numPerms <- 1
	
	pval <- Test.Paired(group.data, numPerms)
	pval
}

DM.Rpart.Test <- function(){
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Create some covariates for our data set
	site <- c(rep("Saliva", nrow(saliva)), rep("Throat", nrow(throat)), 
			rep("Tonsils", nrow(tonsils)))
	covars <- data.frame(Group=site)
	
	### Combine our data into a single object
	data <- rbind(saliva, throat, tonsils)
	
	### For a single rpart tree
	rpartRes <- DM.Rpart(data, covars, numCV=0, numCon=0)
	
	\dontrun{
		### For a cross validated rpart tree
		rpartRes <- DM.Rpart(data, covars, numCV=10, numCon=0)
		
		### For a cross validated rpart tree with consensus
		# Note this is set to 2 for speed and should be at least 100
		rpartRes <- DM.Rpart(data, covars, numCon=2)
	}
}

DM.Rpart.Base.Test <- function(){
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Create some covariates for our data set
	site <- c(rep("Saliva", nrow(saliva)), rep("Throat", nrow(throat)), 
			rep("Tonsils", nrow(tonsils)))
	covars <- data.frame(Group=site)
	
	### Combine our data into a single object
	data <- rbind(saliva, throat, tonsils)
	
	rpartRes <- DM.Rpart.Base(data, covars)
}

DM.Rpart.CV.Test <- function(){
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Create some covariates for our data set
	site <- c(rep("Saliva", nrow(saliva)), rep("Throat", nrow(throat)), 
			rep("Tonsils", nrow(tonsils)))
	covars <- data.frame(Group=site)
	
	### Combine our data into a single object
	data <- rbind(saliva, throat, tonsils)
	
	rpartRes <- DM.Rpart.CV(data, covars)
}

DM.Rpart.CV.Consensus.Test <- function(){
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Create some covariates for our data set
	site <- c(rep("Saliva", nrow(saliva)), rep("Throat", nrow(throat)), 
			rep("Tonsils", nrow(tonsils)))
	covars <- data.frame(Group=site)
	
	### Combine our data into a single object
	data <- rbind(saliva, throat, tonsils)
	
	### This should be at least 100, but we use 2 for speed
	numCon <- 2
	rpartRes <- DM.Rpart.CV.Consensus(data, covars, numCon=numCon)
}

Gen.Alg.Test <- function(){
	\dontrun{
		data(saliva)
		data(throat)
		
		### Combine the data into a single data frame
		group.data <- list(saliva, throat)
		group.data <- formatDataSets(group.data)
		data <- do.call("rbind", group.data)
		
		### Normalize the data by subject
		dataNorm <- t(apply(data, 1, function(x){x/sum(x)}))
		
		### Set covars to just be group membership
		memb <- c(rep(0, nrow(saliva)), rep(1, nrow(throat)))
		covars <- matrix(memb, length(memb), 1)
		
		### We use low numbers for speed. The exact numbers to use depend
		### on the data being used, but generally the higher iters and popSize 
		### the longer it will take to run.  earlyStop is then used to stop the
		### run early if the results aren't improving.
		iters <- 500
		popSize <- 200
		earlyStop <- 250
		
		gaRes <- Gen.Alg(dataNorm, covars, iters, popSize, earlyStop)
	}
}

Gen.Alg.Consensus.Test <- function(){
	\dontrun{
		data(saliva)
		data(throat)
		
		### Combine the data into a single data frame
		group.data <- list(saliva, throat)
		group.data <- formatDataSets(group.data)
		data <- do.call("rbind", group.data)
		
		### Normalize the data by subject
		dataNorm <- t(apply(data, 1, function(x){x/sum(x)}))
		
		### Set covars to just be group membership
		memb <- c(rep(0, nrow(saliva)), rep(1, nrow(throat)))
		covars <- matrix(memb, length(memb), 1)
		
		### We use low numbers for speed. The exact numbers to use depend
		### on the data being used, but generally the higher iters and popSize 
		### the longer it will take to run.  earlyStop is then used to stop the
		### run early if the results aren't improving.
		iters <- 500
		popSize <- 200
		earlyStop <- 250
		numRuns <- 3
		
		gaRes <- Gen.Alg.Consensus(dataNorm, covars, .5, numRuns, FALSE, 3, 
				iters, popSize, earlyStop)
	}
}


### ~~~~~~~~~~~~~~~~~~~~~
### Filter functions
### ~~~~~~~~~~~~~~~~~~~~~
Data.filter.Test <- function(){
	data(saliva) 
	
	### Excludes all samples with fewer than 1000 reads and collapses
	### taxa with 11th or smaller abundance into one category 
	filterDataNum <- Data.filter(saliva, "data", 1000, numTaxa=10) 
	
	### Excludes all samples with fewer than 1000 reads and collapses
	### the least abundant taxa to keep as close to 85% of the data as
	### possible
	filterDataPer <- Data.filter(saliva, "data", 1000, perTaxa=.95) 
	
	dim(saliva)
	dim(filterDataNum)
	dim(filterDataPer)
}

formatDataSets.Test <- function(){
	data(saliva)
	data(throat)
	
	### Set each data set to have 10 different columns
	saliva2 <- saliva[,1:10]
	throat2 <- throat[,11:20]
	
	### Combine the data sets into a single list
	group.data <- list(saliva2, throat2)
	
	formattedData <- formatDataSets(group.data)
	formattedData[[1]][1:5, 1:5]
}



### ~~~~~~~~~~~~~~~~~~~~~
### Generation functions
### ~~~~~~~~~~~~~~~~~~~~~
Dirichlet.multinomial.Test <- function(){
	data(saliva)
	
	### Generate a the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	nrs <- rep(15000, 20) 
	
	### Get gamma from the dirichlet-multinomial parameters
	shape <- dirmult(saliva)$gamma
	
	dmData <- Dirichlet.multinomial(nrs, shape)
	dmData[1:5, 1:5]
}

Multinomial.Test <- function(){
	### Generate the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	nrs <- rep(15000, 25)
	
	### Create a probability vector
	probs <- c(0.4, 0.3, 0.2, .05, 0.04, .01)
	
	mData <- Multinomial(nrs, probs)
	mData[1:5, 1:5]
}



### ~~~~~~~~~~~~~~~~~~~~~
### Sample functions
### ~~~~~~~~~~~~~~~~~~~~~
Xdc.sevsample.Test <- function(){
	data(saliva) 
	data(throat)
	
	### Combine the data sets into a single list
	group.data <- list(saliva, throat)
	
	xdc <- Xdc.sevsample(group.data)
	xdc
}

Xmc.sevsample.Test <- function(){
	data(saliva) 
	data(throat)
	data(tonsils)
	
	### Get pi from the dirichlet-multinomial parameters
	pi0 <- dirmult(saliva)$pi
	
	### Combine the data sets into a single list
	group.data <- list(throat, tonsils)
	
	xmc <- Xmc.sevsample(group.data, pi0)
	xmc
}

Xmcupo.sevsample.Test <- function(){
	data(saliva) 
	data(tonsils)
	data(throat)
	
	### Combine the data sets into a single list
	group.data <- list(saliva, throat, tonsils)
	
	xmcupo <- Xmcupo.sevsample(group.data)
	xmcupo
}

Xoc.sevsample.Test <- function(){
	data(saliva) 
	data(tonsils)
	
	### Combine the data sets into a single list
	group.data <- list(saliva, tonsils)
	
	\dontrun{
		xoc <- Xoc.sevsample(group.data)
		xoc
	}
}

Xsc.onesample.Test <- function(){
	data(saliva)
	data(throat)
	
	### Get pi from the dirichlet-multinomial parameters
	pi0 <- dirmult(saliva)$pi
	
	xsc <- Xsc.onesample(throat, pi0)
	xsc
}
























