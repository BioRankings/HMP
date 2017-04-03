

source("HMP Functions.R")


saliva <- read.csv("saliva.csv", row.names=1)
throat <- read.csv("throat.csv", row.names=1)
tonsils <- read.csv("tonsils.csv", row.names=1)


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
	pval2 <- MC.Xmc.statistics(group.Nrs, numMC, pi0, group.pi, group.theta, "ha")
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
	pval2 <- MC.Xmcupo.statistics(group.Nrs, numMC, pi0, group.pi, group.theta)
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

Plot.Abundance.Test <- function(){
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Combine the data sets into a single list
	group.data <- list(saliva, throat, tonsils)
	
	Plot.Abundance(group.data)
}

Est.PI.Test <- function(){
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Combine the data sets into a single list
	group.data <- list(saliva, throat, tonsils)
	
	\dontrun{
		### Get and plot PI using MLE
		mle <- Est.PI(group.data)
		
		### Get and plot PI using MoM
		mom <- Est.PI(group.data, FALSE)
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

pioest.Test <- function(){
	data(saliva)
	data(throat)
	
	### Combine the data sets into a single list
	group.data <- list(saliva, throat)
	
	pio <- pioest(group.data)
	pio
}

Get.Diversities.Test <- function(){
	data(saliva)
	
	salDiv <- Get.Diversities(saliva)
	salDiv
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



### ~~~~~~~~~~~~~~~~~~~~~
### Filter functions
### ~~~~~~~~~~~~~~~~~~~~~
Data.filter.Test <- function(){
	data(saliva) 
	
	### Excludes all samples with fewer than 1000 reads and collapses
	### taxa with 11th or smaller abundance into one category 
	filterDataNum <- Data.filter(saliva, "data", 1000, 10) 
	
	### Excludes all samples with fewer than 1000 reads and collapses
	### the least abundant taxa to keep as close to 85% of the data as
	### possible
	filterDataPer <- Data.filter(saliva, "data", 1000, perTaxa=.99) 
	
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
























