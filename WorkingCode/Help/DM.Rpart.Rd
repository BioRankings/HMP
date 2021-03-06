\name{DM.Rpart}
\alias{DM.Rpart}
\alias{DM.Rpart.Base}
\alias{DM.Rpart.CV}
\alias{DM.Rpart.CV.Consensus}
\title{Dirichlet-Multinomial RPart}
\description{
	This function combines recursive partitioning and the Dirichlet-Multinomial distribution to identify homogeneous 
	subgroups of microbiome taxa count data.
}

\usage{DM.Rpart(data, covars, plot = TRUE, minsplit = 1, minbucket = 1, cp = 0, numCV = 10, 
	numCon = 100, parallel = FALSE, cores = 3, use1SE = FALSE, lowerSE = TRUE)}

\arguments{
  \item{data}{A matrix of taxonomic counts(columns) for each sample(rows).}
  \item{covars}{A matrix of covariates(columns) for each sample(rows).}
  \item{plot}{When 'TRUE' a tree plot of the results will be generated.}
  \item{minsplit}{The minimum number of observations to split on, see \link[rpart]{rpart.control}.}
  \item{minbucket}{The minimum number of observations in any terminal node, see \link[rpart]{rpart.control}.}
  \item{cp}{The complexity parameter, see \link[rpart]{rpart.control}.}
  \item{numCV}{The number folds for a k-fold cross validation. A value less than 2 will return the rpart result without any cross validation.}
  \item{numCon}{The number of cross validations to repeat to achieve a consensus solution.}
  \item{parallel}{When this is 'TRUE' it allows for parallel calculation of consensus. Requires the package \code{doParallel}.}
  \item{cores}{The number of parallel processes to run if parallel is 'TRUE'.}
  \item{use1SE}{See details.}
  \item{lowerSE}{See details.}
}

\value{
	The 3 main things returned are:
	\item{fullTree}{An rpart object without any pruning.}
	\item{bestTree}{A pruned rpart object based on use1SE and lowerSE's settings.}
	\item{cpTable}{Information about the fullTree rpart object and how it splits.}
	
	The other variables returned include surrogate/competing splits, error rates and a plot of the bestTree if plot is TRUE.
}

\details{
	There are 3 ways to run this function.  The first is setting numCV to less than 2, which will run rpart once
	using the DM distribution and the specified minsplit, minbucket and cp.  This result will not have any kind 
	of branch pruning and the objects returned 'fullTree' and 'bestTree' will be the same.
	
	The second way is setting numCV to 2 or greater (we recommend 10) and setting numCon to less than 2.  This will 
	run rpart several times using a k-fold cross validation to prune the tree to its optimal size.  This is the best method to use.
	
	The third way is setting both numCV and numCon to 2 or greater (We recommend at least 100 for numCon).  This will
	repeat the second way numCon times and build a consensus solution.  This method is ONLY needed for low sample sizes.
	
	When the argument 'use1SE' is 'FALSE', the returned object 'bestTree' is the pruned tree with the lowest MSE. 
	When it is 'TRUE', 'bestTree' is either the biggest pruned tree (lowerSE = FALSE) or the smallest pruned tree (lowerSE = TRUE),
	that is within 1 standard error of the lowest MSE.
}

\examples{
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
	numCV <- 0
	numCon <- 0
	rpartRes <- DM.Rpart(data, covars, numCV=numCV, numCon=numCon)
	
	\dontrun{
		### For a cross validated rpart tree
		numCon <- 0
		rpartRes <- DM.Rpart(data, covars, numCon=numCon)
		
		### For a cross validated rpart tree with consensus
		numCon <- 2 # Note this is set to 2 for speed and should be at least 100
		rpartRes <- DM.Rpart(data, covars, numCon=numCon)
	}
}
