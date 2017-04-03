\name{Kullback.Leibler}
\alias{Kullback.Leibler}
\title{Kullback Leibler}
\description{Calculates Kullback Leibler divergence for all pairwise combinations of the datasets given.}

\usage{Kullback.Leibler(group.data, plot = TRUE, parallel = FALSE, cores = 3)}

\arguments{
  \item{group.data}{A list where each element is a matrix of taxonomic counts(columns) for each sample(rows).}
  \item{plot}{When 'TRUE' a heatmap of the results will also be generated.}
  \item{parallel}{When this is 'TRUE' it allows for parallel calculation of the bootstraps. (See Note 1 in Details).}
  \item{cores}{The number of parallel processes to run if parallel is 'TRUE'.}
}

\details{
	Notes: 
	\enumerate{
		\item Enabling parallel calculation requires the package \code{doParallel}.
	}
}

\value{A matrix of Kullback Leibler divergences and a heatmap if plot is TRUE.}

\author{Sharina Carter, Elena Deych, Berkley Shands, William D. Shannon}

\examples{
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
