\name{Est.PI}
\alias{Est.PI}
\title{Estimate the Pi Vector}
\description{Calculates PI for every group and plots the results.}

\usage{Est.PI(group.data, useMLE = TRUE, main = NULL, plot = TRUE)}

\arguments{
  \item{group.data}{A list of matrices of taxonomic counts(columns) for each sample(rows).}
  \item{useMLE}{A boolean to use the MLE or MoM calculation for PI.}
  \item{main}{A string to be used as the plots title.}
  \item{plot}{A boolean to display the plot or not.}
}

\value{
	A data frame containing the PI, SE and upper/lower bounds for every taxa and
	a plot if 'plot = TRUE'.
}

\author{Sharina Carter, Elena Deych, Berkley Shands, William D. Shannon}

\examples{
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