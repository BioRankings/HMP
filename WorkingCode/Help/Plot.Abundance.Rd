\name{Plot.Abundance}
\alias{Plot.Abundance}
\title{Plot Abundance of Microbiome Samples}
\description{Plots the mean, max and min for every taxa for every sample.}

\usage{Plot.Abundance(group.data, main="Group Abundance")}

\arguments{
  \item{group.data}{A list of matrices of taxonomic counts(columns) for each sample(rows).}
  \item{main}{A string to be used as the plots title.}
}

\value{A plot which shows the mean, max and min for every data set given.}

\author{Sharina Carter, Elena Deych, Berkley Shands, William D. Shannon}

\examples{
	data(saliva)
	data(throat)
	data(tonsils)
	
	### Combine the data sets into a single list
	group.data <- list(saliva, throat, tonsils)
	
	Plot.Abundance(group.data)
}