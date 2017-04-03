\name{Data.filter}
\alias{Data.filter}
\title{A Data Filter}
\description{
	This function creates a new dataset from an existing one that collapses less-abundant taxa into one category as 
	specified by the user and excludes samples with a total number of reads fewer than the user-specified value.
}

\usage{
	Data.filter(data, order.type = "data", minReads = 0, numTaxa = NULL, 
	perTaxa = NULL, K = NULL, reads.crit = NULL)
}

\arguments{
  \item{data}{A matrix of taxonomic counts(columns) for each sample(rows).}
  \item{order.type}{
		If \code{"sample"}: Rank taxa based on its taxonomic frequency. \cr
		If \code{"data"}: Rank taxa based on cumulative taxonomic counts across all samples.}
  \item{minReads}{Samples with a total number of reads less than read.crit's value will be deleted.}
  \item{numTaxa}{The number of the most abundant taxa to keep, while collapsing the others.}
  \item{perTaxa}{The combined percentage of data to keep, while collapsing the least abundant taxa.}
  \item{K}{Deprecated. Replaced with numTaxa for clarity.}
  \item{reads.crit}{Deprecated. Replaced with minReads for clarity.}
}

\value{
	A data frame of taxa and samples with a total number of reads greater than the minimum value. The last taxon labeled 
	'Other' contains the sum of the least abundant taxa collapsed by setting 'numTaxa' or 'perTaxa'.
}

\author{Patricio S. La Rosa, Elena Deych, Berkley Shands, William D. Shannon}

\examples{
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