\name{Get.Diversities}
\alias{Get.Diversities}
\title{Get Sample Diversity Measures}
\description{Computes the Shannon and Simpson diversities along with the richness.}

\usage{Get.Diversities(data)}

\arguments{
  \item{data}{A matrix of taxonomic counts(columns) for each sample(rows).}
}

\value{A dataframe with different diversity measures.}

\author{Sharina Carter, Elena Deych, Berkley Shands, William D. Shannon}

\examples{
	data(saliva)
	
	salDiv <- Get.Diversities(saliva)
	salDiv
}