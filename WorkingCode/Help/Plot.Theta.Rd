\name{Plot.Theta}
\alias{Plot.Theta}
\title{Plot the Theta Value}
\description{Plots the theta value with CI for every group.}

\usage{
	Plot.Theta(estPi, main = "PI Vector")
}

\arguments{
  \item{estPi}{The results for either MLE or MOM from the function 'Est.Pi'.}
  \item{main}{A string to be used as the plots title.}
}

\value{A plot of the theta values for every group.}

\examples{
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