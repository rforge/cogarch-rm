\name{COGgetData}
\alias{COGgetData}
\alias{COGARCH:Data}
\alias{ICA-COGARCH:DATA}
\alias{Risk Measure Cogarch data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Download Data for analysis from a specified data provider}
\description{ The function returns an object of class \code{\link{yuima.data-class}} that contains downloaded data from a specified data provider (default \code{yahoo}). Available data provider are.}
\usage{COGgetData(Ticker, from = NULL, to = NULL, src = "yahoo", type = "Close", delta = NULL, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ticker}{ inserire}
  \item{from}{ inserire}
  \item{to}{ inserire}
  \item{src}{ inserire}
  \item{type}{ inserire}
  \item{delta}{ inserire}
  \item{\dots}{ inserire}
}
\details{Si puo togliere}
\value{The object is an object of yuima.data}
\references{Inserire qualche manuale per i dati o niente}
\author{Nome}
\note{Puoi togliere}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
Ticker <- c("AAL.L", "ABF.L", "ADM.L")
from <- "2011-01-01"
to <- "2017-01-04"
X <- COGgetData(Ticker = Ticker, from = from, to = to)
plot(X)
}
