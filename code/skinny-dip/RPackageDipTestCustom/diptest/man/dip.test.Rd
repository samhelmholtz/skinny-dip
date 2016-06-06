\name{dip.test}
\alias{dip.test}
\title{Hartigans' Dip Test for Unimodality}
\description{
  Compute Hartigans' dip statistic \eqn{D_n}{Dn}, and
  its p-value for the test for unimodality, by interpolating
  tabulated quantiles of \eqn{\sqrt{n} D_n}{sqrt(n) * Dn}.

  For \eqn{X_i \sim F, i.i.d.}{X_i ~ F, i.i.d},
  the null hypothesis is that \eqn{F} is a unimodal distribution.
  Consequently, the test alternative is non-unimodal, i.e., at least
  bimodal.  Using the language of medical testing, you would call the
  test \dQuote{Test for \bold{Multi}modality}.
}
\usage{
dip.test(x, simulate.p.value = FALSE, B = 2000)
}
\arguments{
  \item{x}{numeric vector; sample to be tested for unimodality.}
  \item{simulate.p.value}{a logical indicating whether to compute
    p-values by Monte Carlo simulation.}
  \item{B}{an integer specifying the number of replicates used in the
    Monte Carlo test.}
}
\details{
  If \code{simulate.p.value} is \code{FALSE}, the p-value is computed
  via linear interpolation (of \eqn{\sqrt{n} D_n}{sqrt(n) * Dn}) in the
  \code{\link{qDiptab}} table.
  Otherwise the p-value is computed from a Monte Carlo simulation of a
  uniform distribution (\code{\link{runif}(n)}) with \code{B}
  replicates.
}
\value{
  A list with class \code{"htest"} containing the following
  components:
  \item{statistic}{the dip statistic \eqn{D_n}{Dn}, i.e.,
    \code{\link{dip}(x)}.}
  \item{p.value}{the p-value for the test, see details.}
  \item{method}{character string describing the test, and whether Monte
    Carlo simulation was used.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\seealso{
  For goodness-of-fit testing, notably of continuous distributions,
  \code{\link{ks.test}}.
}
\references{
  see those in \code{\link{dip}}.
}
\author{Martin Maechler}
\note{
  see also the package vignette, which describes the procedure in more details.
}
\examples{
## a first non-trivial case
(d.t <- dip.test(c(0,0, 1,1))) # "perfect bi-modal for n=4" --> p-value = 0
stopifnot(d.t$p.value == 0)

data(statfaculty)
plot(density(statfaculty)); rug(statfaculty)
(d.t <- dip.test(statfaculty))

x <- c(rnorm(50), rnorm(50) + 3)
plot(density(x)); rug(x)
## border-line bi-modal ...  BUT (most of the times) not significantly:
dip.test(x)
dip.test(x, simulate=TRUE, B=5000)

## really large n -- get a message
dip.test(runif(4e5))

}
\keyword{htest}
\keyword{distribution}
