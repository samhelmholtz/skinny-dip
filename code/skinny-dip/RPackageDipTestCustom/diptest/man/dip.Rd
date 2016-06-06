\name{dip}
\alias{dip}
\title{Compute Hartigans' Dip Test Statistic for Unimodality}
\description{
  Computes Hartigans' dip test statistic for testing unimodality,
  and additionally the modal interval.
}
\usage{
dip(x, full.result = FALSE, min.is.0 = FALSE, debug = FALSE)
}
\arguments{
  \item{x}{numeric; the data.}
  \item{full.result}{logical or string; \code{dip(., full.result=TRUE)} returns the full result
    list; if \code{"all"} it additionally uses the \code{mn} and
    \code{mj} components to compute the initial GCM and LCM, see below.}
  \item{min.is.0}{logical indicating if the \bold{min}imal value of the
    dip statistic \eqn{D_n}{Dn} can be zero or not.  Arguably should be
    set to \code{TRUE} for internal consistency reasons, but is false by
    default both for continuity and backwards compatibility reasons, see
    the examples below.}
    % backcompatibility both with earlier
    % versions of the \pkg{diptest} package, and with Hartigan's original
    % implementation.}
  \item{debug}{logical; if true, some tracing information is printed
    (from the C routine).}
}
\value{
  depending on \code{full.result} either a number, the dip statistic, or
  an object of class \code{"dip"} which is a \code{\link{list}} with components
  \item{x}{the sorted \code{\link{unname}()}d data.}
  \item{n}{\code{length(x)}.}
  \item{dip}{the dip statistic}
  \item{lo.hi}{indices into \code{x} for lower and higher end of modal interval}
  \item{xl, xu}{lower and upper end of modal interval}
  \item{gcm, lcm}{(last used) indices for \bold{g}reatest \bold{c}onvex
    \bold{m}inorant and the \bold{l}east \bold{c}oncave \bold{m}ajorant.}
  \item{mn, mj}{index vectors of length \code{n} for the GC minorant and
    the LC majorant respectively.}

  For \dQuote{full} results of class \code{"dip"}, there are
  \code{\link{print}} and \code{\link{plot}} methods, the latter with
  its own \link[=plot.dip]{manual page}.
}
\note{
  For \eqn{n \le 3}{n <= 3} where \code{n <- length(x)}, the dip
  statistic \eqn{D_n}{Dn} is always the same minimum value,
  \eqn{1/(2n)}, i.e., there's no possible dip test.
  Note that up to May 2011, from Hartigan's original Fortran code, \code{Dn}
  was set to zero, when all \code{x} values were identical.  However,
  this entailed discontinuous behavior, where for arbitrarily close
  data \eqn{\tilde x}{x~}, \eqn{D_n(\tilde x) = \frac 1{2n}}{Dn(x~) = 1/(2n)}.

  Yong Lu \email{lyongu+@cs.cmu.edu} found in Oct 2003 that the code
  was not giving symmetric results for mirrored data (and was giving
  results of almost 1, and then found the reason, a misplaced \samp{")"}
  in the original Fortran code.  This bug has been corrected for diptest
  version 0.25-0.

  Nick Cox (Durham Univ.) said (on March 20, 2008 on the Stata-list):\cr
  As it comes from a bimodal husband-wife collaboration, the name
  perhaps should be \emph{\dQuote{Hartigan-Hartigan dip test}}, but that
  does not seem to have caught on.  Some of my less statistical
  colleagues would sniff out the hegemony of patriarchy there, although
  which Hartigan is being overlooked is not clear.

  Martin Maechler, as a Swiss, and politician, would say:\cr
  Let's find a compromise, and call it \emph{\dQuote{Hartigans' dip test}},
  so we only have to adapt orthography (:-).
}
\references{
  P. M. Hartigan (1985)
  Computation of the Dip Statistic to Test for Unimodality;
  \emph{Applied Statistics (JRSS C)} \bold{34}, 320--325.\cr
  Corresponding (buggy!) Fortran code of \sQuote{AS 217} available from Statlib,
  \url{http://lib.stat.cmu.edu/apstat/217}

  J. A. Hartigan and P. M. Hartigan (1985)
  The Dip Test of Unimodality;
  \emph{Annals of Statistics} \bold{13}, 70--84.
}
\author{Martin Maechler \email{maechler@stat.math.ethz.ch}, 1994,
  based on S (S-PLUS) and C code donated from Dario Ringach
  \email{dario@wotan.cns.nyu.edu} who had applied \command{f2c} on the
  original Fortran code available from Statlib.

  In Aug.1993, recreated and improved Hartigans' "P-value" table, which
  later became \code{\link{qDiptab}}.
}
\seealso{
  \code{\link{dip.test}} to compute the dip \emph{and} perform the unimodality test,
  based on P-values, interpolated from \code{\link{qDiptab}};
  \code{\link{isoreg}} for isotonic regression.
}
\examples{
data(statfaculty)
plot(density(statfaculty))
rug(statfaculty, col="midnight blue"); abline(h=0, col="gray")
dip(statfaculty)
(dS <- dip(statfaculty, full = TRUE, debug = TRUE))
plot(dS)
## even more output -- + plot showing "global" GCM/LCM:
(dS2 <- dip(statfaculty, full = "all", debug = 3))
plot(dS2)

data(faithful)
fE <- faithful$eruptions
plot(density(fE))
rug(fE, col="midnight blue"); abline(h=0, col="gray")
dip(fE, debug = 2) ## showing internal work
(dE <- dip(fE, full = TRUE)) ## note the print method
plot(dE, do.points=FALSE)

data(precip)
plot(density(precip))
rug(precip, col="midnight blue"); abline(h=0, col="gray")
str(dip(precip, full = TRUE, debug = TRUE))

##-----------------  The  'min.is.0' option :  ---------------------

##' dip(.) continuity and 'min.is.0' exploration:
dd <- function(x, debug=FALSE) {
   x_ <- x ; x_[1] <- 0.9999999999 * x[1]
   rbind(dip(x , debug=debug),
         dip(x_, debug=debug),
         dip(x , min.is.0=TRUE, debug=debug),
         dip(x_, min.is.0=TRUE, debug=debug), deparse.level=2)
}

dd( rep(1, 8) ) # the 3rd one differs ==> min.is.0=TRUE is *dis*continuous
dd( 1:7 )       # ditto

dd( 1:7, debug=TRUE)
## border-line case ..
dd( 1:2, debug=TRUE)

## Demonstrate that  'min.is.0 = TRUE'  does not change the typical result:
B.sim <- 1000 # or larger
D5  <- {set.seed(1); replicate(B.sim, dip(runif(5)))}
D5. <- {set.seed(1); replicate(B.sim, dip(runif(5), min.is.0=TRUE))}
stopifnot(identical(D5, D5.), all.equal(min(D5), 1/(2*5)))
hist(D5, 64); rug(D5)

D8  <- {set.seed(7); replicate(B.sim, dip(runif(8)))}
D8. <- {set.seed(7); replicate(B.sim, dip(runif(8), min.is.0=TRUE))}
stopifnot(identical(D8, D8.))
}
\keyword{htest}
\keyword{distribution}
