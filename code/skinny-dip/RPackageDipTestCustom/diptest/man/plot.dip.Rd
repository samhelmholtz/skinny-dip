\name{plot.dip}
\alias{plot.dip}
\title{Plot a dip() Result, i.e., Class "dip" Object}
\description{
  Plot method for \code{"dip"} objects, i.e., the result of
  \code{\link{dip}(., full.result=TRUE)} or similar.

  Note: We may decide to enhance the plot in the future, possibly not
  entirely back-compatibly.
}
\usage{
\method{plot}{dip}(x, do.points = (n < 20),
     colG = "red3", colL = "blue3", colM = "forest green",
     col.points = par("col"), col.hor = col.points,
     doModal = TRUE, doLegend = TRUE, \dots)
}
\arguments{
  \item{x}{an \R object of \code{\link{class}} \code{"dip"}, i.e.,
    typically the result of \code{\link{dip}(., full.result= FF)} where
    \code{FF} is \code{TRUE} or a string such as \code{"all"}.}
  \item{do.points}{logical indicating if the ECDF plot should include
    points; passed to \code{\link{plot.ecdf}}.}
  \item{colG, colL, colM}{the colors to be used in the graphics for the
    \bold{G}reatest convex minorant, the
    \bold{L}east concave majorant, and the
    \bold{M}odal interval, respectively.}
  \item{col.points, col.hor}{the color of points or horizontal lines,
    respectively, simply passed to \code{\link{plot.ecdf}}.}
  \item{doModal}{logical indicating if the modal interval \eqn{[x_L, x_U]}{[xL, xU]}
    should be shown.}
  \item{doLegend}{logical indicating if a legend should be shown.}
  \item{\dots}{further optional arguments, passed to \code{\link{plot.ecdf}}.}
}
\author{Martin Maechler}
\seealso{
  \code{\link{dip}}, also for examples; \code{\link{plot.ecdf}}.
}
\keyword{hplot}
