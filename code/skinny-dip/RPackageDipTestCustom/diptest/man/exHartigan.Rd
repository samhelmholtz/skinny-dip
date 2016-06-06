\name{exHartigan}
\alias{exHartigan}
\title{Hartigan's Artificial n-modal Example Data Set}
\description{
 63 (integer) numbers; unimodal or bimodal, that's the question.

 This is now \emph{deprecated}.
 Please use \code{\link{statfaculty}} instead!
}
\examples{
data(exHartigan)
plot(dH <- density(exHartigan))
rug(exHartigan)# should jitter
}
\keyword{data}
