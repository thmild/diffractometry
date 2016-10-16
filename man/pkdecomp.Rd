\name{pkdecomp}
\alias{pkdecomp}

\title{ Decomposition of peaks for the whole data set }
\description{
  Calculates decompositions of peaks for the whole diffractogram
}
\usage{
pkdecomp(baslfit,intnum=0, alpha=0.1, maxiter1=500, maxiter=10000, hmax=5, 
maxsolutions=3,heterosk=TRUE,baselim=c(0.05,5),dispers=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baslfit}{Output of \code{baselinefit}}
  \item{intnum}{Vector of numbers of intervals. If \code{intnum = 0}, all intervals are used}
  \item{alpha}{Test level for residual criterion}
  \item{maxiter1}{Number of attempts to fit a model with 1 component}
  \item{maxiter}{Number of attempts to fit a model with k > 1 components}
  \item{hmax}{Maximum number of components}
  \item{maxsolutions}{Number of solutions with k components}
  \item{heterosk}{If \code{TRUE}, the estimate of noise level given in \code{baslfit} is used (default); otherwise noise level is taken to be proportional to signal height}
  \item{baselim}{Limits for changes in the baseline estimate; first component is given in percent of the baseline height, second in counts/2theta}
  \item{dispers}{Additional dispersion factor; not used if \code{heterosk==T}}

}
\details{
Calls \code{\link{pkdecompint}} to decompose the peaks found by \code{\link{baselinefit}} into Pearson Type VII kernels. For every interval, first one kernel is tried. The number of kernels is increased until either a solution accepted by the residual criterion is found or the maximum number of kernels \code{hmax} is reached. After a solution is accepted, for \code{maxsolutions} greater than 1, further decompositions with the same number of kernels can be produced. 
}
\value{
A vector of lists as given by \code{\link{pkdecompint}}
}
\references{P.L. Davies, U. Gather, M. Meise, D. Mergel, T. Mildenberger (2008): 
"Residual based localization and quantification of peaks in x-ray diffractograms", 
Annals of Applied Statistics, Vol. 2, No. 3, 861-886.. 
\url{http://www.statistik.tu-dortmund.de/fileadmin/user_upload/Lehrstuehle/MSind/Publikationen/2008/2008_-_Davies_Gather_Meise_Mergel_Mildenberger_-_Residual_based_localization_and_quantification_of_peaks_in_x-ray_diffractograms.pdf}

T. Bernholt and T. Hofmeister (2006): "An algorithm for a generalized maximum subsequence
problem", in: J. Correa, A. Hevia, M. Kiwi (editors), "Latin 2006: Theoretical
Informatics", volume 3887 of Lecture notes in Computer Science, pages 178-189, Berlin,
Heidelberg. Springer Verlag
}

\author{ T. Mildenberger; Algorithm for residual criterion by T. Bernholt and T. Hofmeister}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{diffractogram}}, \code{\link{baselinefit}}, \code{\link{pkdecompint}}}
\examples{
## Decomposition of data in peak interval into two components

set.seed(0)

par(mfrow=c(2,1))

data(indiumoxide)
indox<-indiumoxide[1901:2400,]
base<-baselinefit(indox)

ind<-c(base$indlsep[1],base$indrsep[1])

plot(indox[ind[1]:ind[2],1], base$baseline$peaks[ind[1]:ind[2]],xlab="",ylab="")

pks<-pkdecomp(base,intnum=1,maxsolutions=1,maxiter1=50,maxiter=250)

lines(indox[ind[1]:ind[2],1],pks[[2]]$fit,col="red")
plot(indox[ind[1]:ind[2],1],pks[[2]]$fitpk[1,],ylim=c(0,1800),type="l",xlab="",ylab="")
lines(indox[ind[1]:ind[2],1],pks[[2]]$fitpk[2,])


}
