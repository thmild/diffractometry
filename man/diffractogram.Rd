\name{diffractogram}
\alias{diffractogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Complete analysis of diffractograms as described in Davies et al. (2008)}
\description{
  Performs a complete analysis of x-ray diffractogram, i.e. calculation of the baseline and the peak intervals as well as decomposition of the peaks.
}
\usage{
diffractogram(data, tau=2.5,gam=1, scl.factor=1.2, maxwdth=5, intnum=0, 
alpha=0.1, maxiter1=500, maxiter=10000, hmax=5, maxsolutions=3,
heterosk=TRUE, baselim=c(0.05,5), dispers=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A diffractogram given as a matrix, where the first column gives the angles of diffraction in \eqn{2\theta} and the second column gives the corresponding photon counts. Angles of diffraction are taken to be equidistant.}
  \item{tau}{Value of Parameter tau used in the residual criterion.}
  \item{gam}{Factor which is used in the separation between peaks and baseline.}
  \item{scl.factor}{Factor which is used in the heteroscedastic residual criterion.}
  \item{maxwdth}{Maximum width of the peaks. By default set to 5 degrees.}
  \item{intnum}{Vector of numbers of intervals. If \code{intnum = 0}, all intervals are used}
  \item{alpha}{Test level for residual criterion}
  \item{maxiter1}{Number of attempts to fit a model with 1 component}
  \item{maxiter}{Number of attempts to fit a model with k > 1 components}
  \item{hmax}{Maximum number of components}
  \item{maxsolutions}{Number of solutions with k components}
  \item{heterosk}{If \code{TRUE}, the estimate of noise level given in \code{baslfit} is used (default); otherwise noise level is taken to be proportional to signal height}
  \item{baselim}{Limits for changes in the baseline estimate; first component is given in percent of the baseline height, second in counts per \eqn{2\theta}}
  \item{dispers}{Additional dispersion factor; not used if \code{heterosk==T}}
  }
\details{
 \code{diffractogram} first calls \code{baselinefit} with the specified parameters and then \code{pkdecomp} for the output of \code{baselinefit}. For further details, see the corresponding help pages.
}
\value{
  Returns a LIST with components
   \item{basl}{Output of \code{baselinefit}}
   \item{pks}{Output of \code{pkdecomp}}
}
\references{P.L. Davies, U. Gather, M. Meise, D. Mergel, T. Mildenberger (2008): 
"Residual based localization and quantification of peaks in x-ray diffractograms", 
Annals of Applied Statistics, Vol. 2, No. 3, 861-886.. 
\url{http://www.statistik.tu-dortmund.de/fileadmin/user_upload/Lehrstuehle/MSind/Publikationen/2008/2008_-_Davies_Gather_Meise_Mergel_Mildenberger_-_Residual_based_localization_and_quantification_of_peaks_in_x-ray_diffractograms.pdf}
}
\author{P.L. Davies, M. Meise, T.Mildenberger}
\seealso{ \code{\link{baselinefit}}, \code{\link{pkdecomp}}}
\examples{

## Complete Analysis of a segment from indiumoxide data:
## Identification of baseline and peak as well as decomposition of the peaks  
## Plot shows baseline and decomposition of first found peak into two components 

set.seed(0)

par(mfrow=c(3,1))

data(indiumoxide)
indox<-indiumoxide[1901:2400,]

sol<-diffractogram(indox,maxsolutions=1,maxiter1=50,maxiter=250)

ind<-c(sol$bas$indlsep[1],sol$bas$indrsep[1])

plot(indox,xlab="",ylab="")
lines(indox[,1],sol$bas$baseline$basisl,col="red")
points(indox[ind[1]:ind[2],],col="blue")


plot(indox[ind[1]:ind[2],1], sol$bas$baseline$peaks[ind[1]:ind[2]],xlab="",ylab="")

lines(indox[ind[1]:ind[2],1],sol$pks[[2]]$fit,col="red")
plot(indox[ind[1]:ind[2],1],sol$pks[[2]]$fitpk[1,],ylim=c(0,1800),type="l",xlab="",ylab="")
lines(indox[ind[1]:ind[2],1],sol$pks[[2]]$fitpk[2,])



}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
