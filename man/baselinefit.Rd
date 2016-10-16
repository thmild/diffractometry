\name{baselinefit}
\alias{baselinefit}

\title{ Baseline and peak intervals for diffractometry data }
\description{
  Calculates Taut String and Weighted Smoothing Spline approximations of a diffractogram and identifies the baseline and peak intervals
}
\usage{
baselinefit(data,tau=2.5,gam=1, scl.factor=1.2, maxwdth=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A diffractogram given as a matrix, where the first column gives the angles of diffraction in degrees \eqn{2\theta} and the second column gives the corresponding photon counts. Angles of diffraction are taken to be equidistant.}
  \item{tau}{Value of Parameter tau used in the residual criterion.}
  \item{gam}{Factor which is used in the separation between peaks and baseline.}
  \item{scl.factor}{Factor which is used in the heteroscedastic residual criterion.}
  \item{maxwdth}{Maximum width of the peaks. By default set to 5 degrees.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a list with components
  \item{pmg}{Output (list) of the taut string procedure \code{fnpreg()}. Contains the estimated fit \code{pmg$fn}.}
  \item{spl}{Output (list) of the weighted smoothing spline procedure \code{wsspoiss()}. Contains the estimated fit \code{spl$reg}.}
  \item{baseline}{Output (list) of the baseline fit procedure \code{basiserg()}. Contains the baseline fit \code{baseline$basisl}.}
  \item{npks}{Number of peaks.}
  \item{indlsep}{Left indices of the peak intervals.}
  \item{indrsep}{Right indices of the peak intervals.}
  \item{indextr}{Indices of the extreme values.}
  \item{bs}{Output of the weighted smoothing spline precedure for the data with removed peaks. Contains the baseline for non-peak intervals \code{bs$reg}.}      \item{pks}{A vector which contains the data for the peak intervals and zero otherwise.}
  \item{exb}{Output (list) of the routine \code{exber.maxwdth()} which computes the peak intervals. Contains again indlsep and indrsep.}
  \item{x}{x values of the data.}
  \item{y}{y values of the data.} 
}
\references{P.L. Davies, U. Gather, M. Meise, D. Mergel, T. Mildenberger (2008): 
"Residual based localization and quantification of peaks in x-ray diffractograms", 
Annals of Applied Statistics, Vol. 2, No. 3, 861-886.. 
\url{http://www.statistik.tu-dortmund.de/fileadmin/user_upload/Lehrstuehle/MSind/Publikationen/2008/2008_-_Davies_Gather_Meise_Mergel_Mildenberger_-_Residual_based_localization_and_quantification_of_peaks_in_x-ray_diffractograms.pdf}
}
\author{M. Meise and P.L. Davies}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{diffractogram}}, \code{\link{pkdecomp}} }
\examples{
## baseline fit and peak interval identification for a small segment of the indiumoxide data
## baseline plotted in red and data plotted in peak intervals in blue 

data(indiumoxide)
indox<-indiumoxide[1901:2400,]
plot(indox,xlab="",ylab="")
base<-baselinefit(indox)
lines(indox[,1],base$baseline$basisl,col="red")
points(indox[base$indextr,],col="blue") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
