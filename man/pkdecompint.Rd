\name{pkdecompint}
\alias{pkdecompint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Decomposition of peaks in an interval }
\description{
  Decomposition of peaks in an interval of the diffractogram
}
\usage{
pkdecompint(baslfit, intnum, k, thresh=0, alpha=0.1, heterosk=TRUE,
maxiter=10000, dispers=1, baselim=c(0.05,5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baslfit}{Output of baslfit}
  \item{intnum}{Number of interval}
  \item{k}{Number of peak components to fit}
  \item{thresh}{Threshold for residual criterion}
  \item{alpha}{Test level for residual criterion}
  \item{heterosk}{If \code{TRUE}, the estimate of noise level given in \code{baslfit} is used (default); otherwise noise level is taken to be proportional to signal height}
  \item{maxiter}{Number of attempts to fit a model with k components}
  \item{dispers}{Additional dispersion factor; not used if \code{heterosk==T}}
  \item{baselim}{Limits for changes in the baseline estimate; first component is given in percent of the baseline height, second in counts/2theta}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns a LIST with components
  \item{intnr}{Number of interval }
  \item{x}{values of 2theta}
  \item{y}{the diffractogram with baseline removed}
  \item{fit}{the resulting fit, evaluated at all points of \code{x}}
  \item{fitpk}{a matrix with \code{num.ker} rows that contain fits of the individual peak components}
  \item{basl}{the basline estimate as given in \code{baslfit} }  
  \item{baslchg}{chnage of baseline estimate}
  \item{rss}{residual sum of squares, standardized by noise level estimate}
  \item{num.ker}{number of peak components}
  \item{par}{parameter vector as given in section 8 of Davies et al. (2008)}
  \item{parbl}{intercept and slope of the baseline change}
  \item{parpks}{physical characteristics of the peaks}	
  \item{accept}{is the fit accepted by the residual criterion}
  \item{alpha}{test level for residual criterion}
  \item{thresh}{threshold used in residual criterion}
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
\seealso{\code{\link{diffractogram}}, \code{\link{baselinefit}}, \code{\link{pkdecomp}}}
\examples{
## Decomposition of data in peak interval into two components

set.seed(0)

par(mfrow=c(2,1))

data(indiumoxide)
indox<-indiumoxide[1901:2400,]
base<-baselinefit(indox)

ind<-c(base$indlsep[1],base$indrsep[1])

plot(indox[ind[1]:ind[2],1], 
base$baseline$peaks[ind[1]:ind[2]],xlab="",ylab="")

pks<-pkdecompint(base,intnum=1,k=2)

lines(indox[ind[1]:ind[2],1],pks$fit,col="red")
plot(indox[ind[1]:ind[2],1],pks$fitpk[1,],ylim=c(0,1800),type="l",xlab="",ylab="")
lines(indox[ind[1]:ind[2],1],pks$fitpk[2,])

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
