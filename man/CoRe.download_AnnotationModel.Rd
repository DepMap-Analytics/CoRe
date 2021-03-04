\name{CoRe.download_AnnotationModel}
\alias{CoRe.download_AnnotationModel}
\title{Download Cell Passport models annotation file}
\usage{CoRe.download_AnnotationModel(
URL='https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz')}
\arguments{
\item{URL}{url link to download Cell Model Passport annotation file.}
}
\value{
\item{X}{Cell Model Passport annotation file.}
}
\description{Downloading Cell Model Passport annotation file [1].}
\details{
If URL link exists download Binary Dependency Matrix.
}
\examples{
## Downloading Cell Model Passport annotation file
modelAnn<-CoRe.download_AnnotationModel('https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz')

head(modelAnn)
}
\references{
[1] Van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929.
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
