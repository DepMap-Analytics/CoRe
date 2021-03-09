\name{CoRe.download_AnnotationModel}
\alias{CoRe.download_AnnotationModel}
\title{Download Cell Passport models annotation file}
\usage{CoRe.download_AnnotationModel(
URL='https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz')}
\arguments{
\item{URL}{The URL specifying the online location of the annotation file. Default value points to the most up-to-date version of the Cell Model Passports annotation file (which is kept updated).}
}
\value{
A data frame with one row per model and one colum per annotation entry.
}
\description{Downloading Cell Model Passport annotation file [1].}
\examples{
## Downloading Cell Model Passport annotation file
modelAnn<-CoRe.download_AnnotationModel()

head(modelAnn)
}
\references{
[1] Van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929.
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
