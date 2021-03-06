% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phencorr.R
\name{pheno.corr}
\alias{pheno.corr}
\title{Phenotypic Correlation Analysis}
\usage{
pheno.corr(data, genotypes, replication)
}
\arguments{
\item{data}{traits to be analyzed}

\item{genotypes}{vector containing genotypes}

\item{replication}{vector containing replications}
}
\value{
Phenotypic correlation
}
\description{
Phenotypic Correlation Analysis
}
\examples{
data(vardata)
pheno.corr(vardata[3:11],vardata$Genotypes,vardata$Replication)
}
