% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methodRNAseq_findHVG.R
\docType{methods}
\name{findHVG}
\alias{findHVG}
\alias{findHVG,RNAseq-method}
\title{Find highly variable genes}
\usage{
findHVG(object, signThres = 0.1, outputPlots = "./",
  plotting = "pdf", colVarGenes = "deeppink", UseSpike = TRUE,
  minQuantCv2 = 0.2, maxQuantCv2 = 0.8, minQuantMeans = 0.2,
  maxQuantMeans = 0.8)

\S4method{findHVG}{RNAseq}(object, signThres = 0.1, outputPlots = "./",
  plotting = "pdf", colVarGenes = "deeppink", UseSpike = TRUE,
  minQuantCv2 = 0.2, maxQuantCv2 = 0.8, minQuantMeans = 0.2,
  maxQuantMeans = 0.8)
}
\arguments{
\item{object}{\code{RNAseq} object.}

\item{signThres}{Threshold to adjust for multiple testing with the Benjamini-Hochberg 
method. Default: 0.1 (cut at 10 percent).}

\item{outputPlots}{state directory where you want to output the 
plots. Default: Current directory.}

\item{plotting}{It states whether plots are generated or not and you should specify
the type. Options: "pdf", "tiff", "no". Default: "pdf". Default: TRUE}

\item{colVarGenes}{Color for variable genes. Default: deeppink.}

\item{UseSpike}{Logical vector stating whether to use Spike-ins to calculate the fitting
curve or not. Default: NULL (it will depend on whether the SpikeIn slot has spike ins).}

\item{minQuantCv2}{Lower quantile of cv2 values to discard to fit the data. Default: 0.2. This means that those
genes with cv2 (across cells) in the lower 0.2 quantile will be excluded for the fit.}

\item{maxQuantCv2}{Upper quantile of cv2 values to discard to fit the data. Default: 0.8. This means that those
genes with cv2 (across cells) in the top 0.8 quantile will be discarded for the fit.}

\item{minQuantMeans}{Lower quantile of means to discard to fit the data. Default: 0.2. This means that those
genes with mean (across cells) in the 0.2 lower quantile will be discarded for the fit.}

\item{maxQuantMeans}{Upper quantile of means to discard to fit the data. Default: 0.8. This means that those
genes with mean (across cells) in the 0.8 upper quantile will be discarded for the fit.}
}
\value{
vector of highly variable genes. This should be stored in \code{genesHVG}
slot within \code{RNAseq} object.
}
\description{
It finds highly variable genes using the normalised data stored in
the \code{countsNorm} slot within the \code{RNAseq} object provided. It follows
Brennecke et al., Nature Methods, 2013. object@SpikeIn vector is required. 
If any entry is true, it will take spike-ins to calculate the fit. If none is true,
it will take biological genes to calculate the fit.
}
\author{
Blanca Pijuan Sala.
}
