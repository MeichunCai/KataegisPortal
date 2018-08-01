\name{katPoint}
\alias{katPoint}
\title{
function katPoint() to identify potential kataegis events
}
\description{
Identify potential kataegis events
}
\usage{
katPoint(data, sample = "sample", min.mut = 6, max.dis = 1000, txdb = NULL)
}
\arguments{
  \item{data}{An object as returned by the 'mutSNP.input' function.}
  \item{sample}{Name of sample [default: sample].}
  \item{min.mut}{Minimum number of consecutive mutations is used to identify potential kataegis events[default: 6].}
  \item{max.dis}{Maximum number of bases between two mutations is used to identify potential kataegis events [default: 1000].}
  \item{txdb}{An annotation databases, e.g. TxDb.Hsapiens.UCSC.hg19.knownGene [default: NULL].}
}
\examples{
library(KataegisPortal)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
mutData <- "examples/mutData.txt"
mutData <- read.table(mutData, header = TRUE,sep = "\t",as.is = TRUE)
mutSNP = mutSNP.input(mut.data = mutData,
			chr = "chr",
			pos = "pos",
			ref = "ref",
			alt = "alt",
			build = "hg19")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
katPoint(mutSNP,txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)
}

