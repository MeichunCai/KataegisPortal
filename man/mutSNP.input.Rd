\name{mutSNP.input}
\alias{mutSNP.input}
\title{
function mutSNP.input() for data transformation
}
\description{
Makes the contexts useable in katPoint(), mutDis.plot() and baseSpe.plot().
}
\usage{
mutSNP.input(mut.data, chr = "chr", pos = "pos", ref = "ref", alt = "alt", build = NULL, k = 10)
}
\arguments{
  \item{mut.data}{A data frame of mutation(SNV) file.}
  \item{chr}{Column name in the mutation file corresponding to the chromosome.}
  \item{pos}{Column name in the mutation file corresponding to the mutation position.}
  \item{ref}{Column name in the mutation file corresponding to the reference base.}
  \item{alt}{Column name in the mutation file corresponding to the alternate base.}
  \item{build}{Reference build for chromosome sizes. Can be hg18, hg19 or hg38.}
  \item{k}{'K'-mer bases of flanking sequence of mutated base [default: 10].}
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
}
