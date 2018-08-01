\name{baseSpe.plot}
\alias{baseSpe.plot}
\title{
function baseSpe.plot() to plot C>X mutations base spectum
}
\description{
Plot of flanking sequence of all C>X mutations or C>X mutations within a chromosome/chromosome arm. 
}
\usage{
baseSpe.plot(plot.data, sample = "sample", chr = NULL, arm = NULL, color = NULL, k = NULL)
}
\arguments{
  \item{plot.data}{An object as returned by the 'mutSNP.input' function.}
  \item{sample}{Name of sample [default: sample].}
  \item{chr}{Chromosome for which the plot will be generated [default: 1:24].}
  \item{arm}{Chromosome arm for which the plot will be generated. Can be p or q. [default: NULL].}
  \item{color}{Colour to show the bases ("A", "C", "G", "T").}
  \item{k}{'K'-mer bases of flanking sequence of mutated base [default: the same as 'mutSNP.input' function].}
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
baseSpe.plot(plot.data = mutSNP, sample="Test", chr = "chr2",arm="q")
}
