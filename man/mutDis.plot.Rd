\name{mutDis.plot}
\alias{mutDis.plot}
\title{
function mutDis.plot() to plot intermutation distance
}
\description{
Plot intermutation distance to display kataegis or hyper mutated genomic regions.
}
\usage{
mutDis.plot(plot.data, sample = "sample", chr = NULL, color = NULL, min.mut = 6, max.dis = 1000)
}
\arguments{
  \item{plot.data}{An object as returned by the 'mutSNP.input' function.}
  \item{sample}{Name of sample [default: sample].}
  \item{chr}{Chromosome for which the plot will be generated [default: 1:24].}
  \item{color}{Colour to show the mutations ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G").}
  \item{min.mut}{Minimum number of consecutive mutations is used to identify potential kataegis events [default: 6].}
  \item{max.dis}{Maximum number of bases between two mutations is used to identify potential kataegis events [default: 1000].}
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
mutDis.plot(plot.data = mutSNP, sample="Test", chr = c(1:23))
}
