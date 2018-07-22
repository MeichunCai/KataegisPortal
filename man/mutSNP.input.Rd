\name{mutSNP.input}
\alias{mutSNP.input}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function mutSNP.input() for data transformation
}
\description{
Makes the contexts useable in katPoint(), mutDis.plot() and baseSpe.plot().
}
\usage{
mutSNP.input(mut.data, chr = "chr", pos = "pos", ref = "ref", alt = "alt", build = NULL, k = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mut.data}{A data frame of mutation(SNV) file.}
  \item{chr}{Column name in the mutation file corresponding to the chromosome.}
  \item{pos}{Column name in the mutation file corresponding to the mutation position.}
  \item{ref}{Column name in the mutation file corresponding to the reference base.}
  \item{alt}{Column name in the mutation file corresponding to the alternate base.}
  \item{build}{Reference build for chromosome sizes. Can be hg18, hg19 or hg38.}
  \item{k}{'K'-mer bases of flanking sequence of mutated base [default: 10].}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
