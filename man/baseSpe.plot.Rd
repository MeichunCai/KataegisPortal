\name{baseSpe.plot}
\alias{baseSpe.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function baseSpe.plot() to plot C>X mutations base spectum
}
\description{
Plot of flanking sequence of all C>X mutations or C>X mutations within a chromosome/chromosome arm. 
}
\usage{
baseSpe.plot(plot.data, sample = "sample", chr = NULL, arm = NULL, color = NULL, k = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plot.data}{An object as returned by the 'mutSNP.input' function.}
  \item{sample}{Name of sample [default: sample].}
  \item{chr}{Chromosome for which the plot will be generated [default: 1:24].}
  \item{arm}{Chromosome arm for which the plot will be generated. Can be p or q. [default: NULL].}
  \item{color}{Colour to show the bases ("A", "C", "G", "T").}
  \item{k}{'K'-mer bases of flanking sequence of mutated base [default: the same as 'mutSNP.input' function].}
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
baseSpe.plot(plot.data = mutSNP, sample="Test", chr = "chr2",arm="q")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
