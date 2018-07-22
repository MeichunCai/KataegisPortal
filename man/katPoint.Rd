\name{katPoint}
\alias{katPoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function katPoint() to identify potential kataegis events
}
\description{
Identify potential kataegis events
}
\usage{
katPoint(data, min.mut = 6, max.dis = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{An object as returned by the 'mutSNP.input' function.}
  \item{min.mut}{Minimum number of consecutive mutations is used to identify potential kataegis events[default: 6].}
  \item{max.dis}{Maximum number of bases between two mutations is used to identify potential kataegis events [default: 1000].}
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
katPoint(mutSNP)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
