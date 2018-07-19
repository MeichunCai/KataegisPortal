library(KataegisPortal)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
mutData <- "examples/mutData.txt"
mutData <- read.table(mutData, header = TRUE,sep = "\t",as.is = TRUE)
head(mutData)
mutSNP = mutSNP.input(mut.data = mutData,
			chr = "chr",
			pos = "pos",
			ref = "ref",
			alt = "alt",
			build = "hg19",
			k=20)
pdf("examples/KataegisPortal.pdf", width=12, height=8)
layout(matrix(c(1,2,1,2,1,2,1,3,1,3,0,3),2,6))
mutDis.plot(plot.data = mutSNP, sample="Test", chr = c(1:23))
baseSpe.plot(plot.data = mutSNP, sample="Test")
baseSpe.plot(plot.data = mutSNP, sample="Test", chr = "chr2",arm="q")
dev.off()
