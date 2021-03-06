library(KataegisPortal)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
mutData <- system.file("examples", "mutData.txt", package="KataegisPortal")
mutData <- read.table(mutData, header = TRUE,sep = "\t",as.is = TRUE)
head(mutData)
mutSNP = mutSNP.input(mut.data = mutData,
			chr = "chr",
			pos = "pos",
			ref = "ref",
			alt = "alt",
			build = "hg19")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
katPoint(mutSNP,txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)
pdf("KataegisPortal.pdf", width=12, height=12)
layout(matrix(c(1,2,3,1,2,3,1,2,3,1,2,4,1,2,4,0,0,4),3,6))
mutDis.plot(plot.data = mutSNP, sample="Test")
mutDis.plot(plot.data = mutSNP, sample="Test", chr = "chr2")
baseSpe.plot(plot.data = mutSNP, sample="Test")
baseSpe.plot(plot.data = mutSNP, sample="Test", chr = "chr2",arm="q")
dev.off()
