katPoint <- function(data, sample = "sample", min.mut = 6, max.dis = 1000, 
    txdb = NULL) {
    build = data$build[1]
    genome.opts = c("hg19", "hg18", "hg38")
    if (!build %in% genome.opts) {
        stop("Available reference builds: hg18, hg19, hg38")
    }
    if (build == "hg19") {
        chr.arm = c(1.25e+08, 93300000, 9.1e+07, 50400000, 48400000, 
            6.1e+07, 59900000, 45600000, 4.9e+07, 40200000, 53700000, 
            35800000, 17900000, 17600000, 1.9e+07, 36600000, 2.4e+07, 
            17200000, 26500000, 27500000, 13200000, 14700000, 60600000, 
            12500000)
    } else if (build == "hg18") {
        chr.arm = c(124300000, 93300000, 91700000, 50700000, 47700000, 
            60500000, 59100000, 45200000, 51800000, 40300000, 52900000, 
            35400000, 1.6e+07, 15600000, 1.7e+07, 38200000, 22200000, 
            16100000, 28500000, 27100000, 12300000, 11800000, 59500000, 
            11300000)
    } else if (build == "hg38") {
        chr.arm = c(123400000, 93900000, 90900000, 5e+07, 48800000, 
            59800000, 60100000, 45200000, 4.3e+07, 39800000, 53400000, 
            35500000, 17700000, 17200000, 1.9e+07, 36800000, 25100000, 
            18500000, 26200000, 28100000, 1.2e+07, 1.5e+07, 6.1e+07, 
            10400000)
    } else {
        stop("Available reference builds: hg18, hg19, hg38")
    }
    num = dim(data)[1] - 5
    katPoint <- matrix(nrow = num, ncol = 8)
    i = 1
    mutnum = 1
    Cmutnum = 0
    for (i in 1:num) {
	if (data$ref[i] %in% c("C", "G")){
            Cmutnum = Cmutnum + 1
	}
        if (data$dis[i + 1] <= max.dis) {
            mutnum = mutnum + 1
        } else {
            if (mutnum >= min.mut) {
                len = data$pos[i] - data$pos[i - mutnum + 1] + 1
                chr.n = gsub(pattern = "chr", replacement = "", x = data$chr[i], 
                  fixed = TRUE)
                chr.n = gsub(pattern = "X", replacement = "23", x = chr.n, 
                  fixed = TRUE)
                chr.n = gsub(pattern = "Y", replacement = "24", x = chr.n, 
                  fixed = TRUE)
                chr.n = as.numeric(chr.n)
                if (data$pos[i] <= chr.arm[chr.n]) {
                  arm = paste(chr.n, "p", sep = "")
                } else if (data$pos[i - mutnum + 1] >= chr.arm[chr.n]) {
                  arm = paste(chr.n, "q", sep = "")
                } else {
                  arm = paste(chr.n, "p, ", chr.n, "q", sep = "")
                }
                katPoint[i, 1:8] = c(sample, data$chr[i], data$pos[i - mutnum + 
                  1], data$pos[i], arm, len, mutnum, round(Cmutnum/mutnum,3))
            }
            mutnum = 1
	    Cmutnum = 0
        }
    }
    katPoint.out = data.frame(na.omit(katPoint))
    names(katPoint.out) = c("sample", "chrom", "start", "end", "chrom.arm", "length", "number.mut", 
        "weight.C>X")
    for (i in 1:dim(katPoint.out)[1]) {
        if (as.numeric(as.character(katPoint.out$"weight.C>X"[i])) < 0.8) {
            katPoint.out$confidence[i] = 0
	} else {
            katPoint.out$confidence[i] <- length(which(subset(katPoint.out,
               as.numeric(as.character(katPoint.out$"weight.C>X")) >= 0.8)$chrom == katPoint.out$chrom[i]))
            if (katPoint.out$confidence[i] > 3) {
                katPoint.out$confidence[i] = 3
            }
	}
    }
    if (!is.null(txdb)) {
        gr <- GRanges(seqnames = Rle(katPoint.out$chrom), ranges=IRanges(start = 
            as.numeric(as.character(katPoint.out$start)), end =as.numeric(as.character(katPoint.out$end))))
        peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
        katPoint.out$annotation <- peakAnno@anno$annotation
        katPoint.out$distanceToTSS <- peakAnno@anno$distanceToTSS
        katPoint.out$geneName <- peakAnno@anno$SYMBOL
        katPoint.out$geneID <- peakAnno@anno$geneId
    } 
    message(paste(dim(katPoint.out)[1], "potential kataegis events identified", 
        sep = " "))
    return(katPoint.out)
}
