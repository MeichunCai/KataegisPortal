baseSpe.plot <- function(plot.data, sample = "sample", chr = NULL, 
    arm = NULL, color = NULL, k = NULL) {
    build = plot.data$build[1]
    genome.opts = c("hg19", "hg18", "hg38","CHM13")
    if (!build %in% genome.opts) {
        stop("Available reference builds: hg18, hg19, hg38","CHM13")
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
    } else if (build == "CHM13") {
            chr.lens = c(248387328,242696752,201105948,193574945,182045439,
            172126628,160567428,146259331,150617247,134758134,135127769,
            133324548,113566686,101161492,99753195,96330374,84276897,
            80542538,61707364,
            66210255,45090682,51324926,154259566,62460029)
        bsg = BSgenome.Hsapiens.CHM13
    } else {
        stop("Available reference builds: hg18, hg19, hg38","CHM13")
    }
    if (is.null(color)) {
        col = c("darkgreen", "darkblue", "grey", "darkred")
    } else {
        col = color
    }
    names(col) = c("A", "C", "G", "T")
    plot.data = plot.data[which(plot.data$ref %in% c("C", "G")), ]
    if (is.null(chr)) {
        seq = c(1:24)
        seq0 = "whole genome"
    } else {
        plot.data$pos.updated = plot.data$pos
        seq0 = gsub(pattern = "chr", replacement = "", x = chr, fixed = TRUE)
        seq = gsub(pattern = "X", replacement = "23", x = seq0, fixed = TRUE)
        seq = gsub(pattern = "Y", replacement = "24", x = seq, fixed = TRUE)
        seq = as.numeric(seq)
        plot.data = plot.data[which(plot.data$seq %in% seq), ]
        if (!is.null(arm)) {
            if (arm == "p") {
                plot.data = plot.data[which(plot.data$pos <= chr.arm[seq]), 
                  ]
            } else {
                plot.data = plot.data[which(plot.data$pos > chr.arm[seq]), 
                  ]
            }
        }
    }
    if (is.null(k)) {
        k = (nchar(as.character(plot.data$context[1])) - 1)/2
    } else {
        build = plot.data$build[1]
        genome.opts = c("hg19", "hg18", "hg38","CHM13")
        if (!build %in% genome.opts) {
            stop("Available reference builds: hg18, hg19, hg38","CHM13")
        }
        if (build == "hg19") {
            bsg = BSgenome.Hsapiens.UCSC.hg19
        } else if (build == "hg18") {
            bsg = BSgenome.Hsapiens.UCSC.hg18
        } else if (build == "hg38") {
            bsg = BSgenome.Hsapiens.UCSC.hg38
        } else if (build == "CHM13"){
            bsg = BSgenome.Hsapiens.CHM13
        } else {
            stop("Available reference builds: hg18, hg19, hg38","CHM13")
        }
        conv.start = plot.data$pos - k
        conv.end = plot.data$pos + k
        context = getSeq(bsg, plot.data$chr, start = conv.start, end = conv.end)
        if (TRUE) {
            idx = DNAStringSet(plot.data$ref) %in% c("A", "G")
            context[idx] = reverseComplement(context[idx])
        }
        plot.data$context = context
    }
    base <- rbind(data.frame(strsplit(as.character(plot.data$context), 
        "")))
    n = 2 * k + 1
    nBase <- matrix(0, nrow = 4, ncol = n, dimnames = list(c("A", "C", 
        "G", "T"), c(-k:k)))
    for (i in 1:n) {
        nBase[1, i] = length(base[i, ][base[i, ] == "A"])
        nBase[2, i] = length(base[i, ][base[i, ] == "C"])
        nBase[3, i] = length(base[i, ][base[i, ] == "G"])
        nBase[4, i] = length(base[i, ][base[i, ] == "T"])
    }
    par(mai = c(1, 1, 1, 1.5))
    if (seq0 == "whole genome") {
        barplot(nBase, col = col, space = 0, yaxt = "n", main = paste("C>X mutations in", 
            sample, "on whole genome", sep = " "), xlab = "Flanking bases", 
            ylab = "Number of bases")
    } else {
        barplot(nBase, col = col, space = 0, yaxt = "n", main = paste("C>X mutations in", 
            sample, "on chr", seq0, arm, sep = " "), xlab = "Flanking bases", 
            ylab = "Number of bases")
    }
    box()
    axis(2, las = 1, lwd.tick = 0.5, mgp = c(2, 1, 0))
    legend("right", names(col), col = col, pch = 15, inset = c(-0.12, 
        0), bty = "n", xpd = TRUE)
}
