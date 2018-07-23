mutDis.plot <- function(plot.data, sample = "sample", chr = NULL, color = NULL, 
    min.mut = 6, max.dis = 1000) {
    if (is.null(color)) {
        col = RColorBrewer::brewer.pal(n = 6, name = "Set1")
    } else {
        col = color
    }
    names(col) = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    build = plot.data$build[1]
    genome.opts = c("hg19", "hg18", "hg38")
    if (!build %in% genome.opts) {
        stop("Available reference builds: hg18, hg19, hg38")
    }
    if (build == "hg19") {
        chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 
            171115067, 159138663, 146364022, 141213431, 135534747, 
            135006516, 133851895, 115169878, 107349540, 102531392, 
            90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 
            51304566, 155270560, 59373566)
        chr.arm = c(1.25e+08, 93300000, 9.1e+07, 50400000, 48400000, 
            6.1e+07, 59900000, 45600000, 4.9e+07, 40200000, 53700000, 
            35800000, 17900000, 17600000, 1.9e+07, 36600000, 2.4e+07, 
            17200000, 26500000, 27500000, 13200000, 14700000, 60600000, 
            12500000)
    } else if (build == "hg18") {
        chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 
            170899992, 158821424, 146274826, 140273252, 135374737, 
            134452384, 132349534, 114142980, 106368585, 100338915, 
            88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 
            49691432, 154913754, 57772954)
        chr.arm = c(124300000, 93300000, 91700000, 50700000, 47700000, 
            60500000, 59100000, 45200000, 51800000, 40300000, 52900000, 
            35400000, 1.6e+07, 15600000, 1.7e+07, 38200000, 22200000, 
            16100000, 28500000, 27100000, 12300000, 11800000, 59500000, 
            11300000)
    } else if (build == "hg38") {
        chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 
            170805979, 159345973, 145138636, 138394717, 133797422, 
            135086622, 133275309, 114364328, 107043718, 101991189, 
            90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 
            50818468, 156040895, 57227415)
        chr.arm = c(123400000, 93900000, 90900000, 5e+07, 48800000, 
            59800000, 60100000, 45200000, 4.3e+07, 39800000, 53400000, 
            35500000, 17700000, 17200000, 1.9e+07, 36800000, 25100000, 
            18500000, 26200000, 28100000, 1.2e+07, 1.5e+07, 6.1e+07, 
            10400000)
    } else {
        stop("Available reference builds: hg18, hg19, hg38")
    }
    if (is.null(chr)) {
        seq = c(1:24)
        seq0 = c(1:22, "X", "Y")
    } else {
        plot.data$pos.updated = plot.data$pos
        seq0 = gsub(pattern = "chr", replacement = "", x = chr, fixed = TRUE)
        seq = gsub(pattern = "X", replacement = "23", x = seq0, fixed = TRUE)
        seq = gsub(pattern = "Y", replacement = "24", x = seq, fixed = TRUE)
        seq = as.numeric(seq)
        plot.data = plot.data[which(plot.data$seq %in% seq), ]
        for (i in 1:length(chr.lens)) {
            if (!i %in% seq) {
                chr.lens[i] = 0
            }
        }
    }
    xlim = c(0, sum(chr.lens))
    chr.lens.sum = cumsum(chr.lens)
    chr.lens.sum = c(0, chr.lens.sum)
    plot.data$pos.updated = plot.data$pos + chr.lens.sum[plot.data$seq]
    chr.abline = chr.lens[seq]
    chr.abline.sum = cumsum(chr.abline)
    chr.abline.sum = c(0, chr.abline.sum)
    chr.axis = seq
    for (i in 1:length(seq)) {
        chr.axis[i] = (chr.abline.sum[i] + chr.abline.sum[i + 1])/2
    }
    par(mai = c(1, 1, 1, 1.5))
    plot.data$dis.log10 = log10(plot.data$dis)
    plot(plot.data$pos.updated, plot.data$dis.log10, col = col[plot.data$alteration], 
        pch = 16, cex = 0.5, xlim = xlim, axes = F, box(), mgp = c(2, 
            1, 0), main = paste("Rainfall plot for", sample, sep = " "), 
        xlab = "Chromosomes", ylab = "Intermutation distance (bp, log10)")
    abline(h = 2, col = "red", lty = 2)
    abline(v = chr.abline.sum, col = "grey")
    axis(2, las = 1, lwd.tick = 0.5)
    axis(1, at = chr.axis, labels = c(seq0), tick = FALSE, mgp = c(1, 
        0.5, 0))
    legend("right", names(col), col = col, pch = 16, inset = c(-0.15, 
        0), bty = "n", xpd = TRUE)
    arrows.point = c()
    
    num = dim(plot.data)[1] - 5
    katPoint <- matrix(nrow = num, ncol = 7)
    i = 1
    mutnum = 1
    Cmutnum = 0
    a1 = 0
    a2 = 0
    for (i in 1:num) {
	if (plot.data$ref[i] %in% c("C", "G")){
            Cmutnum = Cmutnum + 1
	}
        if (plot.data$dis[i + 1] <= max.dis) {
            mutnum = mutnum + 1
        } else {
            if (mutnum >= min.mut) {
                len = plot.data$pos[i] - plot.data$pos[i - mutnum + 
                  1] + 1
                chr.n = gsub(pattern = "chr", replacement = "", x = plot.data$chr[i], 
                  fixed = TRUE)
                chr.n = gsub(pattern = "X", replacement = "23", x = chr.n, 
                  fixed = TRUE)
                chr.n = gsub(pattern = "Y", replacement = "24", x = chr.n, 
                  fixed = TRUE)
                chr.n = as.numeric(chr.n)
                if (plot.data$pos[i] <= chr.arm[chr.n]) {
                  arm = paste(chr.n, "p", sep = "")
                } else if (plot.data$pos[i - mutnum + 1] >= chr.arm[chr.n]) {
                  arm = paste(chr.n, "q", sep = "")
                } else {
                  arm = paste(chr.n, "p, ", chr.n, "q", sep = "")
                }
		wei = round(Cmutnum/mutnum,3)
                katPoint[i, 1:7] = c(plot.data$chr[i], plot.data$pos[i - 
                  mutnum + 1], plot.data$pos[i], arm, len, mutnum, wei)
                if (wei>=0.8) {
                  a1 = a1 + 1
                  arrows.point[a1] = (plot.data$pos.updated[i] + plot.data$pos.updated[i - 
                    mutnum + 1])/2
		} else {
                  a2 = a2 + 1
                  arrows.point[a2] = (plot.data$pos.updated[i] + plot.data$pos.updated[i - 
                    mutnum + 1])/2
		}
            }
            mutnum = 1
	    Cmutnum = 0
        }
    }
    katPoint.out = data.frame(na.omit(katPoint))
    names(katPoint.out) = c("chr", "start", "end", "chr.arm", "length", "number.mut", 
        "weight.C>X")
    message(paste(dim(katPoint.out)[1], "potential kataegis events identified", 
        sep = " "))
    print(katPoint.out)
    if (a1 > 0) {
        for (k in 1:a1) {
            arrows(arrows.point[k], 0, arrows.point[k], 0.5, length = 0.1, 
                angle = 30, code = 2, col="black")
        }
    }
    if (a2 > 0) {
        for (k in 1:a2) {
            arrows(arrows.point[k], 0, arrows.point[k], 0.5, length = 0.1, 
                angle = 30, code = 2, col="grey")
        }
    }
}
