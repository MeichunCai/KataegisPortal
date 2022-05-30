mutSNP.input <- function(mut.data, chr = "chr", pos = "pos", ref = "ref", 
    alt = "alt", build = NULL, k = 10) {
    if (exists("mut.data", mode = "list")) {
        mut.data <- mut.data
    } else {
        if (file.exists(mut.data)) {
            mut.data <- utils::read.table(mut.data, sep = "\t", header = TRUE, 
                as.is = FALSE, check.names = FALSE)
        } else {
            stop("mut.data is neither a file nor a loaded data frame")
        }
    }
    mut.data <- mut.data[, c(chr, pos, ref, alt)]
    genome.opts = c("hg19", "hg18", "hg38","CHM13")
    if (!build %in% genome.opts || is.null(build)) {
        stop("Available reference builds: hg18, hg19, hg38,CHM13")
    }
    if (build == "hg19") {
        chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 
            171115067, 159138663, 146364022, 141213431, 135534747, 
            135006516, 133851895, 115169878, 107349540, 102531392, 
            90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 
            51304566, 155270560, 59373566)
        bsg = BSgenome.Hsapiens.UCSC.hg19
    } else if (build == "hg18") {
        chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 
            170899992, 158821424, 146274826, 140273252, 135374737, 
            134452384, 132349534, 114142980, 106368585, 100338915, 
            88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 
            49691432, 154913754, 57772954)
        bsg = BSgenome.Hsapiens.UCSC.hg18
    } else if (build == "hg38") {
        chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 
            170805979, 159345973, 145138636, 138394717, 133797422, 
            135086622, 133275309, 114364328, 107043718, 101991189, 
            90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 
            50818468, 156040895, 57227415)
        bsg = BSgenome.Hsapiens.UCSC.hg38
    } else if (build == "CHM13") {
            chr.lens = c(248387328,242696752,201105948,193574945,182045439,
            172126628,160567428,146259331,150617247,134758134,135127769,
            133324548,113566686,101161492,99753195,96330374,84276897,
            80542538,61707364,
            66210255,45090682,51324926,154259566,62460029)
        bsg = BSgenome.Hsapiens.CHM13
    } else {
        stop("Available reference builds: hg18, hg19, hg38,CHM13")
    }
    mut.data$build = build
    if (!all(mut.data$ref %in% DNA_BASES & mut.data$alt %in% DNA_BASES)) {
        stop("Only SNV substitutions are currently supported.")
    }
    ref_base = DNAStringSet(mut.data$ref)
    alt_base = DNAStringSet(mut.data$alt)
    conv.start = mut.data$pos - k
    conv.end = mut.data$pos + k
    context = getSeq(bsg, mut.data$chr, start = conv.start, end = conv.end)
    if (TRUE) {
        idx = mut.data$ref %in% c("A", "G")
        context[idx] = reverseComplement(context[idx])
        ref_base[idx] = reverseComplement(ref_base[idx])
        alt_base[idx] = reverseComplement(alt_base[idx])
    }
    mut.data$alteration = paste(ref_base, alt_base, sep = ">")
    mut.data$context = context
    # Replace chr x and y with numeric value (23 and 24) for better
    # ordering
    seq = gsub(pattern = "chr", replacement = "", x = mut.data$chr, 
        fixed = TRUE)
    seq = gsub(pattern = "X", replacement = "23", x = seq, fixed = TRUE)
    seq = gsub(pattern = "Y", replacement = "24", x = seq, fixed = TRUE)
    mut.data$seq = as.numeric(seq)
    mut.data = mut.data[order(mut.data$seq, mut.data$pos), ]
    chr.lens.sum = cumsum(chr.lens)
    chr.lens.sum = c(0, chr.lens.sum)
    mut.data$dis = c(mut.data$pos[1], diff(mut.data$pos + chr.lens.sum[mut.data$seq]))
    return(mut.data)
}
