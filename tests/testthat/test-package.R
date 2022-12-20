
# ------------------- Prepare original KataegisPortal results -----------------

# # Load packages
# library(KataegisPortal)
# library(BSgenome)
# library(BSgenome.Hsapiens.UCSC.hg19)
#
# # Read in snv data
# mutData <- system.file("examples", "mutData.txt", package="KataegisPortal")
# mutData <- read.table(mutData, header = TRUE,sep = "\t",as.is = TRUE)
# head(mutData)
#
# # Covert to KataegisPortal input
# mutSNP = mutSNP.input(mut.data = mutData,
#                       chr = "chr",
#                       pos = "pos",
#                       ref = "ref",
#                       alt = "alt",
#                       build = "hg19")
#
# # Identify potential kataegis events
# kataegis_expected <- katPoint(mutSNP)
# kataegis_expected %>%
#   readr::write_tsv("tests/testdata/kataegis_expected.tsv")


# ---------------------- Test our changes -------------------------------------

test_that("Package results are identical to the original ones", {
  kataegis_expected <- "../testdata/kataegis_expected.tsv" |>
    readr::read_tsv(col_types = "ccccccccd") |>
    as.data.frame()

  mutData <- system.file("examples", "mutData.txt", package="KataegisPortal")
  mutData <- read.table(mutData, header = TRUE,sep = "\t",as.is = TRUE)

  # Covert to KataegisPortal input
  mutSNP <- mutSNP.input(mut.data = mutData,
                        chr = "chr",
                        pos = "pos",
                        ref = "ref",
                        alt = "alt",
                        build = "hg19")

  # Identify potential kataegis events
  kataegis_res <- katPoint(mutSNP)
  expect_identical(kataegis_res, kataegis_expected)
})
