library(epiAneufinder)

setwd("/sharedFolder")

download.file(
  url = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz",
  destfile = "hg38-blacklist.v2.bed.gz",
  mode = "wb"
)

R.utils::gunzip(
  filename = "hg38-blacklist.v2.bed.gz",
  destname = "hg38-blacklist.v2.bed",
  overwrite = TRUE,
  remove = TRUE
)

file.exists("/sharedFolder/hg38-blacklist.v2.bed")

library(epiAneufinder)

fragment_file <- "/sharedFolder/Data/dv90_atac_fragments.tsv.gz"
blacklist_file <- "/sharedFolder/hg38-blacklist.v2.bed"
outdir <- "/sharedFolder/Results/epiAneufinder/DV90"

file.exists(fragment_file)
file.exists(blacklist_file)

con <- gzfile(fragment_file, "rt")
readLines(con, n = 5)
close(con)

epiAneufinder(
  input = fragment_file,
  outdir = outdir,
  blacklist = blacklist_file,
  windowSize = 1e6,
  genome = "BSgenome.Hsapiens.UCSC.hg38",
  exclude = c("chrX", "chrY", "chrM"),
  reuse.existing = FALSE,
  title_karyo = "DV90 epiAneufinder",
  ncores = 4,
  minFrags = 1000,
  minsizeCNV = 0,
  k = 4,
  plotKaryo = FALSE
)

print("ciao")


