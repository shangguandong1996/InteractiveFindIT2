# InteractiveFindIT2
InteractiveFindIT2 is an R package that help user of FindIT2 explore FindIT2 result 
interactively.

# Install

```
library(devtools)
install_github("shangguandong1996/InteractiveFindIT2")
```

# Usage

```
library(TxDb.Athaliana.BioMart.plantsmart28)
library(FindIT2)
library(InteractiveFindIT2)

Txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))

data("RNA_normCount")
data("ATAC_normCount")
data("test_geneSet")
peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
peak_GR <- loadPeakFile(peak_path)

ChIP_peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
ChIP_peak_GR <- loadPeakFile(ChIP_peak_path)
ChIP_peak_GR$TF_id <- "AT1G28300"

```

```
mmAnno <- mm_geneScan(peak_GR, Txdb)
```

**shinyParse_peakGeneCor**

```
ATAC_colData <- data.frame(
    row.names = colnames(ATAC_normCount),
    type = gsub("_R[0-9]", "", colnames(ATAC_normCount))
)
ATAC_normCount_merge <- integrate_replicates(ATAC_normCount, ATAC_colData)
RNA_colData <- data.frame(
    row.names = colnames(RNA_normCount),
    type = gsub("_R[0-9]", "", colnames(RNA_normCount))
)
RNA_normCount_merge <- integrate_replicates(RNA_normCount, RNA_colData)
mmAnnoCor <- peakGeneCor(
    mmAnno = mmAnno,
    peakScoreMt = ATAC_normCount_merge,
    geneScoreMt = RNA_normCount_merge,
    parallel = FALSE
)
shinyParse_peakGeneCor(mmAnnoCor)
```


**shinyParse_findIT_regionRP**

```
regionRP <- calcRP_region(
    mmAnno = mmAnno,
    peakScoreMt = ATAC_normCount,
    Txdb = Txdb,
    Chrs_included = "Chr5"
)
set.seed(20160806)
result_findIT_regionRP <- findIT_regionRP(
    regionRP = regionRP,
    Txdb = Txdb,
    TF_GR_database = ChIP_peak_GR,
    input_genes = test_geneSet,
    background_number = 3000
)
merge_result <- c(regionRP, result_findIT_regionRP)
shinyParse_findIT_regionRP(merge_result)
```
