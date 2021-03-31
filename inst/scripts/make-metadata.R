### =========================================================================
### GSE13015_GPL6106 metadata
### -------------------------------------------------------------------------
###

meta <- data.frame(
  Title = "GSE13015_GPL6106_Microaray_SEPSIS",
  Description = paste0("Microarray expression data",
                       "for 67 septicemic patients, represented as an ",
                       "SummarizedExperiment R data representation derived from ",
                       "GEO accession GSE13015, GPL6106."),
  BiocVersion = "3.13",
  Genome = "hg19",
  SourceType = "tar.gz",
  SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13015",
  SourceVersion = "Jul 1 2009",
  Species = "Homo sapiens",
  TaxonomyId = NA,
  Coordinate_1_based = NA,
  DataProvider = "GEO",
  Maintainer = "Darawan Rinchai <drinchai@gmail.com>",
  RDataClass = "SummarizedExperiment",
  DispatchClass = "Rda",
  RDataPath = "GSE13015/GSE13015_GPL6106_SummarizedExperiment.Rda" ,
  Tags = "",
  Notes = "Example data for BloodGen3Module Package"
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
