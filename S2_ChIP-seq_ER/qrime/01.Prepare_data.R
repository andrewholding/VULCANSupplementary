suppressMessages(library(dplyr, quietly = TRUE))

metadata <- read.delim("SampleData/metadata.txt", check.names = FALSE)
metadata <- metadata %>%
  mutate(Sample = paste(Run, Group, sep = ":"))

peptideData <- read.delim("SampleData/200116_peptide_intensities_new_reprocessing.txt", stringsAsFactors = FALSE, check.names = FALSE)
peptideData <- peptideData %>%
  rename(Sequence = `Annotated Sequence`, Protein = `Master Protein Accessions`)

retainedColumns <- c("Sequence", "Modifications", "Protein")

runs <- levels(metadata$Run)

for (run in runs)
{
  cat(run, "\n", sep = "")

  runMetadata <- metadata %>% filter(Run == run)

  samples <- runMetadata$Sample
  groups <- runMetadata$Group

  runData <- peptideData %>%
    select(one_of(retainedColumns), one_of(samples)) %>%
    mutate(Missing = apply(select(., one_of(samples)), 1, function(x) sum(is.na(x)))) %>%
    filter(Missing != length(samples)) %>%
    select(-Missing)

  colnames(runData) <- c(retainedColumns, as.character(groups))

  write.table(runData, file = paste("ProcessedData/", run, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
}

