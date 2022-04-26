library(tercen)
library(dplyr)

ctx = tercenCtx()

if (!any(ctx$cnames == "documentId")) stop("Column factor documentId is required.")

# extract files
df <- ctx$cselect()

docId = df$documentId[1]
doc = ctx$client$fileService$get(docId)
filename = tempfile()
writeBin(ctx$client$fileService$download(docId), filename)
on.exit(unlink(filename))

# unzip if archive
if(length(grep(".zip", doc$name)) > 0) {
  tmpdir <- tempfile()
  unzip(filename, exdir = tmpdir)
  f.names <- list.files(tmpdir, full.names = TRUE, recursive = TRUE)
  f.names <- f.names[grep("barcodes.tsv|features.tsv|matrix.mtx", f.names)]
} else {
  f.names <- filename
}

folders <- unique(dirname(f.names))

proportion = ctx$op.value('proportion', as.double, 1)

# check matrix, barcode and gene file
d_out <- lapply(folders, function(i) {
  f.names_sample <- f.names[grep(i, f.names)]
  
  # convert them to one matrix file
  matrix_table <- read.delim(
    file = f.names_sample[grepl("matrix.mtx.gz", f.names_sample)],
    sep = " ",
    header = FALSE,
    skip = 2
  ) %>%
    rename("feature_idx" = V1, "barcode_idx" = V2, "count" = V3)
  
  barcode_table <- read.delim(
    file = f.names_sample[grepl("barcodes", f.names_sample)],
    sep = "\t",
    header = FALSE
  ) %>% 
    mutate(row_idx = 1:nrow(.)) %>% 
    rename(barcode = V1) %>%
    slice_sample(prop = proportion)
  
  if(nrow(barcode_table) == 0) stop("No cells in the dataset. Try to increase the proportion of sampled cells.")
  
  feature_table <- read.delim(
    file = f.names_sample[grepl("genes|features", f.names_sample)],
    sep = "\t",
    header = FALSE
  ) %>% 
    mutate(row_idx = 1:nrow(.)) %>% 
    rename(feature = V1)
  
  matrix_table <- matrix_table %>%
    inner_join(barcode_table, by= c("barcode_idx" = "row_idx")) %>%
    left_join(feature_table, by= c("feature_idx" = "row_idx")) %>%
    rename(gene_name1 = feature, gene_name2 = V2)
  
  matrix_table <- matrix_table %>%
    select(-ends_with("_idx")) %>%
    select(count, barcode, gene_name1, gene_name2) %>%
    mutate(sample_id = gsub(tmpdir, "", i))
  
  return(matrix_table)
})

df_out <- do.call(rbind, d_out)

df_out %>%
  mutate_if(is.integer, as.double) %>%
  mutate(filename = basename(tmpdir)) %>%
  mutate(.ci = 0) %>%
  ctx$addNamespace() %>%
  ctx$save()
