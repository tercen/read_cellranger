library(tercen)
library(dplyr)


ctx = tercenCtx()

if (!any(ctx$cnames == "documentId")) stop("Column factor documentId is required")

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
} else {
  f.names <- filename
}

# check matrix, barcode and gene file

# convert them to one matrix file
percentage = as.double(ctx$op.value('percentage'))

print(percentage)

stop("check the percentage")

matrix_table <- read.delim(file = f.names[grepl("matrix.mtx",f.names)], sep = " ", header = FALSE, skip=2) %>%
  rename("feature_idx"= V1, "barcode_idx" =V2, "count" =V3)

barcode_table <- read.delim(file = f.names[grepl("barcodes.tsv",f.names)], sep = "\t", header = FALSE) %>% 
  mutate(row_idx = 1:nrow(.)) %>% 
  rename(barcode = V1) %>% slice_sample(prop = percentage)

feature_table <- read.delim(file = f.names[grepl("genes.tsv",f.names)], sep = "\t", header = FALSE) %>% 
  mutate(row_idx = 1:nrow(.)) %>% 
  rename(feature = V1)

matrix_table <- matrix_table %>%
  right_join(barcode_table, by= c("barcode_idx" = "row_idx"))

matrix_table <- matrix_table %>%
  left_join(feature_table, by= c("feature_idx" = "row_idx")) %>%
  rename(gene_name1 = feature, gene_name2 = V2)

matrix_table <- matrix_table %>%
  select(-ends_with("_idx"))

matrix_table %>%
  mutate_if(is.integer, as.double) %>%
  mutate(filename = basename(tmpdir)) %>%
  mutate(.ci = 0) %>%
  ctx$addNamespace() %>%
  ctx$save()
