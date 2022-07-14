suppressPackageStartupMessages({
  library(tercen)
  library(tercenApi)
  library(dplyr)
  library(Seurat)
})
ctx = tercenCtx()

if (!any(ctx$cnames == "documentId")) stop("Column factor documentId is required.")

# extract files
df <- ctx$cselect()

docId = df$documentId[1]
doc = ctx$client$fileService$get(docId)
filename = tempfile()
writeBin(ctx$client$fileService$download(docId), filename)
on.exit(unlink(filename))

tmpdir <- tempfile()
unzip(filename, exdir = tmpdir)
f.names <- list.files(tmpdir, full.names = TRUE, recursive = TRUE)
f.names <- f.names[grep("barcodes.tsv|features.tsv|matrix.mtx", f.names)]

samp_names <- gsub("/barcodes.tsv*|/features.tsv*|/matrix.mtx*|.gz", "", f.names) %>%
  unique

nms <- do.call(rbind, strsplit(samp_names, "/"))
idx <- apply(nms, 2, function(x) length(unique(x)) > 1)
if(sum(idx) == 1) {
  sample_names <- nms[, idx]
} else {
  sample_names <- apply(nms[, idx], 1, paste0, collapse="-")
} 

names(samp_names) <- sample_names
names(sample_names) <- samp_names

seurat_list <- lapply(samp_names, function(x){
  raw <- Seurat::Read10X(x)
  if(length(raw) > 1) {  # For output from CellRanger >= 3.0
    raw = CreateSeuratObject(counts = raw$`Gene Expression`, project = sample_names[x])
    # raw[['Protein']] = CreateAssayObject(counts = raw$`Antibody Capture`)
  } else {
    raw = CreateSeuratObject(counts = raw, project = sample_names[x])
  }
  raw
})

names(seurat_list) <- sample_names
mat_tmp <- do.call(rbind, lapply(seurat_list, dim)) 
colnames(mat_tmp) <- c("n_genes", "n_cells")
df_tmp <- mat_tmp %>%
  as_tibble() %>%
  mutate(.ci = 0, sample_name = sample_names)

merged_seurat <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = sample_names,
  project = "Project"
)

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(
  merged_seurat,
  pattern = "^MT-"
)

## sparse matrix to data frame
df_out <- as.data.frame(summary(GetAssayData(merged_seurat))) %>% 
  as_tibble() %>%
  rename(gene_id = i, cell_id = j, value = x) %>%
  mutate(.ci = 0L) %>%
  ctx$addNamespace()

ctx$save(df_out)
