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
f.names <- f.names[grep("barcodes.tsv|features.tsv|genes.tsv|matrix.mtx", f.names)]

samp_names <- gsub("/barcodes.tsv*|/features.tsv*|/genes.tsv*|/matrix.mtx*|.gz", "", f.names) %>%
  unique

nms <- do.call(rbind, strsplit(samp_names, "/"))
idx <- apply(nms, 2, function(x) length(unique(x)) > 1)
if(sum(idx) == 0) {
  sample_names <- tail(c(nms), 1)
} else if(sum(idx) == 1) {
  sample_names <- nms[, idx]
} else {
  sample_names <- apply(nms[, idx], 1, paste0, collapse="-")
} 

names(samp_names) <- sample_names
names(sample_names) <- samp_names

seurat_list <- lapply(samp_names, function(x){
  raw <- Seurat::Read10X(data.dir = x)
  if(!class(raw) == "dgCMatrix") {  # For output from CellRanger >= 3.0
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

if(length(seurat_list) > 1) {
  merged_seurat <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = sample_names,
    project = "Project"
  )
} else {
  merged_seurat <- seurat_list[[1]]
}

## sparse matrix to data frame
spm <- GetAssayData(merged_seurat)
df_out <- as.data.frame(summary(spm)) %>% 
  as_tibble() %>%
  rename(gene_id = i, cell_id = j, value = x) %>%
  mutate(.gene_id = gene_id, .cell_id = cell_id) %>%
  mutate(.ci = 0L) %>%
  ctx$addNamespace()

gene_names <- dimnames(spm)[[1]]
df_gene <- tibble(.gene_id = seq_along(gene_names), gene_names = gene_names)

cell_names <- dimnames(spm)[[2]]
df_cell <- tibble(.cell_id = seq_along(cell_names), cell_names = cell_names)

data_relation <- df_out %>% as_relation()
gene_relation <- df_gene %>% ctx$addNamespace() %>% as_relation()
cell_relation <- df_cell %>% ctx$addNamespace() %>% as_relation()

rel_out <- data_relation %>%
  left_join_relation(gene_relation, ".gene_id", ".gene_id") %>%
  left_join_relation(cell_relation, ".cell_id", ".cell_id") %>%
  as_join_operator(list(), list())

save_relation(rel_out, ctx)
