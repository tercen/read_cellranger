# Read Cell Ranger

##### Description

The `Read Cell Ranger operator` reads Cell Ranger output files into Tercen.

##### Usage

Input projection|.
---|---
`column`        | factor, documentId of the ZIP file containing Cell Ranger output files

Input parameters|.
---|---
`proportion`        | numeric, proportion of data (barcodes) to downsample

Output relations|.
---|---
`count`         | count
`barcode`       | barcode ID
`gene_name1`    | gene ID (Ensembl)
`gene_name2`    | gene ID (Symbol)
`sample_id`     | sample ID (subfolder name)

