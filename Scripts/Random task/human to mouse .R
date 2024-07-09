library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = 'uswest.ensembl.org', ensemblRedirect = FALSE)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = 'uswest.ensembl.org', ensemblRedirect = FALSE)
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = a$VERHAAK_GLIOBLASTOMA_CLASSICAL, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

mart <- useMart("ensembl","mmusculus_gene_ensembl", host = 'uswest.ensembl.org', ensemblRedirect = FALSE)
xxx <- getBM( attributes= c("mgi_symbol",'ensembl_gene_id'),values=genesV2$MGI.symbol,mart= mart,filters= 'mgi_symbol')