###
BiocManager::install("simplifyEnrichment")
library(simplifyEnrichment)


###from DEG-GO_analysis.Rmd


DAsptable3 <- DAspgostres3$result %>%
  dplyr::select(c(source, term_id,
                  term_name,intersection_size,
                  term_size, p_value))# %>%
DAspGO <- DAsptable3$term_id
mat =GO_similarity(DAspGO)
simplifyGO(mat)
DAsp3genelist <- DAspgostres3$meta$genes_metadata$query$query_1
write_tsv(as.data.frame(DAsp3genelist),"output/DAsp3genelist.txt")
