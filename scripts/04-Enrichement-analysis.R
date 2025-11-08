# 5- Gene Set Enrichemenet Analysis Using GO and KEGG----

gost.res_up <- gost(rownames(myModule_up), organism = "hsapiens", correction_method = "fdr")
Upregulated_erichemnet_plot <-gostplot(gost.res_up, interactive = FALSE, capped = T) 
ggsave("results/plots/`Upregulated_erichemnet.png", 
       plot = Upregulated_erichemnet_plot,
       width = 12, height = 8, dpi = 300, bg= "white")


gost.res_down <- gost(rownames(myModule_down), organism = "hsapiens", correction_method = "fdr")
Downregulated_erichemnet_plot <-gostplot(gost.res_down, interactive = FALSE, capped = T) 
ggsave("results/plots/`Downregulated_erichemnet.png", 
       plot = Downregulated_erichemnet_plot,
       width = 12, height = 8, dpi = 300, bg= "white")



hs_gsea_c2 <- msigdbr(species = "Homo sapiens", 
                      category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol) 

mydata.df <- log2.cpm.filtered.norm.df %>% 
  mutate(healthy.AVG = (HS01 + HS02 + HS03 + HS04 + HS05)/5,
         disease.AVG = (CL08 + CL10 + CL11 + CL12 + CL13)/5,
         LogFC = (disease.AVG - healthy.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)


# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)
write.csv(myGSEA.df, "results/tables/GSEA_results.csv", row.names = FALSE)

datatable(myGSEA.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)


GSE_result <- gseaplot2(myGSEA.res, 
          geneSetID = 47, 
          pvalue_table = FALSE, 
          title = myGSEA.res$Description[47]) 

ggsave("results/plots/`GSE_result.png", 
       plot = GSE_result,
       width = 12, height = 8, dpi = 300, bg= "white")

myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

Enreichment_bubble_plot <-ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()
ggsave("results/plots/`Enreichment_bubble_plot.png", 
       plot = Enreichment_bubble_plot,
       width = 12, height = 8, dpi = 300, bg= "white")


