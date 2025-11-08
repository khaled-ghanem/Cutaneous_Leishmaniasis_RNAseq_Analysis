# 0- Import el packages----
library(tidyverse)
library(tximport) 
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(matrixStats)
library(cowplot)
library(DT)#interactive and searchable tables of our GSEA results
library(gt)
library(plotly)
library(ggrepel)
library(limma)
library(gplots)
library(RColorBrewer)
library(IsoformSwitchAnalyzeR)
library(GSEABase) 
library(Biobase) 
library(GSVA) 
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots

# 1- Read el Data-----
targets <- read_tsv("data/studydesign.txt")
path <- file.path("results/kallisto/", targets$sample,"abundance.tsv" )
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))# Take tx2gene from it 
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, 
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)




sampleLabels <- targets$sample

# 2- Apply Filtering and Normalization----
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, 
                                  cols = HS01:CL13, 
                                  names_to = "samples", 
                                  values_to = "expression") 


# Plot to clearfy the effect of TMM normalization and filtering low cpm data 
p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=5 # define this to match your needs
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, 
                                           cols = HS01:CL13, 
                                           names_to = "samples", 
                                           values_to = "expression") 

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = HS01:CL13, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

combined_plot <- plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

ggsave("results/plots/CPM_violin_TMM_filtering_comparison.png", 
       plot = combined_plot,
       width = 12, height = 8, dpi = 300, bg= "white")



# Plot PCA result to clearfy two groups
group <- targets$group
group <- factor(group)

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale. = FALSE, retx = TRUE)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var / sum(pc.var) * 100, 1)

pca.res.df <- as_tibble(pca.res$x)
pca.res.df$sra_accession <- targets$sra_accession
pca.res.df$group <- factor(targets$group)

pca.plot <- ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, color = group, label = sra_accession) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  stat_ellipse() +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(title = "PCA Plot: TMM-Normalized Expression",
       caption = paste0("Produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggsave("results/plots/PCA_TMM_disease_VS_healthy.png", 
       plot = pca.plot,
       width = 12, height = 8, dpi = 300, bg= "white")


