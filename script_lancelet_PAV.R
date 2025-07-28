library(tidyverse)
library(ggplot2)
library(ggpubr)
library(smplot2)
library(data.table)
library(ggrepel)
library(VennDiagram)
library(UpSetR)
library(FactoMineR)
library(ggcorrplot)
library(corrr)
library(factoextra)
library(seriation)
library(heatmaply)
library(ggridges)
#packages = c('seriation', 'heatmaply', 'tidyverse')

setwd("~/R_working/LANCELET") 
data <- read.csv2("meta.csv") # LOAD METADATA

################ LOAD, COUNT & SELECT PAV GENES AND REGIONS (B. belcheri) ########
Bb_PAV <- read.csv2("Bb_PAV.txt") # LOAD THE LIST OF MISSING GENES IN ALL THE SAMPLES
Bb_PAV_5kb <- read.csv2("Bb_PAV_5kb.txt") # LOAD THE LIST OF MISSING REGIONS IN ALL THE SAMPLES
Bb_PAV_count <- Bb_PAV %>% count(Gene) # COUNT THE OCCURRENCES PER GENE
Bb_PAV_count_count <- Bb_PAV_count %>% count(n) # COUNT THE DISTRIBUTION OF OCCURRENCES
Bb_PAV_sel<-Bb_PAV_count[Bb_PAV_count$n>10,] # FILTER DISPENSABLE GENES (<90%)
Bb_PAV_5kb_count <- Bb_PAV_5kb %>% count(Gene) # COUNT THE OCCURRENCES PER REGION
Bb_PAV_5kb_count_count <- Bb_PAV_5kb_count %>% count(n)
Bb_PAV_5kb_sel<-Bb_PAV_5kb_count[Bb_PAV_5kb_count$n>10,]

################ LOAD, COUNT & SELECT PAV GENES AND REGIONS (B. floridae) ########
Bf_PAV <- read.csv2("Bf_PAV.txt")
Bf_PAV_5kb <- read.csv2("Bf_PAV_5kb.txt")
Bf_PAV_d <- read.csv2("Bf_PAV_doughter.txt")
Bf_genes <- read.csv2("Bf_genes.txt")
Bf_PAV_count <- Bf_PAV %>% count(Gene)
Bf_PAV_count_count <- Bf_PAV_count %>% count(n)
Bf_PAV_sel<-Bf_PAV_count[Bf_PAV_count$n>4,]
Bf_PAV_count_d <- Bf_PAV_d %>% count(Gene)
Bf_PAV_sel_d<-Bf_PAV_count_d[Bf_PAV_count_d$n>9,]
Bf_PAV_sel_d50<-Bf_PAV_count_d[Bf_PAV_count_d$n>9,]
Bf_PAV_5kb_count <- Bf_PAV_5kb %>% count(Gene)
Bf_PAV_5kb_sel<-Bf_PAV_5kb_count[Bf_PAV_5kb_count$n>4,]

empVec <- character()
vec1 = c(empVec, Bf_PAV_sel_d50$Gene, na.rm=TRUE)

#empVec <- character()
#vec2 = c(empVec, Bf_PAV_sel$Gene, na.rm=TRUE)

############################### HEATMAP & PCA (B. belcheri) #########################

Bb_matrix_REGION <- read.csv2("Bb_matrix_5kb.csv")  # LOAD THE MATRIX OF MISSING REGIONS IN ALL THE SAMPLES 
Bb_matrix_PAV_R <- inner_join(Bb_matrix_REGION, Bb_PAV_5kb_sel, by = c("Gene"))
row.names(Bb_matrix_PAV_R) <- Bb_matrix_PAV_R$Gene
Bb_matrix_PAV_R_s <- select(Bb_matrix_PAV_R, c(2, 2:101))
Bb_matrix_REGION <- data.matrix(Bb_matrix_PAV_R_s)

#heatmaply(Bb_matrix_REGION, colors = c("orange", "white"), k_row = 8, k_col = 3, column_text_angle = 90, cexRow = 0.8, margins = c(50,50,50,50), 
#          fontsize_row =6, fontsize_col = 10, main="", xlab = "Samples", ylab = "PAV regions (5kb)", width = 1200,
#          height = 750, seriate = "OLO", plot_method = "plotly", file = "HEATMAPS/bb_matrix_regions.html" )

Bb_matrix_gene <- read.csv2("Bb_matrix.csv")
Bb_matrix_PAV_g <- inner_join(Bb_matrix_gene, Bb_PAV_sel, by = c("Gene"))
row.names(Bb_matrix_PAV_g) <- Bb_matrix_PAV_g$Gene
Bb_matrix_PAV_g_s <- select(Bb_matrix_PAV_g, c(2, 2:101))
Bb_matrix_gene <- data.matrix(Bb_matrix_PAV_g_s)

#heatmaply(Bb_matrix_gene, colors = c("darkgreen", "white"), k_row = 8, k_col = 3, column_text_angle = 90, cexRow = 0.8, margins = c(50,50,50,50), 
#        fontsize_row =6, fontsize_col = 10, main="", xlab = "Samples", ylab = "Dispensable Genes", width = 1500,
#         height = 1000, seriate = "OLO", plot_method = "plotly", file = "HEATMAPS/bb_matrix_gene.html")


colSums(is.na(Bb_matrix_gene))
corr_matrix <- cor(Bb_matrix_gene)
CORR_Bb <- ggcorrplot(corr_matrix,tl.cex = 4,
                      tl.col = "black",
                      tl.srt = 45,)
CORR_Bb
data.pca <- princomp(corr_matrix)
summary(data.pca)
fviz_eig(data.pca, addlabels = TRUE)
# Graph of the variables
fviz_pca_var(data.pca, col.var = "black")
fviz_cos2(data.pca, choice = "var", axes = 1:2)
fviz_pca_var(data.pca, col.var = "cos2", geom = c("point"), #, "text"
             gradient.cols = c("grey", "orange", "blue"),
             repel = TRUE)

fviz_pca_biplot(data.pca, label="var", 
                addEllipses=FALSE, ellipse.level=0.7)
fviz_pca_ind(data.pca, geom = c("point"),repel = TRUE)
fviz_pca_var(data.pca, geom = c("point", "text"),repel = TRUE)


metadata_bb <- data[data$Species == "Bb",]
aa <- metadata_bb %>% column_to_rownames("Run")
aa<-aa[names(data.pca$scale),]

PCA_Bb <- fviz_pca_ind(data.pca,  habillage=aa$MALACO, 
             addEllipses=TRUE, ellipse.level=0.7, geom = c("point"),repel = TRUE)
PCA_Bb
names(data.pca$scale)

############################### HEATMAP & PCA (B. floridae) #########################

Bf_matrix_gene <- read.csv2("Bf_matrix.csv")
Bf_matrix_PAV_g <- inner_join(Bf_matrix_gene, Bf_PAV_sel, by = c("Gene"))
row.names(Bf_matrix_PAV_g) <- Bf_matrix_PAV_g$Gene
Bf_matrix_PAV_g_s <- select(Bf_matrix_PAV_g, c(2, 2:141))
Bf_matrix_gene <- data.matrix(Bf_matrix_PAV_g_s)

#heatmaply(Bf_matrix_gene, colors = c("darkgreen", "white"), k_row = 8, k_col = 3, column_text_angle = 90, cexRow = 0.8, margins = c(50,50,50,50), 
#        fontsize_row =6, fontsize_col = 10, main="", xlab = "Samples", ylab = "Dispensable Genes", width = 1500,
#         height = 1000, seriate = "OLO", plot_method = "plotly", file = "HEATMAPS/bb_matrix_gene.html")


colSums(is.na(Bf_matrix_gene))
corr_matrix <- cor(Bf_matrix_gene)
CORR_Bf <- ggcorrplot(corr_matrix,tl.cex = 4,
                      tl.col = "black",
                      tl.srt = 45,)
CORR_Bf
data.pca <- princomp(corr_matrix)
summary(data.pca)
fviz_eig(data.pca, addlabels = TRUE)
# Graph of the variables
fviz_pca_var(data.pca, col.var = "black")
fviz_cos2(data.pca, choice = "var", axes = 1:2)
fviz_pca_var(data.pca, col.var = "cos2", geom = c("point"), #, "text"
             gradient.cols = c("grey", "orange", "blue"),
             repel = TRUE)

fviz_pca_biplot(data.pca, label="var", 
                addEllipses=FALSE, ellipse.level=0.7)
fviz_pca_ind(data.pca, geom = c("point"),repel = TRUE)
fviz_pca_var(data.pca, geom = c("point", "text"),repel = TRUE)


metadata_bf <- read.csv2("metadata_bf.csv")
ab <- metadata_bf %>% column_to_rownames("Run")
ab<-ab[names(data.pca$scale),]

fviz_pca_biplot(data.pca, label="var", 
                addEllipses=FALSE, ellipse.level=0.7)
fviz_pca_ind(data.pca)
fviz_pca_var(data.pca, geom = c("point", "text"),repel = TRUE)

PCA_Bf <-fviz_pca_ind(data.pca,  habillage=ab$INHER, 
             addEllipses=TRUE, ellipse.level=0.7, geom = c("point"), repel = TRUE)

FIG_PCA <- ggarrange(PCA_Bb,CORR_Bb,PCA_Bf,CORR_Bf, ncol = 2, nrow = 2, labels = c("auto"), widths = c(1,1), 
                     heights = c(1,1), common.legend = TRUE, legend = c("none"), font.label = list(size = 20, color = "black"))
FIG_PCA

Bf_matrix_doughter <- read.csv2("Bf_matrix_genes_doughter.csv") ## LOAD THE MATRIX REALTED TO PARENT-OFFSPRING PAV
Bf_10_d <- inner_join(Bf_matrix_doughter, Bf_PAV_sel_d, by = c("Gene"))
row.names(Bf_10_d) <- Bf_10_d$Gene
Bf_s_d <- select(Bf_10_d, c(2, 2:97))
Bf_matrix_d <- data.matrix(Bf_s_d)

#heatmaply(Bf_matrix_d, colors = Blues, k_row = 12, k_col = 2, cexRow = 0.8, margins = c(100,100,100,100), 
#          fontsize_row =1, fontsize_col = 8, main="lupin PAV", xlab = "Samples", ylab = "Dispensable Genes", 
#          seriate = "OLO", plot_method = "plotly", file = "bf_gene_matrix_d.html")

############################ ENRICHMENT ANALYSIS: PFAM & GO ########################

Bb_vol <- read.csv2("Bb_vol_MOD.csv")
vec3 <- Bb_vol[Bb_vol$diffexpressed== "UP",]
vec3 <- as.vector(vec3$ID)

Bb_vol_f <- ggplot(Bb_vol, aes(x=log2(value), y=-log10(FDR), col=diffexpressed, label=delabel)) + geom_point() + theme_minimal()+ geom_text_repel() +
  theme(legend.position="none")+ xlim(-3,5)+
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 14))+
  geom_vline(xintercept=c(-2, 2), col="black", linetype="dashed", lwd=1) +
  geom_hline(yintercept=-log10(0.01), col="black", linetype="dashed", lwd=1)+
  scale_color_manual(values=c("grey", "#FF6666", "green")) + xlab("Enrichment FC (log2)")
Bb_vol_f


Bf_vol <- read.csv2("bf_vol_MOD.csv")
vec4 <- Bf_vol[Bf_vol$diffexpressed== "UP",]
vec4 <- as.vector(vec4$ID)

Bf_vol_f <- ggplot(Bf_vol, aes(x=log2(value), y=-log10(FDR), col=diffexpressed, label=delabel)) + geom_point() + theme_minimal()+ geom_text_repel() +
  theme(legend.position="none")+ xlim(-3,5)+
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 14))+
  geom_vline(xintercept=c(-2, 2), col="black", linetype="dashed", lwd=1) +
  geom_hline(yintercept=-log10(0.01), col="black", linetype="dashed", lwd=1)+
  scale_color_manual(values=c("grey", "#33CCCC", "green")) + xlab("Enrichment FC (log2)")
Bf_vol_f

Bb_GO <- read.csv2("bb_GO_MOD.csv")
vec5 <- Bb_GO[Bb_GO$diffexpressed== "UP" & Bb_GO$class== "biological_process",]
vec5 <- as.vector(vec5$GOID)
vec7 <- Bb_GO[Bb_GO$diffexpressed== "UP" & Bb_GO$class== "molecular_function",]
vec7 <- as.vector(vec7$GOID)

Bb_GO_f <- ggplot(Bb_GO, aes(x=log2(value), y=-log10(FDR), , col=class, label=delabel)) + geom_point() + theme_minimal()+ geom_text_repel() +
  theme(legend.position="none")+ xlim(0,7)+
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 14))+
  geom_vline(xintercept=c(-2, 1.41), col="black", linetype="dashed", lwd=1) +
  geom_hline(yintercept=-log10(0.01), col="black", linetype="dashed", lwd=1)+
  scale_color_manual(values=c("orange", "#FF6666", "blue")) + xlab("Enrichment FC (log2)")
Bb_GO_f

Bf_GO <- read.csv2("bf_GO_MOD.csv")
vec6 <- Bf_GO[Bf_GO$diffexpressed== "UP" & Bf_GO$class == "biological_process",]
vec6 <- as.vector(vec6$GOID)
vec8 <- Bf_GO[Bf_GO$diffexpressed== "UP" & Bf_GO$class == "molecular_function",]
vec8 <- as.vector(vec8$GOID)

Bf_GO_f <- ggplot(Bf_GO, aes(x=log2(value), y=-log10(FDR), col=class, label=delabel)) + geom_point() + theme_minimal()+ geom_text_repel() +
  theme(legend.position="none")+ xlim(0,7)+
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 14))+
  geom_vline(xintercept=c(-2, 1.41), col="black", linetype="dashed", lwd=1) +
  geom_hline(yintercept=-log10(0.01), col="black", linetype="dashed", lwd=1)+
  scale_color_manual(values=c("orange", "#FF6666", "blue")) + xlab("Enrichment FC (log2)")
Bf_GO_f

venn.diagram(x = list(vec3, vec4), category.names = c("B. belcheri" , "B. floridae"),
  filename = 'common PFAM.png', main = "Pfam", main.fontfamily = "sans", main.col = "black",
  main.cex = 2.5,
  output=TRUE, imagetype="png" ,  height = 1000,   width = 1000,   resolution = 300,  compression = "lzw",
  lwd = 2,  col=c("#FF6666", '#33CCCC'),
  fill = c(alpha("#FF6666",0.3), alpha('#33CCCC',0.3)),
  cex = 1.5,  fontfamily = "sans", cat.fontface = "italic", cat.cex = 1.2,  cat.default.pos = "outer",
  cat.pos = c(-27, 27),  cat.dist = c(0.055, 0.055),  cat.fontfamily = "sans",  cat.col = c("black", 'black'))

venn.diagram(x = list(vec5, vec6), category.names = c("B. belcheri" , "B. floridae"),
             filename = 'common GO_BP.png', main = "Biological process",
             main.fontfamily = "sans", main.col = "orange", main.cex = 1.5,
             output=TRUE, imagetype="png" ,  height = 650 ,   width = 650 ,   resolution = 300,  compression = "lzw",
             lwd = 2,  col=c("#FF6666", '#33CCCC'),
             fill = c(alpha("#FF6666",0.3), alpha('#33CCCC',0.3)),
             cex = 1.5,  fontfamily = "sans", cat.fontface = "italic", cat.cex = 0.8,  cat.default.pos = "outer",
             cat.pos = c(-27, 27),  cat.dist = c(0.055, 0.055),  cat.fontfamily = "sans",  cat.col = c("black", 'black'))

venn.diagram(x = list(vec7, vec8), category.names = c("B. belcheri" , "B. floridae"),
             filename = 'common GO_MF.png', main = "Molecular function",
             main.fontfamily = "sans", main.col = "blue", main.cex = 1.5,
             output=TRUE, imagetype="png" ,  height = 650,   width = 650, resolution = 300,  compression = "lzw",
             lwd = 2,  col=c("#FF6666", '#33CCCC'),
             fill = c(alpha("#FF6666",0.3), alpha('#33CCCC',0.3)),
             cex = 1.5,  fontfamily = "sans", cat.fontface = "italic", cat.cex = 0.8,  cat.default.pos = "outer",
             cat.pos = c(-20, 20),  cat.dist = c(0.055, 0.055),  cat.fontfamily = "sans",  cat.col = c("black", 'black'))

############ DISPENSABILITY OF SELECTED GENES ######################################
freq <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

IDEN_Bb <- read.csv2("Bb.all.csv") # LOAD THE ANNOTATION FILE AS PREPARED BY INTERPROSCAN
IDEN_Bb_count <- IDEN_Bb %>% count(annotation)

IDEN_PAV_Bb <- inner_join(Bb_PAV_sel, IDEN_Bb, by = c("Gene"))
IDEN_PAV_Bb_count <- IDEN_PAV_Bb %>% count(annotation)

PERF_Bb <-IDEN_PAV_Bb[IDEN_PAV_Bb$descr == "MAC/Perforin domain" |IDEN_PAV_Bb$descr == "C-terminal domain of apextrin" | 
                     IDEN_PAV_Bb$descr =="RDRP" | IDEN_PAV_Bb$descr =="Death domain" |
                     IDEN_PAV_Bb$descr =="Lectin C-type domain" | IDEN_PAV_Bb$descr =="Big defensin" |
                     IDEN_PAV_Bb$descr == "TIR domain" |  IDEN_PAV_Bb$descr == "AIG1 family",]
PERF_Bb <- PERF_Bb[!duplicated(PERF_Bb$Gene), ]
PERF_Bb$Species <- "B. belcheri"

#vec_PERF <- IDEN_PAV[IDEN_PAV$descr== "MAC/Perforin domain",]
#PERF_matrix <- inner_join(vec_PERF, bb_gene, by = c("Gene"))
#PERF_count <- inner_join(vec_PERF, bb_PAV_sel, by = c("Gene"))

IDEN_Bf <- read.csv2("Bf_all.csv")
IDEN_Bf_count <- IDEN_Bf %>% count(annotation)

IDEN_PAV_Bf <- inner_join(Bf_PAV_sel, IDEN_Bf, by = c("Gene"))
IDEN_PAV_Bf_count <- IDEN_PAV_Bf %>% count(annotation)

PERF_Bf <-IDEN_PAV_Bf[IDEN_PAV_Bf$descr == "MAC/Perforin domain" | IDEN_PAV_Bf$descr == "C-terminal domain of apextrin" | IDEN_PAV_Bf$descr == "Death domain" | 
                        IDEN_PAV_Bf$descr == "RDRP" | IDEN_PAV_Bf$descr == "Death domain" |
                        IDEN_PAV_Bf$descr == "Lectin C-type domain" | IDEN_PAV_Bf$descr =="Big defensin" |
                        IDEN_PAV_Bf$descr == "TIR domain" |  IDEN_PAV_Bf$descr == "AIG1 family", ]
PERF_Bf <- PERF_Bf[!duplicated(PERF_Bf$Gene), ]

PERF_Bf$Species <- "B. floridae"
PERF_Bf$C <- 41
PERF_Bb$C <- 91
PERF <- rbind(PERF_Bf, PERF_Bb)

library(dplyr)
PERF1 <- PERF %>% 
  group_by(Species) %>%  mutate(freq = n/C*100) %>%
  ungroup

PERF_FIG_bef <- ggboxplot(PERF1, y = "freq", x= "descr",  color = "Species", add = "jitter",  legend = "bottom", lwd=1) +
  theme_minimal() + ylab("Disp. (%)") + theme(text = element_text(size = 20), axis.text.x = element_text(size = 10))+
  theme(legend.position="none", plot.title = element_text(size=14)) +
  rotate_x_text(angle = 0) + xlab("") +
  scale_x_discrete(labels = c("MAC/Perforin domain" = "MAC/Perf", "C-terminal domain of apextrin" = "Apex C_term", "Death domain" = "Death domain",
                              "RDRP" = "RDRP", "Lectin C-type domain" = "C-type lectin", "Big defensin" = "Big defensin",
                              "TIR domain" = "TIR", "AIG1 family" = "AIG")) 
PERF_FIG_bef
PERF_FIG <- ggarrange(ggplot() + theme_void(),PERF_FIG_bef,
                    ncol = 2, nrow = 1, labels = c("g",""), widths = c(0.05,1), 
                    common.legend = TRUE, legend = c("none"), font.label = list(size = 20, color = "black"))
PERF_FIG
VOL_up <- ggarrange(Bb_vol_f,Bf_vol_f,ggplot() + theme_void(), 
                    Bb_GO_f,Bf_GO_f,ggplot() + theme_void(),
                    ncol = 3, nrow = 2, labels = c("auto"), widths = c(1,1,0.5), 
                 heights = c(1,1), common.legend = TRUE, legend = c("none"), font.label = list(size = 20, color = "black"))
VOL_up

##################################### DE-NOVO & FIGURE 2##########################

Bb_denovoPFAM <- read.csv2("Bb_Pfam_denovo.csv") # LOAD ANNOTATIONS FOR THE DE-NOVO PAV GENES
Bf_denovoPFAM <- read.csv2("Bf_Pfam_denovo.csv")

Bb_denovoPFAM_count <- Bb_denovoPFAM %>% count(Domain)
Bf_denovoPFAM_count <- Bf_denovoPFAM %>% count(Domain)

Bb_denovoPFAM_count_top20 <- top_n(Bb_denovoPFAM_count, 20, n)
Bf_denovoPFAM_count_top20 <- top_n(Bf_denovoPFAM_count, 20, n)

denovo_Pfam_Bb <- ggplot(Bb_denovoPFAM_count_top20, aes(x=reorder(Domain,-n), y=n)) + 
  geom_bar(stat="identity", lwd=1, color= "#FF6666", fill="transparent") + xlab("Sample ID")+
  theme_minimal()+ ylab("No. of domains") +
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 14))+
  scale_fill_brewer(palette = "Spectral") + xlab("Pfam domain")+
  rotate_x_text(angle = 90) 
denovo_Pfam_Bb
denovo_Pfam_Bf <- ggplot(Bf_denovoPFAM_count_top20, aes(x=reorder(Domain,-n), y=n)) + 
  geom_bar(stat="identity", lwd=1, color= "#33CCCC", fill="transparent") + xlab("Sample ID")+
  theme_minimal()+ ylab("No. of domains") +
  theme(text = element_text(size = 18), axis.text.x = element_text(size = 14))+
  scale_fill_brewer(palette = "Spectral") + xlab("Pfam domain")+
  rotate_x_text(angle = 90) 
denovo_Pfam_Bf

FIG2_bottom <- ggarrange(ggplot() + theme_void(),denovo_Pfam_Bb, ggplot() + theme_void(), denovo_Pfam_Bf,
                         ncol = 4, nrow = 1, labels = c("h","", "i", ""), widths = c(0.1,1,0.1,1),
                         common.legend = TRUE, legend = c("none"), font.label = list(size = 20, color = "black"))
FIG2_bottom


FIG2 <- ggarrange(VOL_up, PERF_FIG, FIG2_bottom, ncol = 1, nrow = 3, labels = c("","g",""),
                    heights = c(3,1,1.5), common.legend = TRUE, legend = c("none"), font.label = list(size = 20, color = "black"))
FIG2

################################## FIGURE 5 ###########################################
tab <- read.csv2("vecRES_tab.csv") # LOAD COVERAGES OF BELCHERIHV-1 COMPUTED USINGE GENES AND REGIONS AND THEIR RATIOS
tab2 <- read.csv2("SNP_freq.csv") # LOAD THE FREQUENCIES OF SNPs per GENOME nt
tab3 <- read.csv2("BelcheriHV_SNP.csv") # LOAD THE FREQUENCIES OF SNPs FOR BELCHERIHV1

#FIG5b <- ggplot(tab, mapping = aes(x = region, y = gene)) +
#  sm_statCorr(fill = '#0f993d', linetype = "dotted", color = "#0f993d", corr_method = 'spearman')+
#   theme_minimal() + theme(legend.position="none") +
#  geom_point(shape = 21, size = 5, aes(fill="black")) + theme(text = element_text(size = 16))+
#  xlab("Cov. BelcheriHV-1/Bf (region)")+   ylab("Cov. BelcheriHV-1/Bf (genes)") + xlim(0.5,1.4) + ylim(0.5,1.4)

FIG5b <-ggpaired(tab2, x = "type", y = "SNP.nt", color = "type", line.color = "gray", line.size = 0.4, lwd = 5, 
                 palette = "npg")+ stat_compare_means(paired = TRUE, label.y = 4.5) +
  theme_minimal() + ylab("Polymorphism rate (nt / SNP)") + xlab("") + theme(text = element_text(size = 20))+
  theme(legend.position="none", plot.title = element_text(size=20))
FIG5b <- FIG5b + scale_y_continuous(trans='log10')
FIG5b

snp_counts <- tab3 %>% group_by(Sample) %>% summarise(SNP_count = n())
FIG5c <- ggplot(tab3, aes(x = Frequency, y = Sample, fill = Sample)) +
  geom_density_ridges(scale = 1.2, alpha = 0.8, color = "white", rel_min_height = 0.1) +
  geom_text(data = snp_counts, aes(x = 20,  y = Sample,label = paste0("n = ", SNP_count)), 
            hjust = 0, vjust = 0, size = 4, inherit.aes = FALSE) +
  theme_minimal() +  xlab("SNP frequency (%)") + xlim(20,100) +   ylab("") +
  theme(text = element_text(size = 20),  axis.text.y = element_text(size = 14),  legend.position = "none",
        plot.title = element_text(size = 20))
FIG5c

FIG5d <-ggpaired(tab, x = "type", y = "coverage", color = "type", line.color = "gray", line.size = 0.4, lwd = 5, 
                 palette = c("#8A2BE2", "#B8860B"))+ ylim(0.5,1.25)+ 
  theme_minimal() + ylab("Coverage (BelcheriHV-1 / Bf)") + xlab("") + theme(text = element_text(size = 20))+
  theme(legend.position="none", plot.title = element_text(size=20))
FIG5d

Figure5_left <- ggarrange(ggplot() + theme_void(),ggplot() + theme_void(),
                     ggplot() + theme_void(),ggplot() + theme_void(), ggplot() + theme_void(),FIG5b, 
                     ncol = 2, nrow = 3,labels = c("a","","","", "b",""), heights = c(2,0.1,1.2), widths = c(0.05,1),
                     common.legend = TRUE, legend = c("none"), font.label = list(size = 20, color = "black"))
Figure5_left
Figure5_right <- ggarrange(ggplot() + theme_void(),FIG5c,
                           ggplot() + theme_void(),ggplot() + theme_void(), ggplot() + theme_void(),FIG5d, 
                           ncol = 2, nrow = 3,labels = c("c","","","", "d",""), heights = c(2,0.1,1.2), widths = c(0.05,1),
                           common.legend = TRUE, legend = c("none"), font.label = list(size = 20, color = "black"))
Figure5_right

Figure5 <- ggarrange(Figure5_left, Figure5_right, 
                     ncol = 2, nrow = 1,labels = c("",""), widths = c(1,1), 
                     common.legend = TRUE, legend = c("none"), font.label = list(size = 24, color = "black"))
Figure5

##################### Expression analysis of PAV genes (B. belcheri) ###############

Bb_RNA <- read.csv2("Bb_RNA.csv")
Bb_RNA_t <- read.csv2("Bb_RNA_t.csv")
Bb_dispensable <- Bb_PAV_count["Gene"]

Bb_RNA_PAV <- inner_join(Bb_PAV_sel, Bb_RNA, by = c("Gene"))
Bb_RNA_noPAV <- anti_join(Bb_RNA, Bb_dispensable, by = "Gene")
Bb_a <- melt(Bb_RNA_PAV, id.vars="Gene")
Bb_b <- Bb_a %>%  filter(!variable=='n')
Bb_b$cond <- "shell/cloud"

Bb_RNA_PAV$Total_0s<-rowSums(Bb_RNA_PAV==0)
Bb_RNA_noPAV$Total_0s<-rowSums(Bb_RNA_noPAV==0)
mean(Bb_RNA_noPAV$Total_0)
mean(Bb_RNA_PAV$Total_0)

Bb_zero_PAV <- select(Bb_RNA_PAV, c(1, 94))
Bb_zero_PAV$cond <- "shell/cloud"
Bb_zero <- select(Bb_RNA_noPAV, c(1, 93))
Bb_zero$cond <- "core"
Bb_double_zero <- rbind(Bb_zero_PAV, Bb_zero)

Bb_c <- melt(Bb_RNA_t, id.vars="Gene")
Bb_c$cond <- "core"
Bb_d <- rbind(Bb_c, Bb_b)
Bb_cld <- Bb_d %>% group_by(variable) %>% summarise(pvalue = t.test(value[cond=="shell/cloud"], value[cond=="core"], paired= FALSE)$p.value)

geni_virus <- Bb_PAV_sel  %>%  filter(n=='80')
RNA_virus <- inner_join(geni_virus, Bb_RNA, by = c("Gene"))
e <- melt(RNA_virus, id.vars="Gene")
e$cond <- "BbHV"
e <- e %>%  filter(!variable=='n')
e <- e %>%  filter(!variable=='Total_0s')
f <- rbind(Bb_c,e)

f$value <- replace(f$value, f$value < 3, NA)

FIG_RNA_virus <- ggboxplot(f, y = "value", x= "variable",  color = "cond", legend = "bottom", lwd=1) +
  theme_minimal() + ylab("Expression level (TPMs)") + xlab("") + theme(text = element_text(size = 14))+
  theme(legend.position="bottom", plot.title = element_text(size=14)) +
  rotate_x_text(angle = -90) 
FIG_RNA_virus <- FIG_RNA_virus + scale_y_continuous(trans='log10')
FIG_RNA_virus

FIG_RNA_virus_gene <- ggboxplot(e, y = "value", x= "Gene",  color = "lightsalmon", 
                           add = "jitter",    legend = "bottom", lwd=1) +
  theme_minimal() + ylab("Expression level (TPMs)") + xlab("") + theme(text = element_text(size = 14))+
  theme(legend.position="none", plot.title = element_text(size=14)) +
  rotate_x_text(angle = -90) +
  geom_hline(yintercept=5, col="black", linetype="dashed", lwd=1)
FIG_RNA_virus_gene <- FIG_RNA_virus_gene + scale_y_continuous(trans='log10')
FIG_RNA_virus_gene

############################### Expression analysis of PAV genes (B. floridae) #####

Bf_RNA <- read.csv2("Bf_RNA.csv")
Bf_RNA_PAV <- inner_join(Bf_PAV_sel, Bf_RNA, by = c("Gene"))
Bf_a <- melt(Bf_RNA_PAV, id.vars="Gene")
Bf_b <- Bf_a %>%  filter(!variable=='n')
Bf_b$cond <- "shell/cloud"
Bf_dispensable <- Bf_PAV_count["Gene"]
Bf_RNA_noPAV <- anti_join(Bf_RNA, Bf_dispensable, by = "Gene")

Bf_RNA_PAV$Total_0s<-rowSums(Bf_RNA_PAV==0)
Bf_RNA_noPAV$Total_0s<-rowSums(Bf_RNA_noPAV==0)
mean(Bf_RNA_noPAV$Total_0)
mean(Bf_RNA_PAV$Total_0)

Bf_zero_PAV <- select(Bf_RNA_PAV, c(1, 64))
Bf_zero_PAV$cond <- "shell/cloud"
Bf_zero <- select(Bf_RNA_noPAV, c(1, 63))
Bf_zero$cond <- "core"
Bf_double_zero <- rbind(Bf_zero_PAV, Bf_zero)
Bf_double_zero$Species <- "B. floridae"
Bb_double_zero$Species <- "B. belcheri"

Bb_double_zero$cond <- factor(Bb_double_zero$cond, levels = c("core", "shell/cloud"))
Bf_double_zero$cond <- factor(Bf_double_zero$cond, levels = c("core", "shell/cloud"))

Bf_FIG_zero <- ggboxplot(Bf_double_zero, y = "Total_0s", x= "cond",  color = "#33CCCC", legend = "bottom", lwd=1) +
  theme_minimal() + ylab("Count of '0'") + xlab("") + theme(text = element_text(size = 14))+
  theme(legend.position="none", plot.title = element_text(size=14)) +
  rotate_x_text(angle = -90) +
  stat_compare_means(label = "p.signif", method = "t.test",label.x = 1.5, size = 10, label.y = 50) 
Bf_FIG_zero
Bb_FIG_zero <- ggboxplot(Bb_double_zero, y = "Total_0s", x= "cond",  color = "#FF6666", legend = "bottom", lwd=1) +
  theme_minimal() + ylab("Count of '0'") + xlab("") + theme(text = element_text(size = 14))+
  theme(legend.position="none", plot.title = element_text(size=14)) +
  rotate_x_text(angle = -90) +
  stat_compare_means(label = "p.signif", method = "t.test",label.x = 1.5, size = 10, label.y = 75) 
Bb_FIG_zero

Bf_c <- melt(Bf_RNA, id.vars="Gene")
Bf_c$cond <- "core"

Bf_b <- Bf_b %>%  filter(!variable=='Total_0s')
Bf_c <- Bf_c %>%  filter(!variable=='Total_0s')

Bf_d <- rbind(Bf_c, Bf_b)
Bf_cld <- Bf_d %>% group_by(variable) %>% summarise(pvalue = t.test(value[cond=="PAV"], value[cond=="total"], paired= FALSE)$p.value)

Bf_d$Species <- "B. floridae"
Bb_d$Species <- "B. belcheri"
d <- rbind(Bf_d, Bb_d)

Bb_FIG_RNA_tot <- ggboxplot(Bb_d, y = "value", x= "cond",  color = "#FF6666", legend = "none", lwd=1) +
  theme_minimal() + ylab("Expression level (TPMs)") + xlab("") + theme(text = element_text(size = 14))+
  theme(legend.position="none", plot.title = element_text(size=14)) +
  rotate_x_text(angle = -90) +
  stat_compare_means(label = "p.signif", method = "t.test", size = 10, label.x = 1.5, label.y = 4)  
Bb_FIG_RNA_tot <- Bb_FIG_RNA_tot + scale_y_continuous(trans='log10')
Bb_FIG_RNA_tot

Bf_FIG_RNA_tot <- ggboxplot(Bf_d, y = "value", x= "cond",  color = "#33CCCC", legend = "none", lwd=1) +
  theme_minimal() + ylab("Expression level (TPMs)") + xlab("") + theme(text = element_text(size = 14))+
  theme(legend.position="none", plot.title = element_text(size=14)) +
  rotate_x_text(angle = -90) +
  stat_compare_means(label = "p.signif", method = "t.test", size = 10, label.x = 1.5, label.y = 4) 
Bf_FIG_RNA_tot <- Bf_FIG_RNA_tot + scale_y_continuous(trans='log10')
Bf_FIG_RNA_tot

########################### PAV ANALYSIS PARENT --> OFFSPRING #####################

SRR6162896_F <- read.csv2("SRR6162896_mapping_depthgenes.coverage.csv", sep = ",", dec = ".")
mean_cov_F = c(mean(SRR6162896_F$X0, na.rm=TRUE))
F_cov_PAV <- inner_join(SRR6162896_F, Bf_PAV_sel_d, by = c("Gene"))
F_cov_noPAV <- inner_join(SRR6162896_F, Bf_genes, by = c("Gene"))
mean_cov_F_nopav = c(mean(F_cov_noPAV$X0, na.rm=TRUE))

F_cov_low <- inner_join(SRR6162896_F, Bf_PAV_sel_d50, by = c("Gene"))
mean_cov_F_cov_low = c(mean(F_cov_low$X0, na.rm=TRUE))

F_cov_fig <- ggplot(F_cov_low, aes(x=X0))  + geom_density(alpha=0.5, fill="blue")  +
  theme_minimal()+   theme(text=element_text(size=14),axis.text = element_text(size = 14)) + 
  theme(legend.position = "bottom") + xlim(0,75)+ xlab("Coverage") +
  geom_vline(xintercept = 7.5, linetype="dashed", color = "violet", lwd =1) +
  geom_vline(xintercept = 15, linetype="dashed", color = "blue", lwd =1) +
  geom_vline(xintercept = 45, linetype="dashed", color = "blue", lwd =1) +
  geom_vline(xintercept = 60, linetype="dashed", color = "orange", lwd =1)
F_cov_fig
F_cov <- ggplot(SRR6162896_F, aes(x=X0))  + geom_density(alpha=0.5, fill="blue") +
  theme_minimal()+   theme(text=element_text(size=14),axis.text = element_text(size = 14)) + 
  theme(legend.position = "bottom") + xlim(0,75)+ xlab("Coverage")+
  geom_vline(xintercept = 60, linetype="dashed", color = "orange", lwd =1)
F_cov

f_absent <- F_cov_PAV[F_cov_PAV$X0<7.5,]
f_present <- F_cov_PAV[F_cov_PAV$X0>45,]
f_hemy <- F_cov_PAV[F_cov_PAV$X0>15 & F_cov_PAV$X0<45 ,]

SRR6162897_M <- read.csv2("SRR6162897_mapping_depthgenes.coverage.csv", sep = ",", dec = ".")
mean_cov_M = c(mean(SRR6162897_M$X0, na.rm=TRUE))
M_cov_PAV <- inner_join(SRR6162897_M, Bf_PAV_sel_d, by = c("Gene"))
M_cov_noPAV <- inner_join(SRR6162897_M, Bf_genes, by = c("Gene"))
mean_cov_M_nopav = c(mean(M_cov_noPAV$X0, na.rm=TRUE))

M_cov_low <- inner_join(SRR6162897_M, Bf_PAV_sel_d50, by = c("Gene"))
mean_cov_M_cov_low = c(mean(M_cov_low$X0, na.rm=TRUE))

M_cov_fig <- ggplot(M_cov_low, aes(x=X0))  + geom_density(alpha=0.5, fill="green")  + 
  theme_minimal()+   theme(text=element_text(size=14),axis.text = element_text(size = 14)) + 
  theme(legend.position = "bottom") + xlim(0,75)+ xlab("Coverage") +
  geom_vline(xintercept = 7, linetype="dashed", color = "violet", lwd =1) +
  geom_vline(xintercept = 14, linetype="dashed", color = "blue", lwd =1) +
  geom_vline(xintercept = 42, linetype="dashed", color = "blue", lwd =1) +
  geom_vline(xintercept = 56, linetype="dashed", color = "orange", lwd =1)
M_cov_fig

M_cov <- ggplot(SRR6162897_M, aes(x=X0))  + geom_density(alpha=0.5, fill="green") +
  theme_minimal()+   theme(text=element_text(size=14),axis.text = element_text(size = 14)) + 
  theme(legend.position = "bottom") + xlim(0,75)+ xlab("Coverage")+
  geom_vline(xintercept = 56, linetype="dashed", color = "orange", lwd =1)
M_cov

m_absent <- M_cov_PAV[M_cov_PAV$X0<7,]
m_present <- M_cov_PAV[M_cov_PAV$X0>42,]
m_hemy <- M_cov_PAV[M_cov_PAV$X0>14 & M_cov_PAV$X0<42,]

FIGS5 <- ggarrange(F_cov, F_cov_fig, M_cov, M_cov_fig, ncol = 2, nrow = 2, labels = c("auto"), 
                     common.legend = TRUE, legend = c("none"), font.label = list(size = 20, color = "black"), hjust = 0)
FIGS5

venn.diagram(
  x = list(vec1, vec2),
  category.names = c("parent/off" , "PAV"),
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff'))

NOPAV <- inner_join(M_cov_noPAV, F_cov_noPAV, by = c("Gene"))
ALWAYS_present <- inner_join(f_present, m_present, by = c("Gene"))
ALWAYS_absent <- inner_join(f_absent, m_absent, by = c("Gene"))
ALWAYS_hemy <- inner_join(f_hemy, m_hemy, by = c("Gene"))
Mpres_Fabs <- inner_join(m_present, f_absent, by = c("Gene"))
Fpres_Mabs <- inner_join(f_present, m_absent, by = c("Gene"))
pres_abs <- merge(Mpres_Fabs, Fpres_Mabs, by = "Gene", all = TRUE)
F_HEMY_Mabsent <- inner_join(f_hemy, m_absent, by = c("Gene"))
M_HEMY_Fabsent <- inner_join(m_hemy, f_absent, by = c("Gene"))
HEMY_ABS <- merge(F_HEMY_Mabsent, M_HEMY_Fabsent, by = "Gene", all = TRUE)

F_HEMY_Mpres <- inner_join(f_hemy, m_present, by = c("Gene"))
M_HEMY_Fpres <- inner_join(m_hemy, f_present, by = c("Gene"))
HEMY_PRES <- merge(F_HEMY_Mpres, M_HEMY_Fpres, by = "Gene", all = TRUE)

Bf_doughter <- read.csv2("Bf_matrix_genes_doughter.csv")
Bf_10_d <- inner_join(Bf_doughter, NOPAV, by = c("Gene"))
row.names(Bf_10_d) <- Bf_10_d$Gene
Bf_s_d <- select(Bf_10_d, c(2, 2:97))
Bf_doughter <- data.matrix(Bf_s_d)

#heatmaply(Bf_doughter, colors = Blues, k_row = 12, k_col = 2, cexRow = 0.8, margins = c(100,100,100,100), 
#          fontsize_row =1, fontsize_col = 8, main="lupin PAV", xlab = "Samples", ylab = "Dispensable Genes", 
#          seriate = "OLO", plot_method = "plotly", file = "bf_gene_matrix_d.html")

frac <- read.csv2("fraction.csv")  # LOAD THE FRACTIONS OF GENES IN PAV FOR THE DIFFERENT CLASSES (TABLE 2)
#my_comparisons <- list( c("II", "I"),  c("II", "III") )
FIG8 <- ggviolin(frac, x = "Class", y = "perc",  add = "boxplot", 
                  xlab = "PAV gene classes", ylab = "Absent in progeny (%)",
                  short.panel.labs = TRUE,  size = 1, color = "black", fill = "white") + theme_minimal() +
                  theme(legend.position="none", text = element_text(size = 14)) +
                  stat_compare_means(method = "anova",label.y = 95, label.x = 2.5)
FIG8

##################################### ORTHOLOGY ANALYSIS ##########################

Bf_ortho_all <- read.csv2("Bf_selected_2.csv")
Bb_ortho_all <- read.csv2("Bb_selected_2.csv")
lancelet_ortho_all <- read.csv2("lancelet_selected_2.csv")
Orthogroups <- read.csv2("Orthogroups.csv")
Orthogroups_sum <- read.csv2("Orthogroups.GeneCount.csv")

Bf_ORTHO_count <- Bf_ortho_all %>% count(Orthogroup)
Bf_CONTE <- left_join(Bf_ORTHO_count, Orthogroups_sum, by = c("Orthogroup"))
Bb_ORTHO_count <- Bb_ortho_all %>% count(Orthogroup)
Bb_CONTE <- left_join(Bb_ORTHO_count, Orthogroups_sum, by = c("Orthogroup"))

Bf_duplicated_rows <- duplicated(Bf_ortho_all)
Bb_duplicated_rows <- duplicated(Bb_ortho_all)
lancelet_duplicated_rows <- duplicated(lancelet_ortho_all)

# Extracting unique rows from the matrix 
Bf_ortho_unique <- Bf_ortho_all[!Bf_duplicated_rows, ]
Bb_ortho_unique <- Bb_ortho_all[!Bb_duplicated_rows, ]
lancelet_ortho_unique <- lancelet_ortho_all[!lancelet_duplicated_rows, ]

lancelet_ortho_unique_number <- left_join(lancelet_ortho_unique, Orthogroups_sum, by = c("Orthogroup"))
#write.table(lancelet_ortho_unique_number, "lancelet_ortho_unique_number.csv", sep = ";",
#            dec = ",", col.names = NA)

ortho <- read.csv2("ortho_count_2.csv")
ortho_PAV <- read.csv2("lancelet_ortho_unique_number.csv")

upset(ortho, nsets = 13, point.size = 3.5, line.size = 1.5, nintersects = 10, 
      mb.ratio = c(0.5, 0.5),
      mainbar.y.label = "No. of orthologs", sets.x.label = "No. of species orthologs", 
      text.scale = c(2, 2, 1, 1, 1.5, 2), order.by = "freq")

Bb_ortho_unique_ID <- as.vector(Bb_ortho_unique$Orthogroup)
Bf_ortho_unique_ID <- as.vector(Bf_ortho_unique$Orthogroup)

venn.diagram(x = list(Bb_ortho_unique_ID, Bf_ortho_unique_ID), category.names = c("B. belcheri" , "B. floridae"),
             filename = 'common ORTHO.png', main = "Orthogroups",
             main.fontfamily = "sans", main.col = "black", main.cex = 1.5,
             output=TRUE, imagetype="png" ,  height = 900 ,   width = 900 ,   resolution = 300,  compression = "lzw",
             lwd = 2,  col=c("#FF6666", '#33CCCC'), fill = c(alpha("#FF6666",0.3), alpha('#33CCCC',0.3)),
             cex = 1,  fontfamily = "sans", cat.fontface = "italic", cat.cex = 0.8,  cat.default.pos = "outer",
            cat.pos = c(-27, 27),  cat.dist = c(0.055, 0.055),  cat.fontfamily = "sans",  cat.col = c("black", 'black'))

Bb_CONTE$fraction <- Bb_CONTE$n / Bb_CONTE$Bb
Bf_CONTE$fraction <- Bf_CONTE$n / Bf_CONTE$Bf
Bb_ortho_fraction <- as.vector(Bb_CONTE$fraction)
Bf_ortho_fraction <- as.vector(Bf_CONTE$fraction)
plot(Bb_ortho_fraction)

bb_fraction <- select(Bb_CONTE, c('fraction'))
bb_fraction$Species <- "Bb"
bf_fraction <- select(Bf_CONTE, c('fraction'))
bf_fraction$Species <- "Bf"
fraction <- rbind(bb_fraction, bf_fraction)
Bb_PAV_fraction <- ggviolin(fraction, x = "Species", y = "fraction",  add = "boxplot", 
                            xlab = "", ylab = "Fraction of dispensable genes per orthogroup",
                            short.panel.labs = TRUE,  size = 1, color = "Species", fill = "white") + theme_minimal() +
  theme(legend.position="none", text = element_text(size = 12)) 
Bb_PAV_fraction


FIG3_down <- ggarrange(ggplot() + theme_void(), ggplot() + theme_void(), ggplot() + theme_void(), Bb_PAV_fraction,  ncol = 4, nrow = 1, labels = c("b","","c",""), 
                  widths = c(0.1,1,0.1,1), common.legend = TRUE, 
                  legend = c("none"), font.label = list(size = 20, color = "black"), hjust = 0) 
FIG3_mid <- ggarrange(ggplot() + theme_void(), FIG3_down,  ncol = 1, nrow = 2, labels = c("a",""), 
                  heights = c(2,1), common.legend = TRUE, 
                  legend = c("none"), font.label = list(size = 20, color = "black"), hjust = 0) 
FIG3 <- ggarrange(FIG3_mid, ggplot() + theme_void(), ncol = 1, nrow = 2, labels = c("","d"), 
                          heights = c(2,1), common.legend = TRUE, 
                          legend = c("none"), font.label = list(size = 20, color = "black"), hjust = 0) 
FIG3

ortho_PAV_denovo <- read.csv2("Orthogroups.GeneCount_denovo.csv")
upset(ortho_PAV_denovo, nsets = 13, point.size = 3.5, line.size = 1.5, nintersects = 10, 
      mb.ratio = c(0.5, 0.5), matrix.color = "#0f993d", 
      main.bar.color = c("#0f993d"), 
      sets.bar.color = "#0f993d", 
      mainbar.y.label = "No. of orthologs", sets.x.label = "No. of species orthologs", 
      text.scale = c(2, 2, 1, 1, 1.5, 2), order.by = "freq")

##################################### FIGURE 1 ######################################
data_forFIG1a <- data[!data$Run %in% c("SRR2757687", "SRR2757688", "SRR12010277"), ]
Fig1a <- ggboxplot(data_forFIG1a, y = "No_PAV_gene", color = "Species", palette = c("#FF6666", "#33CCCC"),
                   add = "jitter", legend = "bottom") +
  theme_minimal() + ylab("No . of absent genes") + xlab("") + theme(text = element_text(size = 14))+
  theme(legend.position="bottom", plot.title = element_text(size=14)) + 
  rotate_x_text(angle = 0) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
Fig1a

data_bb <- subset(data, Species == "Bb")
data_bf <- subset(data, Species == "Bf")

Fig1b_Bb <- ggplot(data_bb, mapping = aes(x = No_PAV_region, y = No_PAV_gene, fill=Species)) +
  geom_point(shape = 21, size = 5) +
  sm_statCorr(corr_method = 'spearman')+ scale_fill_manual(values = c("#FF6666"))+
  theme_minimal() + theme(legend.position = "none", text = element_text(size = 14), 
  plot.title = element_text(size = 14, hjust = 0.5, face = "italic")) +
  theme(text = element_text(size = 14))+ theme(plot.title = element_text(size=14)) +
  xlab("No. of absent fragments (5kb)")+  ylab("No. of absent genes") + xlim(1000,2000) + ylim(300,650) 
Fig1b_Bb
Fig1b_Bf <- ggplot(data_bf, mapping = aes(x = No_PAV_region, y = No_PAV_gene, fill=Species)) +
  geom_point(shape = 21, size = 5) +
  sm_statCorr(corr_method = 'spearman')+ scale_fill_manual(values = c("#33CCCC"))+
  theme_minimal() + theme(legend.position = "none", text = element_text(size = 14), 
  plot.title = element_text(size = 14, hjust = 0.5, face = "italic")) +
  theme(text = element_text(size = 14))+ theme(plot.title = element_text(size=14)) +
  xlab("No. of absent fragments (5kb)")+  ylab("No. of absent genes") + xlim(3750,4500) + ylim(350,650) 
Fig1b_Bf

FIG1_sopra <- ggarrange(ggplot() + theme_void(),Fig1a, ggplot() + theme_void(),Fig1b_Bb,Fig1b_Bf, ncol = 5, nrow = 1, labels = c("a","","b","",""), 
                        widths = c(0.05,0.6,0.05,0.6,0.6), common.legend = TRUE, heights = c(1,1),
                        legend = c("right"), font.label = list(size = 20, color = "black"), hjust = 0) 
FIG1_sopra
FIG1_sotto <- ggarrange(ggplot() + theme_void(), Bb_FIG_RNA_tot,ggplot() + theme_void(), Bf_FIG_RNA_tot, ggplot() + theme_void(),Bb_FIG_zero, ggplot() + theme_void(),Bf_FIG_zero,
                        ncol = 8, nrow = 1, labels = c("c","","d","","e","","f"), 
                        widths = c(0.05,1,0.05,1,0.05,1,0.05,1),common.legend = TRUE, legend = c("none"), font.label = list(size = 20, color = "black"), hjust = 0) 
FIG1_sotto
FIG1 <- ggarrange(FIG1_sopra, ggplot() + theme_void(),FIG1_sotto, ncol = 1, nrow = 3, labels = c(""), heights = c(2,0.05,1.5),
                  common.legend = TRUE, legend = c("none"), hjust = 0)
FIG1

##################### FIG S1 - DISPENSABILITY LEVEL################################
Bb_disp <- read.csv2("bb_1c.csv")
Fig1c <- ggplot(Bb_disp, mapping = aes(x = dispensability, y = nn)) +
  sm_statCorr(fill = '#0f993d', linetype = "dotted", color = "#0f993d", corr_method = 'spearman')+
  theme(legend.position="bottom") +  theme_minimal() +
  geom_point(shape = 18, size = 3, color="#FF6666") + theme(text = element_text(size = 14))+ ylab("No. of occurences")+
  xlab("Dispensability (% of individuals)")
Fig1c <- Fig1c + scale_y_continuous(trans='log10')
Fig1c

Bf_disp <- read.csv2("bb_1d.csv")
Fig1d <- ggplot(Bf_disp, mapping = aes(x = dispensability, y = nn)) +
  sm_statCorr(fill = '#0f993d', linetype = "dotted", color = "#0f993d", corr_method = 'spearman')+
  theme(legend.position="bottom") +  theme_minimal() +
  geom_point(shape = 18, size = 3, color="#33CCCC") + theme(text = element_text(size = 14))+ ylab("No. of occurences")+
  xlab("Dispensability (% of individuals)")
Fig1d <- Fig1d + scale_y_continuous(trans='log10')
Fig1d
FIGS1 <- ggarrange(ggplot() + theme_void(),Fig1c, ggplot() + theme_void(),Fig1d,
                   ncol = 4, nrow = 1, labels = c("a","","b",""), widths = c(0.05,1,0.05,1), 
                   common.legend = TRUE, legend = c("none"), font.label = list(size = 20, color = "black"))
FIGS1
