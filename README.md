# HNRNPC_OSCC
Code of HNRNPC drives ITGA1+ myofibroblast-dependent T cell exhaustion via metabolic reprogramming in oral squamous cell carcinoma 
#umap
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(DoubletFinder)
library(harmony)
set.seed(1)
scRNA_harmony  <- NormalizeData(scRNA_harmony ) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution =1.3)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:30)
gene=c("KRT14","KRT16","CDH1",
       "EMCN","INSR","GRB10",
       "RELN","NRXN1","STMN2",
       "COL1A1","COL1A2","COL11A1",
       "DMXL2","MSR1","SLC11A1",
      	'RIPOR2','CD247',"CD2","CD8A","CD8B",
       "ACTA1","MYH1","TPM2")
DotPlot(scRNA_harmony,features =gene,dot.scale=10,col.min = 0)+RotatedAxis()+ 
  scale_color_gradientn(colors = c("white", "#D75455","#C3272B"))+
  theme(axis.text.x = element_text(angle = 90))
scRNA_harmony=RenameIdents(scRNA_harmony,"0"="Epithelial","7"="Epithelial","1"="Epithelial","9"="Epithelial","18"="Epithelial","21"="Epithelial","2"="Epithelial",
                           "3"="Endothelial","8"="Endothelial","14"="Endothelial","13"="Endothelial",
                           "4"="Neuron","17"="Neuron",
                           "5"="Stromal cell","15"="Stromal cell","16"="Stromal cell","6"="Myeloid cell","12"="T cell",
                           "10"="Muscle cell","11"="Muscle cell","19"="Muscle cell","20"="Muscle cell")
scRNA_harmony$maintype=scRNA_harmony@active.ident
DimPlot(scRNA_harmony , reduction = "umap",label = T) 

#cell ratio
library(reshape2)
library(ggplot2)
library(RColorBrewer)
#maintype
col=c("Endothelial"="#E41A1C" ,"Epithelial"="#377EB8" ,"Myeloid cell"="#4DAF4A", "Muscle cell"="#984EA3","Neuron"="#FF7F00","Stromal cell"="#A65628","T cell"="#80B1D3")
pB2_df=read("maintype-group-ratio.csv")
pB4 <- ggplot(data = pB2_df, aes(x =group, y =Number, fill =  celltype)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(y="Ratio",x="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
pB4

#stromal
COLI=c("Fibroblast"="#8CABC5","myCAF1"="#FFCC4F","myCAF2"="#8D5A85","myCAF3"="#e06625","ADSC"="#5CB392","myCAF4"="#F5B7D2","myCAF5"="#C3272B")
pB3_df=write.csv("stromal-group-ratio.csv")
pB5 <- ggplot(data = pB3_df, aes(x =group, y =Number, fill =  celltype)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=COLI) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(y="Ratio",x="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
pB5

#T cell
col=c("Effector_CD8+_T_cell"="#AAB1D9","Exhausted_CD8+_T_cell"="#E09345","CD4+_T_cell"="#F5B7D2")
pB4_df=write.csv("T-cell-celltype-group-ratio.csv"")
pB6 <- ggplot(data = pB4_df, aes(x =group, y =Number, fill =  celltype)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(y="Ratio",x="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45) )
pB6

#heatmap of marker genes
library(Seurat)
library(ggplot2)
markers <- FindAllMarkers(object = stromal, test.use="wilcox" ,only.pos = TRUE, logfc.threshold = 0.25)   
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n =15, wt = avg_log2FC)
stromal2=subset(stromal ,downsample=100)
top10=write.csv("stromal_top10marker.csv")
DoHeatmap(stromal2,features = top10$gene,group.colors = COLI)+ 
  scale_fill_gradientn(colors = c("#377EB8", "white", "#C3272B"))

#monocle2
ibrary(monocle)
table(stromal$group)
Idents(stromal)="group"
hi_mycaf=subset(stromal,idents="HNRNPC-High")
Idents(hi_mycaf)="celltype"
table(hi_mycaf$celltype)
data <- as(as.matrix(hi_mycaf@assays$RNA@counts),'sparseMatrix') 
pd <- new('AnnotatedDataFrame', data = hi_mycaf@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, 
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id    #过滤
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
COLI=c("Fibroblast"="#8CABC5","myCAF1"="#FFCC4F","myCAF2"="#8D5A85","myCAF3"="#e06625",
              "ADSC"="#5CB392","myCAF4"="#F5B7D2","myCAF5"="#C3272B")
plot_cell_trajectory(HSMM, color_by = "celltype") +
  facet_wrap(~celltype, nrow = 3)+#  
  scale_colour_manual(
    values = COLI) )

#GSVA
library('GSEABase')
library(GSVA)
exp=AverageExpression(stromal) 
exp=exp[["RNA"]]
GSVA_hall <- gsva(expr=as.matrix(exp), 
                  gset.idx.list=geneSets, 
                  mx.diff=T, 
                  kcdf="Poisson",
                  parallel.sz=4) 
b=read.csv("myCAF-cluster-gsva")
pheatmap::pheatmap(b,
                   cluster_rows = F,
                   cluster_cols =F,
                   show_colnames=T,
                   scale = "row", 
                   color =colorRampPalette(c("#377EB8","white",  "#C3272B"))(100))
