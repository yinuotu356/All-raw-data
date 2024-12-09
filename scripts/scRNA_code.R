rm(list = ls())
suppressMessages(library(survminer))
suppressMessages(library(clusterProfiler))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(magrittr))
suppressMessages(library(gtools))
suppressMessages(library(stringr))
suppressMessages(library(Matrix))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))
suppressMessages(library(scales))
suppressMessages(library(ggsci))
suppressMessages(library(sctransform))
suppressMessages(library(harmony))
suppressMessages(library(tidydr))
suppressMessages(library(celldex))
suppressMessages(library(pheatmap))
suppressMessages(library(clustree))
suppressMessages(library(xlsx))
suppressMessages(library(gridExtra))
suppressMessages(library(ggthemes))
suppressMessages(library(ggnewscale))
suppressMessages(library(CellChat))
suppressMessages(library(ggpubr))
suppressMessages(library(patchwork))
suppressMessages(library(monocle))
options(stringsAsFactors = F)
colors <- c("#4ea64a","#8e4c99","#e88f18","#e47faf","#b698c5","#a05528","#58a6d6","#1f2d6f","#279772","#add387","#d9b71a","#d5231d","#3777ac")

## 
# 
## 
setwd("00_origin_datas/GEO")
dir <- dir("./")
samples_name = dir
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-") 
  }
}
### 
names(scRNAlist) <- samples_name
sce <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
sum(table(sce@meta.data$orig.ident))##  93431
###############################################################################
saveRDS(sce,file = "01_scRNA/sce.RDS")
##############################################################################################################
###   
scRNA = subset(sce, subset=nFeature_RNA>300&nCount_RNA>3000&percent.mito<25)
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA)
scRNA <- RunPCA(scRNA, features=VariableFeatures(object=scRNA), verbose=FALSE)
#scRNA = SCTransform(scRNA, vars.to.regress="percent.mito", verbose=FALSE)
colnames(scRNA@meta.data)[1] = "Sample"
sum(table(scRNA@meta.data$Sample))##  47550
scRNA = RunHarmony(scRNA, group.by.vars="Sample")
###############################################################################
saveRDS(scRNA,file = "01_scRNA/scRNA.RDS")
#scRNA=readRDS("01_scRNA/scRNA.RDS")
ElbowPlot(scRNA,ndims=30)
scRNA_umap <- RunUMAP(scRNA, dims=1:15, reduction="harmony")
scRNA_umap <- FindNeighbors(scRNA_umap, dims=1:15, reduction="harmony")
scRNA_umap <- FindClusters(scRNA_umap  , resolution=0.2)##这里的resolution=0.1
mydata<-scRNA_umap

scRNA_tsne <- RunTSNE(scRNA, dims=1:15, reduction="harmony")
scRNA_tsne <- FindNeighbors(scRNA_tsne, dims=1:15, reduction="harmony")
scRNA_tsne <- FindClusters(scRNA_tsne  , resolution=0.5)##resolution=0.1

######################################  
p1 = DimPlot(scRNA_umap , reduction="umap", group.by="seurat_clusters", pt.size=1)+theme(legend.position="right", plot.title=element_blank())
#DimPlot(scRNA_tsne , reduction="tsne", label=T,group.by="seurat_clusters", pt.size=1)+theme(legend.position="right", plot.title=element_blank())

p2= VlnPlot(scRNA_umap, features=c("nFeature_RNA", "nCount_RNA","percent.mito"),  pt.size=0)+NoLegend()+#+theme(text=element_text(family="Times"))
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
FigureS1=mg_merge_plot(p2,p1,ncol=2,nrow=1,labels = c('A','B'),widths=c(2,1),common.legend = T,legend ="right")
ggsave('PDFs/Figure S1.pdf',p1,height = 8,width = 8)
ggsave('PDFs/Figure S1.jpg',p1,height = 8,width = 8)


markers <- FindAllMarkers(scRNA_umap, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "01_scRNA/All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
#markers <- read.table("01_scRNA/All_celltype_DEG.txt",header = T,check.names = F,fill=T,sep = "\t")
markers["pct.diff"]=markers$pct.1-markers$pct.2

#################################################################################
### 
Top5 <- markers %>% group_by(cluster) %>% slice_max(n =5, order_by =avg_logFC )
Top51 <- markers %>% group_by(cluster) %>% slice_max(n =5, order_by =pct.diff )
Cellmarker <- read.table("01_scRNA/Cellmarker2.0.txt",header = T,check.names = F,fill=T,sep = "\t")
Cellmarker <- Cellmarker[,c(9,7)]
colnames(Top51)[7]="marker"
Cellmarker2 <- merge(Cellmarker,Top51,by="marker")
write.table(Cellmarker2, "01_scRNA/Cellmarker2.txt", col.names=T, row.names=F, quote=F, sep="\t")

u_plot_clusters=DimPlot(scRNA_umap, reduction="umap",label = T,label.size = 8,group.by="seurat_clusters", pt.size=0.1)+
  theme(legend.position="right", plot.title=element_blank())+
  theme(text=element_text(family="Times"))+
  #scale_color_manual(values=colors)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))
u_plot_clusters

VlnPlot(scRNA_umap, features=c("FGFBP2","GZMH"), pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())




###########
# 0: Natural killer cell:"GZMB","GZMH","FGFBP2","GNLY","CD160",
# 1:Natural killer cell:"GZMB","GZMH","FGFBP2","GNLY","CD160",
# 2:T cell:"CD2","CD3D","CD3E","IL7R"
# 3:T cell:"CD2","CD3D","CD3E","IL7R"
# 4:Macrophage："AIF1","MS4A7","LYZ"
# 5:Plasma cell:"DERL3","IGHG1","MZB1"
# 6:CD8+ T cell""MKI67","MKI67","STMN1"
# 7:B cell:"CD79A","MS4A1"
# 8:Macrophage："AIF1","MS4A7","LYZ"
# 9:Macrophage："AIF1","MS4A7","LYZ"
# 10:Hepatocyte:"ALB","APOA1","KRT18"
# 

cell_label = c(
   "Natural killer cell",
   "Natural killer cell",
  "T cell",
  "T cell",
  "Macrophage",
  "Plasma cell",
 "CD8+ T cell",
 "B cell",
 "Macrophage",
 "Macrophage",
 "Hepatocyte")
#################################################################################
## 
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
saveRDS(mydata,file = "01_scRNA/mydata.RDS")
#mydata=readRDS("01_scRNA/mydata.RDS")
u_plot=UMAPPlot(mydata, pt.size=0.2, label=T, label.size=4)+NoLegend()+scale_color_manual(values=colors)+theme(text=element_text(family="Times"))
ggsave("01_scRNA/UMAPPlot.pdf", u_plot, width=7, height=6)
genes = unique(c(
  "GZMB","GZMH","FGFBP2","GNLY","CD160",
  "GZMB","GZMH","FGFBP2","GNLY","CD160",
  "CD2","CD3D","CD3E","IL7R",
  "CD2","CD3D","CD3E","IL7R",
  "AIF1","MS4A7","LYZ",
  "DERL3","IGHG1","MZB1",
  "MKI67","MKI67","STMN1",
  "CD79A","MS4A1",
  "AIF1","MS4A7","LYZ",
  "AIF1","MS4A7","LYZ",
  "ALB","APOA1","KRT18"))

diff.marker.dotplot1= DotPlot(object = mydata, features = unique(genes),
                              dot.scale =6,
                              #dot.min = 0,
                              scale =T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()+scale_color_gradientn(colors=c( "grey","snow", "chartreuse4"))
#ggsave("01_scRNA/diff.marker.dotplot1.pdf", diff.marker.dotplot1, width=6, height=6)


####################################################################### 
cell.prop<-as.data.frame(prop.table(table(Idents(mydata), as.vector(mydata$Sample))))
colnames(cell.prop)<-c("Cell_type","Samples","proportion")
unique(cell.prop$Samples)

bili=ggplot(cell.prop,aes(Samples,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=colors)+
  ggtitle("")+
  theme_bw()+
  theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")

mydata@meta.data$Type=""
mydata@meta.data$Type[which(mydata@meta.data$Sample%in%c("HCC1","HCC2","HCC3"))]="HCC"
mydata@meta.data$Type[which(mydata@meta.data$Sample%in%c("HD1","HD2","HD3"))]="HD"



cell.prop2<-as.data.frame(prop.table(table(Idents(mydata), as.vector(mydata$Type))))
colnames(cell.prop2)<-c("Cell_type","Sample_type","proportion")
unique(cell.prop2$Sample_type)

bili2=ggplot(cell.prop2,aes(Sample_type,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=colors)+
  ggtitle("")+
  theme_bw()+
  theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12,  hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")


Figure1=mg_merge_plot(u_plot,diff.marker.dotplot1,bili2,bili,ncol=2,nrow=2,labels = c('A','B','C','D'),legend ="right")
ggsave('PDFs/Fig1.pdf',Figure1,height = 10,width = 12)
ggsave('PDFs/Fig1.jpg',Figure1,height = 10,width = 12)




######################


table(mydata$cell_type)
immuce_cells=subset(mydata,cell_type%in%c('Natural killer cell'))

immuce_cells <- NormalizeData(immuce_cells)
immuce_cells <- FindVariableFeatures(immuce_cells)
immuce_cells <- ScaleData(immuce_cells)
immuce_cells <- RunPCA(immuce_cells, verbose=FALSE)
immuce_cells = RunHarmony(immuce_cells, group.by.vars="Sample")

immuce_cells <- RunUMAP(immuce_cells, dims=1:20, reduction="harmony")
#immuce_cells <- RunTSNE(immuce_cells, dims=1:20, reduction="harmony")

DimPlot(immuce_cells, reduction="umap", group.by="Sample", pt.size=1)+theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
#
immuce_cells <- FindNeighbors(immuce_cells, dims=1:20, reduction="harmony")
immuce_cells <- FindClusters(immuce_cells, resolution=0.1)
markers <- FindAllMarkers(immuce_cells, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "02_NKcell/NKcell_All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
markers <- read.table("02_NKcell/NKcell_All_celltype_DEG.txt",header = T,check.names = F,fill=T,sep = "\t")
markers["pct.diff"]=markers$pct.1-markers$pct.2

#################################################################################
### 
Top51 <- markers %>% group_by(cluster) %>% slice_max(n =20, order_by =pct.diff )
Cellmarker <- read.table("01_scRNA/Cell_marker_Human.txt",header = T,check.names = F,fill=T,sep = "\t")
Cellmarker <- Cellmarker[,c(9,7)]
colnames(Top51)[7]="marker"
Cellmarker2 <- merge(Cellmarker,Top51,by="marker")
write.table(Cellmarker2, "02_NKcell/NKcell_Cellmarker2.txt", col.names=T, row.names=F, quote=F, sep="\t")
DimPlot(immuce_cells, reduction="umap", group.by="seurat_clusters", pt.size=1,label = T)+theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)

VlnPlot(immuce_cells, features=c("NKG7","GZMH" ), pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())

#0:FGFBP2+ NK cell  "FGFBP2" , "FCGR3A","GZMB", "GZMH"  
#1:CD160+ NK cell:  "GZMK" ,"CD160", 
#2:IL7R+ NK cell: "IL7R" ,  "SELL"  ,"LMNA"   ,"CD44" 
cell_label =c("FGFBP2+ NK cell","CD160+ NK cell","IL7R+ NK cell")

names(cell_label) <- levels(immuce_cells)
immuce_cells <- RenameIdents(immuce_cells, cell_label)
immuce_cells[["cell_label"]] = Idents(immuce_cells)
color2=c("#EF8080","#B8D2ED","#EDD2ED")

saveRDS(immuce_cells,file = '02_NKcell/immuce_cells.rds')

u_plot_Macrophage=UMAPPlot(immuce_cells, pt.size=1, label=T, label.size=4)+NoLegend()+scale_color_manual(values=color2)+theme(text=element_text(family="Times"))

genes=c("FGFBP2" , "FCGR3A","GZMB", "GZMH"  , "GZMK" ,  "CD160",  "IL7R" ,  "SELL"  ,"LMNA"   ,"CD44"  )
diff.marker.dotplot_subcell= DotPlot(object = immuce_cells, features = as.vector(genes),
                                        dot.scale =6,
                                        #dot.min = 0,
                                        scale =T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()+scale_color_gradientn(colors=c( "grey","snow", "red"))
genes2=c("FGFBP2" ,  "CD160",  "IL7R" )

Vln_plot=VlnPlot(immuce_cells, features=as.vector(genes2), pt.size=0, stack = T, cols=c(color2))+NoLegend()+theme(axis.title.x=element_blank())

marker_heatmap = DoHeatmap(immuce_cells, features=Top51$marker, group.colors=color2, label=FALSE)+scale_fill_gradientn(colors=c("white", "grey", "firebrick3"))+theme(text=element_text(family="Times"))+theme(legend.title = element_blank())



cell.prop3<-as.data.frame(prop.table(table(Idents(immuce_cells), as.vector(immuce_cells$Type))))
colnames(cell.prop3)<-c("Cell_type","Sample_type","proportion")
unique(cell.prop3$Sample_type)

bili3=ggplot(cell.prop3,aes(Cell_type,proportion,fill=Sample_type))+
  geom_bar(stat="identity",position = 'dodge')+
  scale_fill_manual(values=c("#d5231d","#3777ac"))+
  ggtitle("")+
  theme_bw()+
  theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12,  angle = 45,hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")




#################
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(survival)
library(survminer)
tcga.t.exp=readMatrix('00_origin_datas/TCGA/tcga.t.exp.txt')
tcga.n.exp=readMatrix('00_origin_datas/TCGA/tcga.n.exp.txt')
tcga.t.cli=readMatrix('00_origin_datas/TCGA/tcga.t.cli.txt')
sample_del=setdiff(colnames(tcga.t.exp), tcga.t.cli$SampleID)
tcga.t.exp=tcga.t.exp[,-which(colnames(tcga.t.exp)==sample_del)]
tcga.t.cli=tcga.t.cli[match(colnames(tcga.t.exp),tcga.t.cli$SampleID),]
identical(tcga.t.cli$SampleID,colnames(tcga.t.exp))
#0:FGFBP2+ NK cell  "FGFBP2" , "FCGR3A","GZMB", "GZMH"  
#1:CD160+ NK cell:  "GZMK" ,"CD160", 
#2:IL7R+ NK cell: "IL7R" ,  "SELL"  ,"LMNA"   ,"CD44" 

gene_set = read.table("02_NKcell/NKcell_All_celltype_DEG.txt",header = T,check.names = F,fill=T,sep = "\t")
gene_set = subset(gene_set, cluster=="0")
list <- list()
list[['0']] <- as.vector(unique(gene_set$gene))

ssgsea_t_score <- ssGSEAScore_by_muti_group_genes(gene.exp = tcga.t.exp,
                                                  genelist = list)
ssGSEA_score <- as.data.frame(t(ssgsea_t_score))
ssGSEA_score$"SampleID"=rownames(ssGSEA_score)
survival_df = tcga.t.cli[, c("SampleID","OS.time", "OS")]
score_survival = merge(x=survival_df, y=ssGSEA_score, by="SampleID")
colnames(score_survival)[4]="subcell NK cell"
tcga.cutoff <- survminer::surv_cutpoint(data.frame(time=score_survival$OS.time/365, event=score_survival$OS, Score=score_survival$`subcell NK cell`),time = "time", event = "event",variables = c("Score"))
tcga.cutoff.cut=surv_categorize(tcga.cutoff)

#write.table(score_survival, "./score_survival.txt", col.names=T, row.names=T, sep="\t", quote=F)

fit = survfit(Surv(time, event)~Score, data=tcga.cutoff.cut)
sur_plot=ggsurvplot(fit,
           data = tcga.cutoff.cut,
           conf.int = TRUE,
           pval = TRUE,
           surv.median.line = "hv",
           risk.table = TRUE,
           xlab = "Time(years)",
           palette = c("tomato1", "steelblue1"),
           #legend.labs = c("Higher Enrichment", "Lower Enrichment"),
           legend.title = "Risk",
           break.x.by = 3,
           pval.size = 7,
           font.legend = 12
)
fig2a=mg_merge_plot(u_plot_Macrophage,marker_heatmap,ncol=2,nrow=1,widths = c(1,1.2),labels = LETTERS[c(1,2)],legend = "bottom")
fig2b=mg_merge_plot(diff.marker.dotplot_subcell,bili3,sur_plot$plot,ncol=3,nrow=1,labels =  LETTERS[c(3,4,5)],widths = c(1,1,1.2))
fig2=mg_merge_plot(fig2a,fig2b,ncol=1,nrow=2)
#CD160=sur_plot$plot
#IL7R=sur_plot$plot
#figS2=mg_merge_plot(CD160,IL7R,ncol=2,nrow=1)
#ggsave('PDFs/Figure S2.pdf',figS2,height = 5,width = 8)
#ggsave('PDFs/Figure S2.jpg',figS2,height = 5,width = 8)

ggsave('PDFs/Fig2.pdf',fig2,height = 10,width = 12)
ggsave('PDFs/Fig2.jpg',fig2,height = 10,width = 12)


#####################

immuce_cells <- readRDS('02_NKcell/immuce_cells.rds')
monocle_data=subset(immuce_cells,cell_label%in%c('FGFBP2+ NK cell'))


##---------------1.⚠️-----------------------------
##
expr_matrix <- as(as.matrix(monocle_data@assays$RNA@counts), 'sparseMatrix')
##p_data(phenotype_data)
p_data <- monocle_data@meta.data
# p_data$celltype <- monocle_data@active.ident  ##
##
f_data <- data.frame(gene_short_name = row.names(monocle_data),row.names = row.names(monocle_data))
##expr_

##-----------------------
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
#
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
##----------------------------------3.----------------------------------
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)



##-----------------------------------4. -----------------------------------
cds <- detectGenes(cds, min_expr = 0.1) #num_cells_expressed
print(head(fData(cds)))#
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) #
#---------------------------------5. 
#######
diff = differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr="~Type", cores=16)
deg = subset(diff, qval<1e-08)
deg = deg[order(deg$qval, decreasing=F), ]
#######
#~
#######
#######Idents(monocle_data) = monocle_data@meta.data$cell_label
#######deg = FindAllMarkers(monocle_data, logfc.threshold=0.25, min.pct=0.25, only.pos=T)
#######
##

## 
ordergene <- unique(deg$gene_short_name )
#######ordergene <- unique(deg$gene)

# ordergene <-diff
cds <- setOrderingFilter(cds, ordergene)  
#
#[["use_for_ordering"]]，table(cds@featureData@data[["use_for_ordering"]])
pdf("03_Pseudotime/train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
#
#⚠️
#gene 
#ordergene <- row.names(deg)[order(deg$qval)][1:400]

##########Step 2: 
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

######Step 3: 
###1.pseudotime 
cds1 <- orderCells(cds)
####?????
cds2 <- orderCells(cds1, root_state=8)

#⚠️
colors1 <- c("#4ea64a","#8e4c99","#d5231d","#3777ac","#a05528","#1f2d6f","#279772","black","#d9b71a","#EF8080","#B8D2ED")
p1=plot_cell_trajectory(cds2,color_by="Pseudotime", cell_size=0.5,show_backbone=TRUE) 
p4=plot_cell_trajectory(cds2,color_by="State", cell_size=0.5,show_backbone=TRUE)+scale_color_manual(values=colors1) 
p2=plot_cell_trajectory(cds2,color_by="Type", cell_size=0.5,show_backbone=TRUE)+scale_color_manual(values=c("#d5231d","#3777ac"))
#+scale_color_manual(breaks = c('CA','NL'), values=c("#DE582B","#018A67"))
p3=plot_cell_trajectory(cds2, color_by = "cell_label",cell_size=0.5,show_backbone=TRUE)+scale_color_manual(values=color2)#+theme(legend.position = "right")


p2|p1
Time_diff <- differentialGeneTest(cds2[ordergene,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #
write.csv(Time_diff, "03_Pseudotime/Time_diff_all.csv", row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p_all=plot_pseudotime_heatmap(cds2[Time_genes,],num_clusters = 2, show_rownames=F, return_heatmap=T)
dev.off()
ggsave("03_Pseudotime/Time_heatmapAll1.pdf", p_all, width =4, height =5)

clusters <- cutree(p_all$tree_row, k = 2)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.csv(clustering, "03_Pseudotime/Time_clustering_all.csv", row.names = T)
#
# Time_genes <- top_n(Time_diff, n = 50, desc(qval)) %>% pull(gene_short_name) %>% as.character()
Time_diff_res=cbind.data.frame(Clusters=clustering[Time_diff$gene_short_name,],Time_diff)
head(Time_diff_res)
write.csv(Time_diff_res, "03_Pseudotime/Time_diff_res.csv", row.names = T)
# 

Top20 <- Time_diff_res %>% group_by(Clusters) %>% slice_min(n=100,order_by = qval)  
p = plot_pseudotime_heatmap(cds2[Top20$gene_short_name,], num_clusters=2, show_rownames=F, return_heatmap=T)
dev.off()
p_plot=grid::grid.draw(p$gtable)
ggsave("03_Pseudotime/Time_heatmapTop100.pdf", p, width = 4, height = 10)

#
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(Time_diff_sig, "03_Pseudotime/Time_diff_sig.csv", row.names = F)


####
#BEAM_res2 = BEAM(cds[ordergene,], branch_point=2, cores=1)
#BEAM_res1 = BEAM_res[order(BEAM_res2$qval), ]
#BEAM_res1 = BEAM_res1[, c("gene_short_name", "pval", "qval")]
#BEAM_res1 = subset(BEAM_res1, qval<0.01)
#p_BEAM = plot_genes_branched_heatmap(cds[row.names(BEAM_res1),], branch_point=2, num_clusters=2, show_rownames=F, branch_labels=c("Cell fate 1", "Cell fate 2"), return_heatmap=T, branch_colors=c("peachpuff1", "chocolate1", "olivedrab1"))
#branch_clusters = cutree(p_BEAM$ph_res$tree_row, k=3)
#branch_clustering = data.frame(branch_clusters)
#branch_clustering[, 1] = as.character(branch_clustering[, 1])
#colnames(branch_clustering) = "Branch_Clusters"
#write.table(branch_clustering, "03_Pseudotime/BEAM_clusters.txt", col.names=T, row.names=T, sep="\t", quote=F)
####
cluster1.enrich=mg_clusterProfiler(genes =Time_diff_res$gene_short_name[Time_diff_res$Clusters==1])
Clusters1_kegg=enrichplot::dotplot(cluster1.enrich$KEGG)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
ggsave("03_Pseudotime/Clusters1_kegg.pdf", Clusters1_kegg, width =5, height =5)

Fig3=mg_merge_plot(p2,p1,p$gtable,Clusters1_kegg,ncol=2,nrow=2,labels = c('A','B','C','D'))
ggsave('PDFs/Fig3.pdf',Fig3,height = 10,width = 10)
ggsave('PDFs/Fig3.jpg',Fig3,height = 10,width = 10)



############CellChat#########
##1.4 ########
mydata=readRDS("01_scRNA/mydata.RDS")
immuce_cells <- readRDS('01_scRNA/immuce_cells.rds')
sc_df2=data.frame(immuce_cells@meta.data, immuce_cells@reductions$umap@cell.embeddings)

table(mydata@meta.data$cell_type)
table(immuce_cells@meta.data$cell_label)
dim(immuce_cells@meta.data)


mydata$cell_type2=''
mydata$cell_type2=as.vector(mydata$cell_type)
mydata$cell_type2[rownames(sc_df2)]=as.vector(sc_df2$cell_label)
table(mydata$cell_type2)
comm_cells=mydata

table(comm_cells$cell_type2)
colnames(comm_cells@meta.data)

# #
DefaultAssay(comm_cells) <- 'RNA'
cellchat <- createCellChat(object = comm_cells, meta = comm_cells@meta.data, group.by = "cell_type2")
cellchat@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat@data.signaling
#
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
# # #
# #
# # #
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)
# # 
library(NMF)
# selectK(cellchat, pattern = "outgoing")3
# selectK(cellchat, pattern = "incoming")4

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
# save(cellchat,file='results/01.scRNA/CellChat/cellchat.RData')
# load('results/01.scRNA/CellChat/cellchat.RData')

groupSize <- as.numeric(table(cellchat@idents))

cellchat@netP$pathways
length(cellchat@netP$pathways)
#8


netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing") 
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming") 
dev.off()
pdf('PDFs/cellchat_circle.pdf',height = 5,width = 8)
par(mfrow=c(1,3))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F,
                 title.name = "Interaction weights/strength")

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F
                 ,title.name = "FGFBP2+ NK cell",sources.use = c('FGFBP2+ NK cell'))


dev.off()


table(comm_cells$cell_type2)
pdf('PDFs/cellchat_dotplot_1.pdf',height = 5,width = 5)
netVisual_bubble(cellchat, signaling = cellchat@netP$pathways,sources.use = c('FGFBP2+ NK cell'))
dev.off()

pdf('PDFs/cellchat_dotplot_2.pdf',height = 5,width = 5)
netVisual_bubble(cellchat, signaling = cellchat@netP$pathways,targets.use = c('FGFBP2+ NK cell'), remove.isolate = T)
dev.off()




