rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(devtools))
suppressMessages(library(customLayout))
suppressMessages(library(stringr))
suppressMessages(library(ConsensusClusterPlus))
suppressMessages(library(tidydr))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(fgsea))
suppressMessages(library(corrplot))
suppressMessages(library(colorspace))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(maftools))
suppressMessages(library(vegan))
suppressMessages(library(forcats))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(rstatix))
suppressMessages(library(ggstatsplot))
suppressMessages(library(ggcor))
suppressMessages(library(ggstance))
suppressMessages(library(tidyverse))
suppressMessages(library(GOplot))
suppressMessages(library(caret))
suppressMessages(library(writexl))
suppressMessages(library(rcartocolor))
suppressMessages(library(ggcorrplot))
suppressMessages(library(psych))
suppressMessages(library(clusterProfiler))
suppressMessages(library(dplyr))
suppressMessages(library(cols4all))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(scales))
suppressMessages(library(oncoPredict))
suppressMessages(library(gghalves))
suppressMessages(library(cowplot))
suppressMessages(library(IOBR))
suppressMessages(library(estimate))
suppressMessages(library(UpSetR))
suppressMessages(library(ggbiplot))
suppressMessages(library(ggsci))
suppressMessages(library(WGCNA))
suppressMessages(library(circlize))
suppressMessages(library(rJava))
suppressMessages(library(xlsxjars))
suppressMessages(library(xlsx))
suppressMessages(library(glmnet))
suppressMessages(library(tidyr))
suppressMessages(library(pROC))
suppressMessages(library(ROCR))
#01Cluster###########################
tcga.t.exp=readMatrix('00_origin_datas/TCGA/tcga.t.exp.txt')
tcga.n.exp=readMatrix('00_origin_datas/TCGA/tcga.n.exp.txt')
tcga.t.cli=readMatrix('00_origin_datas/TCGA/tcga.t.cli.txt')
sample_del=setdiff(colnames(tcga.t.exp), tcga.t.cli$SampleID)
tcga.t.exp=tcga.t.exp[,-which(colnames(tcga.t.exp)==sample_del)]
tcga.t.cli=tcga.t.cli[match(colnames(tcga.t.exp),tcga.t.cli$SampleID),]
identical(colnames(tcga.t.exp),as.vector(tcga.t.cli$SampleID))

##############lasso##############
#########03######
gene_set = read.table("02_NKcell/NKcell_All_celltype_DEG.txt",header = T,check.names = F,fill=T,sep = "\t")
gene_set = subset(gene_set, cluster=="0")
tcga.deg.sig=as.vector(unique(gene_set$gene))

length(tcga.deg.sig)



tcga.t.exp_use=tcga.t.exp
tcga.subtype.cli=tcga.t.cli
rownames(tcga.subtype.cli)=tcga.subtype.cli$SampleID
colnames(tcga.subtype.cli)[1]="Samples"

tcga.cox=cox_batch(t(scale(t(tcga.t.exp_use[tcga.deg.sig,])))
                   ,time =  tcga.subtype.cli$OS.time/365
                   ,event =tcga.subtype.cli$OS)
dim(tcga.cox)

table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)
writeMatrix(tcga.cox,outpath = '04_Lasso/tcga.cox.txt')



p.cutoff=0.05
tcga.cox_use=tcga.cox
tcga.cox_use$coef=log(tcga.cox_use$HR)
tcga.cox_use$Gene=rownames(tcga.cox_use)
tcga.cox_use$type=rep('None',nrow(tcga.cox_use))
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef>0)]='Risk'
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef<0)]='Protective'
table(tcga.cox_use$type)

######### lasso
tcga.gene.sig=rownames(tcga.cox)[which(tcga.cox$p.value<p.cutoff)]
length(tcga.gene.sig)

table(tcga.cox_use$type)


tcga.cox_forVis=tcga.cox_use
tcga.cox_forVis=tcga.cox_forVis[which(tcga.cox_forVis$type %in% c('Risk','Protective')),]
tcga.cox_forVis$p.value=-log10(tcga.cox_forVis$p.value)
range(tcga.cox_forVis$p.value)


#################### LASSO
table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)

tcga.exp.sig=tcga.t.exp_use[tcga.gene.sig,]
tcga.exp.sig=t(tcga.exp.sig)
dim(tcga.exp.sig)


dim(tcga.exp.sig)
options(ggrepel.max.hnscerlaps = Inf)
tcga.subtype.cli$Samples=as.vector(tcga.subtype.cli$Samples)
identical(rownames(tcga.exp.sig),tcga.subtype.cli$Samples)
tcga.lasso.res=mg_lasso_cox_use(tcga.exp.sig
                                , time = tcga.subtype.cli$OS.time/365
                                , event = tcga.subtype.cli$OS
                                , nfolds = 10
                                , lambda.min = T
                                , figLabels=c('B','C'))
tcga.lasso.res$Genes

tcga.lasso.res$lambda

tcga.lasso.res$plot

tcga.exp.for.cox=tcga.t.exp_use[match(tcga.lasso.res$Genes,row.names(tcga.t.exp_use)),]
dim(tcga.exp.for.cox)
identical(colnames(tcga.exp.for.cox),tcga.subtype.cli$Samples)

lst.modl=createCoxModel_use((t(tcga.exp.for.cox))
                            , time = tcga.subtype.cli$OS.time/365
                            , event = tcga.subtype.cli$OS
                            , isStep = T)
lst.modl$Cox
lst.modl$Genes
lst.modl$fmla

lst.modl.Coef=lst.modl$Coef
names(lst.modl.Coef)=lst.modl$Genes
lst.modl.Coef

tcga.risk.score=lst.modl$Score
#tcga.risk.score=scale(tcga.risk.score)[,1]
tcga.risk.score=mosaic::zscore(tcga.risk.score)

range(tcga.risk.score)

lst.modl$Coef

gene.coef=data.frame(Gene=lst.modl$Genes,Coef=lst.modl$Coef)
gene.coef$Type=ifelse(lst.modl$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
table(gene.coef$Type)

fig3c=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  geom_text(aes(label=round(Coef,digits = 3)),color="black",hjust = "left")+
  
  coord_flip() +
  scale_fill_manual(values=pal_nejm(alpha = 0.9)(8)[c(1,2)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "Lasso Cox coefficient") +
  theme_classic()+theme(legend.position = c(0,1))
# theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position="top")


fig3AB=mg_plot_lasso_use(fit = tcga.lasso.res$Mode1
                         , cv_fit = tcga.lasso.res$Model2
                         , show_text = F
                         , figLabels = c('A', 'B'))
fig3AB
#fig3abc=mg_merge_plot(fig1a,fig3AB,fig3c,nrow = 1,ncol = 3,widths = c(1,2,1))
#savePDF('PDFs/Fig7AB.pdf',fig7A,height = 4,width = 9)
#savePDF('PDFs/fig3abc.pdf',fig3abc,height = 5,width = 15)


tcga.exp.forCox<- cbind(time=tcga.subtype.cli$OS.time/365,
                        status=tcga.subtype.cli$OS,
                        t(tcga.t.exp_use)[rownames(tcga.subtype.cli), lst.modl$Genes])



dim(tcga.exp.forCox)

fmla <- as.formula(paste0("Surv(time, status) ~",paste0(lst.modl$Genes,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga.exp.forCox))
fig3d=survminer::ggforest(cox,data=tcga.exp.forCox,noDigits = 3)
fig3abc=mg_merge_plot(fig3AB,fig3c,nrow = 1,ncol =2,widths = c(2,1),labels = c("","C"))

#savePDF('PDFs/fig3abc.pdf',fig3abcd,height = 5,width = 16)


############### TCGA
#tcga.cutoff <- survminer::surv_cutpoint(data.frame(time=tcga.subtype.cli$OS.time/365, event=tcga.subtype.cli$OS, risk=tcga.risk.score),time = "time", event = "event",variables = c("risk"))
#tcga.cutoff=tcga.cutoff$cutpoint$cutpoint
#tcga.cutoff=median(tcga.risk.score)
tcga.cutoff=0
identical(colnames(tcga.exp.for.cox),tcga.subtype.cli$Samples)
risk.group.color=c("#EE0000FF","#3B4992FF")
names(risk.group.color)=c('High','Low')
tcga.roc=plotCoxModel_Batch_use(riskScore = tcga.risk.score
                                ,dat = t(tcga.exp.for.cox[match(lst.modl$Genes, row.names(tcga.exp.for.cox)), ])
                                , time = tcga.subtype.cli$OS.time/365
                                , event = tcga.subtype.cli$OS
                                , cutoff = tcga.cutoff
                                , labs = c('High','Low')
                                , title = 'RiskType'
                                , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                                
                                , pal = risk.group.color
                                , mks = c(1,3,5))
tcga.roc1=tcga.roc[[1]]
#pdf('PDFs/fig3e2.pdf',height = 6,width = 6)
tcga.roc1

dev.off()

tcga.group=ifelse(tcga.risk.score>tcga.cutoff,'High','Low')
tcga.group=data.frame(tcga.group)
colnames(tcga.group)='group'
table(tcga.group$group)
write.table(cbind(tcga.risk.score,tcga.group),file = '04_Lasso/tcga.group.txt',sep='\t',quote = F)


###########################

icgc.t.exp=readMatrix('00_origin_datas/ICGC/icgc.t.exp.txt')
icgc.t.cli=readMatrix('00_origin_datas/ICGC/icgc.t.cli.txt')
identical(colnames(icgc.t.exp),icgc.t.cli$SampleID)
icgc.HCCDB18.t.exp=icgc.t.exp
match(lst.modl$Genes,row.names(icgc.HCCDB18.t.exp))
length(lst.modl$Genes)
icgc.HCCDB18.t.cli.os=icgc.t.cli
identical(icgc.HCCDB18.t.cli.os$SampleID , colnames(icgc.HCCDB18.t.exp)) ####
icgc.HCCDB18.model.dat=icgc.HCCDB18.t.exp[match(lst.modl$Genes,row.names(icgc.HCCDB18.t.exp)),]

icgc.HCCDB18.risk.score=predictScoreByCoxModel(coxModel = lst.modl
                                               ,(t(icgc.HCCDB18.model.dat)))
icgc.HCCDB18.risk.score=mosaic::zscore(icgc.HCCDB18.risk.score)

lst.modl$fmla
#icgc.HCCDB18.cutoff <- survminer::surv_cutpoint(data.frame(time=icgc.HCCDB18.t.cli.os$OS.time/365, event=icgc.HCCDB18.t.cli.os$OS, risk=icgc.HCCDB18.risk.score),time = "time", event = "event",variables = c("risk"))
#icgc.HCCDB18.cutoff=icgc.HCCDB18.cutoff$cutpoint$cutpoint

icgc.HCCDB18.cutoff=0

icgc.ROC=plotCoxModel_Batch_use(riskScore = icgc.HCCDB18.risk.score
                                , dat = t(icgc.HCCDB18.t.exp[intersect(lst.modl$Genes, row.names(icgc.HCCDB18.t.exp)),])
                                , time = as.numeric(icgc.HCCDB18.t.cli.os$OS.time/365) 
                                , event = as.numeric(icgc.HCCDB18.t.cli.os$OS)
                                , cutoff = icgc.HCCDB18.cutoff
                                , labs = c('High','Low')
                                , title = 'RiskType'
                                , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                                , pal = risk.group.color
                                , mks = c(1,2,3))
icgc.ROC1=icgc.ROC[[1]]

icgc.ROC1
dev.off()
icgc.HCCDB18.group=ifelse(icgc.HCCDB18.risk.score>icgc.HCCDB18.cutoff,'High','Low')
icgc.HCCDB18.group=data.frame(icgc.HCCDB18.group)
colnames(icgc.HCCDB18.group)='group'
table(icgc.HCCDB18.group)

Fig4_ROC=mg_merge_plot(tcga.roc1,icgc.ROC1,ncol=2,nrow=1)
Fig4=ggpubr::ggarrange(fig3abc,Fig4_ROC, ncol = 1, nrow = 2,heights = c(1,2))

ggsave('PDFs/Fig4.pdf',Fig4,height = 12,width = 15)
ggsave('PDFs/Fig4.jpg',Fig4,height = 12,width = 15)
identical(rownames(as.data.frame(tcga.risk.score)),rownames(tcga.subtype.cli))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=tcga.risk.score)
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff,'High','Low')
write.table(tcga.risktype.cli,file = '04_Lasso/tcga.risktype.cli.txt',sep='\t',quote = F)


#  #############
tcga.risktype.cli=read.table('04_Lasso/tcga.risktype.cli.txt',header = T,check.names = F,fill=T,sep = "\t")
tcga.t.exp=readMatrix('00_origin_datas/TCGA/tcga.t.exp.txt')
tcga.t.cli=readMatrix('00_origin_datas/TCGA/tcga.t.cli.txt')
sample_del=setdiff(colnames(tcga.t.exp), tcga.t.cli$SampleID)
tcga.t.exp=tcga.t.exp[,-which(colnames(tcga.t.exp)==sample_del)]
identical(as.vector(tcga.risktype.cli$Samples),colnames(tcga.t.exp))

risk.group.color1=c("#EE0000FF","#3B4992FF")
names(risk.group.color1)=c('High','Low')

library(estimate)
#### ESTIMATE
#tcga.exp.estimate<-deconvo_estimate(eset=tcga.t.exp)
#save(tcga.exp.estimate,file='05_imm/tcga.exp.estimate.RData')
load('05_imm/tcga.exp.estimate.RData')
tcga.exp.estimate=get.IOBR.immu.format(tcga.exp.estimate)


############ TIMER 
#tcga.exp.timer<-deconvo_timer(eset=as.matrix(tcga.t.exp),indications=rep('LIHC',ncol(tcga.t.exp)))
#save(tcga.exp.timer,file='05_imm/tcga.exp.timer.RData')
load('05_imm/tcga.exp.timer.RData')
tcga.exp.timer=get.IOBR.immu.format(tcga.exp.timer)

library(EPIC)
############ EPIC 
#tcga.exp.epic<-deconvo_epic(eset=as.matrix(tcga.t.exp),tumor = TRUE)
#save(tcga.exp.epic,file='05_imm/tcga.exp.epic.RData')
load('05_imm/tcga.exp.epic.RData')
tcga.exp.epic=get.IOBR.immu.format(tcga.exp.epic)

############ MCP-counter 
#tcga.exp.mcp<-deconvo_mcpcounter(eset=as.matrix(tcga.t.exp))
#save(tcga.exp.mcp,file='05_imm/tcga.exp.mcp.RData')
load('05_imm/tcga.exp.mcp.RData')
tcga.exp.mcp=get.IOBR.immu.format(tcga.exp.mcp)

### CIBERSORT
#tcga.exp.cibersort<-deconvo_cibersort(eset=tcga.t.exp,arrays=T)
#save(tcga.exp.cibersort,file='05_imm/tcga.exp.cibersort.RData')
load('05_imm/tcga.exp.cibersort.RData')
tcga.exp.cibersort=get.IOBR.immu.format(tcga.exp.cibersort)

#######sssGSEA#######
#geo.immu.ssgsea=immu_ssgsea(exp = tcga.t.exp)
#save(geo.immu.ssgsea,file='05_imm/geo.immu.ssgsea.RData')
load('05_imm/geo.immu.ssgsea.RData')



#
tcga.t.estimate=tcga.exp.estimate[rownames(tcga.risktype.cli),1:3]
tcga.t.timer=tcga.exp.timer[rownames(tcga.risktype.cli),]
tcga.t.epic=tcga.exp.epic[rownames(tcga.risktype.cli),]
tcga.t.mcp=tcga.exp.mcp[rownames(tcga.risktype.cli),]
tcga.t.cibersort=tcga.exp.cibersort[rownames(tcga.risktype.cli),1:22]
tcga.t.ssGSEA28=as.data.frame(geo.immu.ssgsea[rownames(tcga.risktype.cli),])

fig5a=get_PlotMutiBoxplot(tcga.t.estimate,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')

#fig5a=groupViolin(tcga.t.estimate,
#            tcga.risktype.cli$Risktype,
#            ylab = 'Score',
#            group_col=risk.group.color1)

fig5a

fig5b=get_PlotMutiBoxplot(tcga.t.timer,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')

#fig5b=groupViolin(tcga.t.timer,
#                  tcga.risktype.cli$Risktype,
#                  ylab = 'Score',
#                  group_col=risk.group.color1)


fig5b

fig5c=get_PlotMutiBoxplot(tcga.t.epic,tcga.risktype.cli
                          ,group_cols = risk.group.color
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5c

fig5d=get_PlotMutiBoxplot(tcga.exp.mcp,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          #,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')

#fig5d=groupViolin(tcga.exp.mcp,
#            tcga.risktype.cli$Risktype,
#            ylab = 'Score',
#            group_col=risk.group.color1)
fig5d





fig5e=get_PlotMutiBoxplot(tcga.t.cibersort,tcga.risktype.cli
                          ,group_cols = risk.group.color
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5e


fig5f=get_PlotMutiBoxplot(tcga.t.ssGSEA28,tcga.risktype.cli                          ,group_cols = risk.group.color1,ylab = 'Score' ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')

fig5f

fig6a=mg_merge_plot(fig5a,fig5b,fig5d,nrow = 1,ncol = 3
                   ,labels = LETTERS[1:3])




# ICGs####

tcga.icgs=immu_ICGs(tcga.t.exp_use)
colnames(tcga.icgs)
#select_col=c('CTLA4','PDCD1',"TIGIT",'LGALS9','CD276','TNFSF4','TNFSF9')
select_col=c("CTLA4","LGALS9",'PDCD1',"PDCD1LG2","CD28","CD80","CD244","CD276" , "CD44" )
#select_col=colnames(tcga.icgs)
icg.dat.RS=cbind(tcga.risktype.cli$Riskscore
                 ,t(tcga.t.exp_use)[rownames(tcga.risktype.cli),lst.modl$Genes]
                 ,tcga.icgs[rownames(tcga.risktype.cli),c(select_col)])
#c('CTLA4','PDCD1','PDCD1LG2',"TIGIT",'LGALS9','CD80','CD28','HAVCR2')
colnames(icg.dat.RS)[1]='Riskcsore'

icg_cor_res <- Hmisc::rcorr(as.matrix(icg.dat.RS),type = 'spearman')
icg_cor_res$P[is.na(icg_cor_res$P)] <- 0
icg_cor_res.p=icg_cor_res$P
icg_cor_res.p[1:5,1:5]
icg_cor_res.p<-ifelse(icg_cor_res.p<0.0001,'****',
                      ifelse(icg_cor_res.p<0.001,'***', 
                             ifelse(icg_cor_res.p<0.01,'**',
                                    ifelse(icg_cor_res.p<0.05,'*',''))))
ICG_plot <- pheatmap(icg_cor_res$r[-c(1:(length(lst.modl$Genes)+1)),c(lst.modl$Genes,'Riskcsore')],
                     color = circlize::colorRamp2(c(-1,-0.5, 0,0.5, 1), c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988","#BB4444" )),
                     main="Heatmap", # 
                     display_numbers = icg_cor_res.p[-c(1:(length(lst.modl$Genes)+1)),c(lst.modl$Genes,'Riskcsore')], # 
                     annotation_legend = FALSE ,
                     cluster_cols = F, # 
                     cluster_rows = F,
                     show_rownames = T, #
                     show_colnames = T,
                     fontsize_row = 12, # 
                     fontsize_col = 16)


pdf('PDFs/ICG_plot.pdf',height = 6,width = 6)
ICG_plot
dev.off()

#####TIDE######
tcga.tide<-read.csv('05_imm/TIDE.csv',row.names = 1,stringsAsFactors = F)
tcga.tide=tcga.tide[rownames(tcga.risktype.cli),]
dim(tcga.tide)

tide.selected=c('Exclusion','Dysfunction','TIDE')

tcga.subtype.tide.p.all=c()
for(i in tide.selected){
  df_half=data.frame("RiskType"=tcga.risktype.cli$Risktype,"TIDE_score"=tcga.tide[rownames(tcga.risktype.cli),i])
  p <-  ggplot(data = df_half, aes(x = RiskType, y = TIDE_score, fill = RiskType)) +
    geom_half_violin(side='R',position = position_nudge(x = 0.2, y = 0), alpha = 1) +
    geom_point(aes(y = TIDE_score, color = RiskType), position = position_jitter(width = 0.15), size = 1, alpha = 1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 1) +
    labs(y = i
         , x = 'RiskType') +
    guides(fill = FALSE, color = FALSE) +
    scale_fill_manual(values = risk.group.color) +
    scale_colour_manual(values = risk.group.color) +theme_classic2()+
    ggpubr::stat_compare_means(comparisons = list(c(1,2)),method = 'wilcox.test',label= "p.signif")
  
  tcga.subtype.tide.p.all=c(tcga.subtype.tide.p.all,list(p))
}
length(tcga.subtype.tide.p.all)

TIDE_plot=mg_merge_plot(tcga.subtype.tide.p.all,nrow = 1,ncol = length(tide.selected),common.legend=T,legend =  "left")
TIDE_plot


################imm############


tcga_durg_ic50_res=read.csv('05_imm/calcPhenotype_Output/DrugPredictions.csv',row.names = 1)
tcga.risktype.cli=read.table('04_Lasso/tcga.risktype.cli.txt',header = T,check.names = F,fill=T,sep = "\t")

tcga.risktype.cli=tcga.risktype.cli[rownames(tcga_durg_ic50_res),]
dim(tcga_durg_ic50_res)


IC50.mat=data.frame(Riskscore=tcga.risktype.cli$Riskscore,tcga_durg_ic50_res[as.vector(tcga.risktype.cli$Samples),])

IC50_RS_cor <- corr.test(x =IC50.mat$Riskscore,
                         y = IC50.mat[,-1])


IC50_RS_cor_res=data.frame(drugs=colnames( IC50.mat[,-1]))
IC50_RS_cor_res$cor<-as.numeric(IC50_RS_cor$r)
IC50_RS_cor_res$p.adj<-as.numeric(IC50_RS_cor$p.adj)
head(IC50_RS_cor_res)
table(IC50_RS_cor_res$p.adj<0.05,abs(IC50_RS_cor_res$cor)>0.3)
IC50_RS_cor_res=IC50_RS_cor_res[IC50_RS_cor_res$p.adj<0.05 & abs(IC50_RS_cor_res$cor)>0.3,]
IC50_RS_cor_res=IC50_RS_cor_res[order(IC50_RS_cor_res$cor),]
head(IC50_RS_cor_res)
IC50_plot <- ggplot(data=IC50_RS_cor_res,aes(x=cor,y=reorder(drugs,cor),
                                             color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_colour_gradient(low ='#ffc7c7' ,high = "#8785a2")+
  geom_segment(aes(yend=drugs,xend=0),size=1) +
  labs(x='spearman Correlation',y='')+theme_bw()+
  theme(text = element_text(family = 'Times'),legend.position = "bottom")


fig6b=mg_merge_plot(ICG_plot,TIDE_plot,IC50_plot,nrow = 1,ncol = 3,labels = LETTERS[4:6])
fig6=mg_merge_plot(fig6a,fig6b,nrow = 2,ncol =1,heights = c(1.5,1))

ggsave("PDFs/Fig6.pdf", fig6, width = 14, height = 10)
save.image("project.Rdata")
