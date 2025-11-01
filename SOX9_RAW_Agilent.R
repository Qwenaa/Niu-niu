
rm(list = ls())
##### Raw Data before analysis ----
## Targets
library(limma)
library(tidyverse)
library(GEOquery)
sample_name <- c('con_1.txt',"con_2.txt",'con_3.txt','shsox9_1.txt',"shsox9_2.txt","shsox9_3.txt")
group = rep(c("Con", "Test"), c(3,3))
Target <- data.frame(sample_name = sample_name,
                     group = group,
                     row.names = sample_name)


## read.maimages
shsox9_raw <- read.maimages(files = sample_name,
                               source = "agilent",#source代表是经过哪种程序得到的，有些Agilent芯片是通过genepix处理的
                               path = "./raw_data/",
                               names =sample_name,
                               other.columns = "gIsWellAboveBG",#读取进去是为了判断是否高于背景值
                               green.only = T)#代表是个单色芯片，默认是false

## add targets info
shsox9_raw$targets <- Target

## probe type summary
table(shsox9_raw$genes$ControlType)

## view data
head(shsox9_raw$E)[, 1:6]
dim(shsox9_raw)

## save data
save(shsox9_raw, file = "raw_data/shsox9_raw.Rdata")

##### 1 data normalization ----
rm(list = ls())
library(limma)
library(GEOquery)
library(tidyverse) 

###  1.1 load data ----
load_input <- load("raw_data/shsox9_raw.Rdata")
load_input

dim(shsox9_raw)

###  1.2 targets ----
shsox9_raw_targets <-shsox9_raw$targets

expr <-shsox9_raw$E
class(shsox9_raw$E)
#### 1.3 Plots before QC ----
##箱线图
library(RColorBrewer)#R配色的一个包
colors = rep(brewer.pal(6, "Set2")[1:2], each = 3)
boxplot(log2(shsox9_raw$E),
        ylab = expression(log[2](intensity)),
        las = 2,
        col = colors,
        outline = FALSE)

source('./Custom_Functions.R')
##PCA图
PCA_new(log2(shsox9_raw$E), 
        ntop = nrow(shsox9_raw$E),
        group = shsox9_raw$group,
        show_name = T)

#### 1.4 data normalization ----
shsox9_raw_bgc <-  backgroundCorrect(RG = shsox9_raw, 
                                    method = "normexp",#推荐normexp进行背景矫正
                                    offset = 50,#补偿值50
                                    normexp.method = "mle")#limma包推荐的mle
##数据标准化处理 芯片之间的水平,以进行log2处理
shsox9_norm <- normalizeBetweenArrays(shsox9_raw_bgc, #normalizeBetweenArrays适用于单色芯片
                                         method = "quantile")#normalizeinArrays适用于双色芯片
##提取表达矩阵
shsox9_expr_norm<- as.data.frame(shsox9_norm$E)
colnames(shsox9_expr_norm) <- shsox9_raw_targets$geo_accession

##查看基因表达数据中是否有缺失值
sum(is.na(shsox9_expr_norm))

#### 1.5 Plots after normalization ----
##  boxplot after normalization
boxplot(shsox9_expr_norm,
        ylab = expression(log[2](intensity)),
        las = 2,
        col = colors,
        outline = FALSE)

## PCA after normalization
PCA_new(shsox9_expr_norm, 
        ntop = nrow(shsox9_expr_norm),
        group = shsox9_raw_targets$group,
        show_name = T)

##### 2.Analysis of differences between two sets of unpaired samples Agilent data----
rm(list = ls())
#### 2.1 load data ----
source("./Custom_Functions.R")
load_input <- load("raw_data/shsox9_limma_processed.Rdata") #经过标准化处理的数据
load_input

#### 2.2 Organise the data ----
targets <- shsox9_raw_targets
table(targets$group)

## 过滤数据，探针过滤（Agilent数据中的gIsWellAboveBG）
y <- shsox9_norm #标准化后的数据
dim(y)
Control <- y$genes$ControlType == 1L#是剔除一些不用的探针
IsExpr <- rowSums(y$other$gIsWellAboveBG > 0) >= 3#1是高于背景值，0是和背景值差别不大。因为normal有12个所以大概选了个12
yfilt <- y[!Control & IsExpr, ]#!Control 取相反，非对照的
dim(yfilt)

## expr 提取表达矩阵
expr <- yfilt$E
rownames(expr) <- yfilt$genes$ProbeName #将yfilt$genes中的ProbeName赋值给expr作为行名
colnames(expr) <- shsox9_raw_targets$sample_name
dim(expr)

library(RColorBrewer)#R配色的一个包
colors = rep(brewer.pal(6, "Set2")[1:2], each = 3)
#箱线图
boxplot(expr, las = 2, col = colors, outline = F)
#PCA图
PCA_new(expr, 
        group = shsox9_raw_targets$group,
        show_name = F) 

## 芯片注释
shsox9_P2S <- read.table('./shsox9_annotations.txt',
                           na.strings = "",
                           sep = '', # '\t' | ','   ''代表没有什么分隔符
                           header = T)

expr_anno <- annotate_expr(expr, shsox9_P2S)

is.matrix(expr_anno)
expr_anno <- as.matrix(expr_anno)

#### 2.3 differential expression analysis  ----
library(limma)
design <- model.matrix(~ group, data = shsox9_raw_targets) #设计矩阵
colnames(design)
colnames(design) <- c("Con", "TestvsCon")
head(design)

# 线性拟合，表达矩阵信息和design信息
fit <- lmFit(expr_anno, design)
# 经验贝叶斯方法，计算moderad t的统计量
fit <- eBayes(fit,trend = T)
#多重检验校正和提取差异分析结果
shsox9_DEG <- topTable(fit, coef = "TestvsCon",number = Inf)
##此处矫正后的P值均大于0.05
shsox9_DEG_sig  <- shsox9_DEG  %>% 
  dplyr::filter(abs(logFC) > 0.2, P.Value < 0.05) #筛选出显著差异基因

summary(decideTests(fit))

shsox9_DEG_sig1 <- shsox9_DEG_sig %>% rownames_to_column('symbol')
shsox9_DEG_1 <- shsox9_DEG %>% rownames_to_column('symbol')

save(expr_anno,shsox9_DEG,shsox9_raw_targets,shsox9_DEG_sig,shsox9_DEG_1,shsox9_DEG_sig1,
     file = "raw_data/shsox9_anno_DEG.Rdata")


##### 3 Plot ----
#### 3.1 ARGs ----
library(tidyverse)
## 自噬基因
autophagy_gene <- read.csv('raw_data/AUT_gene.CSV')
rownames(autophagy_gene) <- autophagy_gene[,1]
colnames(autophagy_gene)[1]<-  'symbol'

deg <- shsox9_DEG_sig1
deg$symbol <- gsub("\\s+","",deg$symbol)#去掉字符串间的空格
colnames(deg)[1] <- 'symbol'
deg$symbol <- toupper(deg$symbol)
shsox9_AUT_DEGs <- deg %>% inner_join(autophagy_gene,by = "symbol") # 自噬相关差异表达基因

#筛选出|logFC|>0.5的auto_DEs
shsox9_sutoDEG_sig  <- shsox9_AUT_DEGs  %>% 
  dplyr::filter(abs(logFC) > 0.50, P.Value < 0.05) #筛选出显著差异基因

## AUT_DEG_ALL
shsox9_DEG_1$symbol <- gsub("\\s+","",shsox9_DEG_1$symbol)
shsox9_DEG_1$symbol <- toupper(shsox9_DEG_1$symbol)

shsox9_AUT_DEG_ALL <- shsox9_DEG_1 %>% inner_join(autophagy_gene,by = "symbol") #自噬相关差异表达基因

#### 3.2 volcano plot  ----
rm(list = ls())#一键清空
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
load("raw_data/shsox9_anno_DEG.Rdata")#加载差异分析结果
volcano<-shsox9_DEG 
#文件新建一列，用来显示基因的差异表达情况
volcano$type[(volcano$P.Val > 0.05|volcano$P.Val=="NA")|(volcano$logFC < 0.2)& volcano$logFC > -0.2] <- "none significant"
volcano$type[volcano$P.Val <= 0.05 & volcano$logFC >= 0.2] <- "up-regulated"
volcano$type[volcano$P.Val <= 0.05 & volcano$logFC <= -0.2] <- "down-regulated"
volcano <-volcano %>% rownames_to_column('symbol')
#画图
p = ggplot(volcano,aes(logFC,-1*log10(P.Value),color=type))
p + geom_point(alpha = 0.8) #散点图，alpha就是点的透明度
x_lim <- max(volcano$logFC,-volcano$logFC) #x轴的显示范围

gg=p + geom_point( size=1.3,alpha = 0.8) + xlim(-2,x_lim) +
  ylab(expression(-log[10]("P Value"))) +
  xlab(expression(log[2]("Fold Change"))) +
  scale_color_manual(values =c("#2f5688","#BBBBBB","#CC0000"))+
  theme_bw()+
  ggtitle('')+
  geom_vline(xintercept=c(-0.2,0.2),colour="black", linetype="dashed")+ # 加垂直线
  geom_hline(aes(yintercept=-1*log10(0.05)),colour="black", linetype="dashed")#加水平线
print(gg)


index<-c("Sox9","Ulk1","Pea15","Bag3","Gabarapl2",'Pten','Pelp1','Capns1','Raf1','Eef2')
select_gene<-volcano[which(volcano$symbol%in%index==T),]

kk= gg+geom_label_repel(data = select_gene, 
                 aes(label = as.character(symbol)), 
                 size = 3,box.padding = unit(0.6, 'lines'), 
                 segment.color = 'black',
                 show.legend = FALSE)+
                 theme_bw()  

 print(kk)

#### 3.3 heat map ----
library(pheatmap)
load('raw_data/shsox9_AUT_DEGs.Rdata')
load("raw_data/shsox9_anno_DEG.Rdata")
shsox9_AUT_DEGs=read.table("raw_data/shsox9_AUT_DEGs.txt",sep="\t",
                        header=T,check.names=F)

is.data.frame(expr_anno)
expr_anno <-as.data.frame(expr_anno)
expr_anno <- expr_anno %>% rownames_to_column("symbol")

## 取交集
##数据
library(Hmisc)#将首字母转换成大写
AUT_DEGs <- as.data.frame(shsox9_AUT_DEGs$symbol)
colnames(AUT_DEGs) <- 'symbol'
AUT_DEGs$symbol <- tolower(AUT_DEGs$symbol)#得先转成小写，才能转为大写
AUT_DEGs$symbol <- capitalize(AUT_DEGs$symbol)

AUT_DEGs[17,]<- "Sox9"#某列中添加一个元素

exprset<-expr_anno %>% inner_join(AUT_DEGs,by = "symbol") #取交集
rownames(exprset) <-exprset[,1]
exprset <-exprset[,-1]

#### 绘制热图 
annotation_col =data.frame(Group = factor(c(rep("Control",3),rep("shSox9", 3))))#选择group那一列
rownames(annotation_col) = colnames(exprset)

pheatmap(exprset,   #筛选出的基因表达矩阵
         scale = "row",#每个基因的水平scale一下,以行来标准化，这个功能很不错
         # border_color = NA, #热图的边界
         main = "shSox9_AUT_heatmap",#标题
         annotation_legend = T, #注释的图例
         annotation_col = annotation_col, #列的注释
         # annotation_row = annotation_row,
         clustering_method = "average",# clustering_method参数设定不同聚类方法，默认为"complete",可以设定为'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
         clustering_distance_rows = "euclidean",# clustering_distance_rows = "correlation"参数设定行聚类距离方法为Pearson corralation，默认为欧氏距离"euclidean"
         # color = colorRampPalette(c("navy", "white", "firebrick3"))(100),# color参数自定义颜色
         color = colorRampPalette(c("blue", "white","red"))(50),
         cluster_row = T,
         cluster_col = F,# cluster_COL = FALSE参数设定不对行进行聚类
         legend_labels = c("1.0","2.0","3.0","4.0","5.0"),
         legend = T,
         # show_rownames=T,
         # show_colnames = T,
         treeheight_row =10,
         #treeheight_col = 50,
         fontsize = 6,#基因名显示大小
         cellwidth = 15,
         cellheight = 8,
filename = "Heatmap_AUT_Sig.pdf")
dev.off()

#### 3.4 corrplot  ----
rm(list = ls())
rt=read.table("raw_data/shsox9_AUT_Exprset(top 5 genes).txt",sep="\t",
              header=T,check.names=F,row.names = 1)  #读取输入文件
fix(rt)    #查看数据（不错的方法）
rt=t(rt) #将数据行列转置

rt= round(cor(rt,method = "spearman"), 2) #保留2位小数
library(ggcorrplot)
library(ggthemes)
pdf("corHeatmap1.pdf",height=6,width=6)   
p.mat <- round(cor_pmat(rt),3) #计算p值

ggcorrplot(rt,
           hc.order = T,
           hc.method = "ward.D",
           ggtheme = theme_bw(),
           tl.cex = 15, 
           # type = 'upper',
           colors = c('#6D9EC1',"white","#E46726"),
           lab = T,
           lab_size = 3.5,
           # p.mat = p.mat,insig = "blank",
           outline.color = 'white'
)
dev.off()

#### 3.5 VennDiagram ----
library(VennDiagram)
library(Hmisc)#将首字母转换成大写
load("raw_data/shsox9_AUT_DEGs.Rdata")
autophagy_gene$symbol <- tolower(autophagy_gene$symbol)#得先转成小写，才能转为大写
autophagy_gene$symbol <- capitalize(autophagy_gene$symbol)

venn.diagram(list("ARGs" = autophagy_gene$symbol,
                  "DEGs" = shsox9_DEG_sig1$symbol),
             height=3000,
             width=3200,
             resolution=500,
             imagetype="tiff",#图像的格式
             filename="AUT_DEG_VennPlot.tiff",
             col="transparent",#指定图形的圆周边缘颜色 transparent 透明
             fill = c("cornflowerblue", "red"),
             alpha = 0.50,#透明度
             label.col = c("black", "black", "black"),
             cat.fontface = "bold",#字体的格式
             # cat.fontfamily = "serif",#字体的分类
             lty='blank',#区域边框线类型
             cat.pos = c(180, 190),#分类名称在圆的位置，默认正上方，通过角度进行调整
             cex=2,#圈内文字大小
             cat.cex=1.5,#圈外文字大小
             # cat.dist = 0.09, #分类名称距离边的距离（可以为负数）
             # cat.just = list(c(-1, -1), c(1, 1)), #分类名称的位置
             # ext.pos = 30, #线的角度 默认是正上方12点位置
             ext.dist = -0.04, #外部线的距离
             ext.length = 0.85, #外部线长度
             ext.line.lwd = 2, #外部线的宽度
             ext.line.lty = 'dashed') #外部线为虚线

#### KEGG ----
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

setwd("./KEGG")             
rt=read.table("id.txt",sep="\t",header=T,check.names=F)      
rt=rt[is.na(rt[,"entrezID"])==F,]          
gene=rt$entrezID

#kegg
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   
write.table(kk,file="KEGGId.txt",sep="\t",quote=F,row.names = F)                         

#barplot
pdf(file="barplot.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()

#bubble
pdf(file="bubble.pdf",width = 10,height = 7)
dotplot(kk, showCategory = 30)
dev.off()
#### GO ----
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

setwd("./GO")                  
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                
gene=rt$entrezID

#GO
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)       

#barplot
pdf(file="barplot.pdf",width = 10,height = 8)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#bubble
pdf(file="bubble.pdf",width = 10,height = 8)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()