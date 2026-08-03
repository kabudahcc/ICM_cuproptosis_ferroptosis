rm(list = ls())
#打破下载时间的限制,改前60秒，改后10w秒
options(timeout = 100000) 
options(scipen = 20)#不要以科学计数法表示

# GSE57338: 心衰 vs 健康，共 313 例
# 平台: GPL11532 ([HuGene-1_1-st] Affymetrix Human Gene 1.1 ST Array)
# 平台注释文件: GPL11532_annot.txt（已准备好在工作目录）

gse_id <- "GSE57338"

library(GEOquery)
eSet = getGEO(filename = paste0(gse_id, '_series_matrix.txt.gz'), destdir=".", getGPL = F)

class(eSet)
length(eSet)

#(1)提取表达矩阵exp
exp <- exprs(eSet)
dim(exp)
range(exp)
#exp = log2(exp+1) 
exp <- as.data.frame(exp)
boxplot(exp,las = 2)

#(2)提取临床信息
pd <- pData(eSet)

#(3)让exp列名与pd的行名顺序完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) {
  s = intersect(rownames(pd),colnames(exp))
  exp = exp[,s]
  pd = pd[s,]
}

save(pd,exp,file = paste0(gse_id, "_step1output.Rdata"))

# Group(实验分组)和ids(探针注释)
library(stringr)
library(dplyr)

#提取目标GSE Group-分组 ####
# GSE57338 标题特征: "Non-failing" 为对照，其余（Ischemic / Idiopathic Dilated CMP）为疾病
GSE_group <- pd %>%
  dplyr::select(title, geo_accession) %>%
  mutate(group = if_else(str_detect(title, regex("Non-failing", ignore_case = TRUE)) ,
                         "Normal","Disease")) %>%
  arrange(group, geo_accession)
head(GSE_group)

Group = factor(GSE_group$group,levels = c("Normal","Disease"))
Group

write.csv(GSE_group, paste0(gse_id, '_group.csv'))

#2.探针注释的获取-----------------
library(data.table)
gpl_number <- eSet@annotation;gpl_number
# 平台注释文件（已下载/解析为 ID + Gene_Symbol 两列）
ann <- fread('GPL11532_annot.txt',header = T,sep = '\t')

colnames(ann)
ann <- ann %>% select(ID, Gene_Symbol)
colnames(ann) <- c('id','symbol')
ann$id <- as.character(ann$id)

exp$id <- rownames(exp)

GSE_ann <- inner_join(exp,ann,by = 'id')

GSE_ann <- GSE_ann %>% 
  select(symbol, everything(), -id)

GSE_ann <- GSE_ann %>%
  filter(symbol != '')

x=GSE_ann$symbol
a1=strsplit(x,split = " /// ",fixed = T)
gene.all = sapply(a1,function(x){x[1]})
GSE_ann$symbol=gene.all

GSE_ann <- GSE_ann %>% 
  group_by(symbol) %>%
  summarise_all(mean) %>%
  data.frame()

row.names(GSE_ann) <- GSE_ann$symbol
GSE_ann <- GSE_ann [,-1]

write.csv(GSE_ann, paste0(gse_id, '_symbol.csv'))
save(exp,GSE_group,Group,GSE_ann,file = paste0(gse_id, "_step2_GEO_data.Rdata"))

# 3.PCA 图----
dat=as.data.frame(t(GSE_ann))
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(dat, graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = Group,
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE,
             legend.title = "Groups"
)

# 4.top 1000 sd 热图---- 
g = names(tail(sort(apply(GSE_ann,1,sd)),1000)) 
n = GSE_ann[g,]
library(pheatmap)
annotation_col = data.frame(row.names = colnames(n),
                            Group = Group)
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row",
         breaks = seq(-3,3,length.out = 100)
) 

#5.差异分析 #####
expr <- GSE_ann
expr <- expr[, match(GSE_group$geo_accession, colnames(expr))]

boxplot(expr, outline=F, notch=T, col=Group,las = 2,main = gse_id)

all(colnames(expr) == GSE_group$geo_accession)

contrast <- paste0(rev(levels(Group)), collapse = "-")
contrast

design <- model.matrix( ~ 0 + Group)
colnames(design) <- levels(Group)
design

library(limma)
contrast.matrix <- makeContrasts(contrast, levels = design)
contrast.matrix

fit <- lmFit(expr, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

library(tidyverse)
GSE_DEG <- topTable(fit, coef = 1, n = Inf) %>% 
  rownames_to_column(var = "symbol") 

GSE_DEG_p <- GSE_DEG %>% 
  dplyr::filter(abs(logFC) > 1, P.Value < 0.05)

write.table(GSE_DEG, file = paste0(gse_id, '_diff.txt'), sep = '\t', quote = F,row.names = T, col.names = T)

##6.火山图##
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(data.table)
data <- fread(paste0(gse_id, '_diff.txt'),sep = '\t')
data <-data[,-1]
logFCfilter = 0.5
logFCcolor = 1
colnames(data)[2] <- 'logFC'
colnames(data)[1] <- 'gene'

index = data$P.Value <0.05 & abs(data$logFC) > logFCfilter
data$group <- 0
data$group[index & data$logFC>0] = 1
data$group[index & data$logFC<0] = -1
data$group <- factor(data$group,levels = c(1,0,-1),labels =c("High","NS","Low") )

### 正式画图
ggplot(data=data, aes(x=logFC, y =-log10(P.Value),color=group)) +
  geom_point(alpha=0.8, size=1.2)+
  scale_color_manual(values = c("#EFC000", "grey50", "#0073C2"))+
  labs(x="log2 (fold change)",y="-log10 (P.Value)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(-logFCfilter,logFCfilter),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(legend.position="top")+
  geom_point(data=subset(data, abs(logFC) >= logFCcolor & P.Value <0.05),alpha=0.8, size=1,col="green4")+
  geom_text_repel(data=subset(data, abs(logFC) >= logFCcolor & P.Value <0.05),
                  aes(label=gene),col="black",alpha = 0.8)

#7. top10差异基因热图
if(T){
  library(dplyr)
  dat2 = data %>%
    filter(group!="NS") %>% 
    arrange(logFC) 
  cg = c(head(dat2$gene,10),
         tail(dat2$gene,10))
}else{
  cg = data$gene[data$group !="NS"]
  length(cg)
}

n=expr[cg,]
dim(n)

library(pheatmap)
annotation_col = GSE_group["group"]
heatmap_plot <- pheatmap(n,show_colnames =F,
                         scale = "row",
                         annotation_col=annotation_col,
                         breaks = seq(-3,3,length.out = 100)
) 
heatmap_plot
