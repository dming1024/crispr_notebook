#复现GSE126486中的分析
#对sgRNA进行QC
setwd("D:\\writing\\20240301_CRISPR\\GSE126486")
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

library(pheatmap)
library('factoextra')
#library("FactoMineR")
library(ggrepel)

#ET_r1: D0 replicate1
#ET_r2: D0 replicate2
#D18_DMSO_r1: DMSO D18 replicate1
#D18_DMSO_r2: DMSO D18 replicate2
#D18_CompB_r1: compoubB D18 replicate1
#D18_CompB_r2: compoubB D18 replicate2
#PP_r1: plasmid replicate1
#PP_r2: plasmid replicate2


sgRNA_count<-fread("GSE126486_sgRNA_read_counts.txt/GSE126486_sgRNA_read_counts.txt") %>% 
  dplyr::select(c(5:12))

p1=sgRNA_count %>% 
  tidyr::pivot_longer(everything(),names_to = 'Samples',values_to = 'Reads') %>% 
  ggplot(aes(x=Samples,y=log2(Reads+1)))+
  geom_boxplot(fill="#1874CD")+
  labs(x="",y="Log2 transformed sgRNA counts")+
  theme_bw(
    base_size = 15
  )+
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
  )
ggsave('results/fig1.png',p1,width = 10,height = 5)

p2=sgRNA_count %>%
  tidyr::pivot_longer(everything(),names_to = 'Samples',values_to = 'Reads') %>% 
  ggplot(aes(Reads,colour=Samples))+
  stat_density(geom='line',position = 'identity')+
  labs(x="Reads Count",y='Frequency',title = 'Density plot of sgRNA reads')+
  theme_bw(
    base_size = 15
  )+
  theme(
    legend.title=element_blank(),
    legend.position =  c(0.8,0.5),
    legend.background = element_blank()
  )
ggsave('results/fig2.png',p2,width = 10,height = 5)


p3=sgRNA_count %>% 
  tidyr::pivot_longer(everything(),names_to = 'Samples',values_to = 'Reads') %>% 
  ggplot(aes(Reads,colour=Samples))+
  stat_ecdf()+
  labs(x="Reads Count",y='Likelihood Probability',title = 'Cumulative distribution')+
  theme_bw(
    base_size = 15
  )+
  theme(
    legend.title=element_blank(),
    legend.position =  c(0.8,0.5),
    legend.background = element_blank()
  )+
  xlim(c(1,500))

ggsave('results/fig3.png',p3,width = 10,height = 5)

#样本相关性分析
mx=cor(sgRNA_count)
pheatmap(
  mx,
  cluster_rows = T,
  cluster_cols = T,
  treeheight_row = 0,
  show_rownames = T,
  show_colnames = F,
  display_numbers = T,
  number_format = "%.2f",
  number_color = 'white',
  fontsize_number = 15,
  filename = 'results/fig7.png'
)

#PCA分析
sample_info=data.frame(Samples=colnames(sgRNA_count)) %>% 
  mutate(Group=gsub("(.*)\\_.*","\\1",Samples))
res.pca <- prcomp(t(sgRNA_count))
#分组作图分析
p1=fviz_pca_ind(
  res.pca,
  geom.ind = c("point","text"),
  alpha.ind = 0.5,
  repel = T
)+geom_point(aes(colour=factor(sample_info$Group)))+
  guides(colour = guide_legend(title = "Group"))
ggsave('results/fig8.png',p1,width = 6,height = 5)



#使用chronos计算各样本中的gene effect值
fread("GSE126486_sgRNA_read_counts.txt/GSE126486_sgRNA_read_counts.txt") %>% select(c(2,1)) %>% 
  setnames(c("sgrna","gene")) %>% 
  write.csv("guide_map.csv",row.names = F)
tmp=fread("GSE126486_sgRNA_read_counts.txt/GSE126486_sgRNA_read_counts.txt") %>% select(c(5:10)) %>% as.data.frame()
rownames(tmp)=fread("GSE126486_sgRNA_read_counts.txt/GSE126486_sgRNA_read_counts.txt") %>% pull(2)
tmp %>% t() %>% as.data.frame() %>% write.csv(.,"sgRNA_counts.csv")
tmp %>%select(c(3:6)) %>%  t() %>% as.data.frame() %>% write.csv(.,"sgRNA_counts_2_vs_2.csv")
#D:\writing\20240301_CRISPR\GSE126486\results，使用python脚本实现

#essential QC
#geneTypes=readxl::read_xlsx("results//4_types_of_genes.xlsx") %>% setnames(c("Gene","geneType"))
geneTypes=fread("results/human.CRISPR_gene_type.depmap22Q2.csv")%>% setnames(c("Gene","geneType"))
tmp1=fread("results/AchillesCommonEssentialControls.csv") %>% mutate(geneType="Essential")
tmp2=fread("results/AchillesNonessentialControls.csv") %>% mutate(geneType="Nonessential")
geneTypes=rbind(tmp1,tmp2)
chronos_gene_effect<-fread("results/screen/run/gene_effects.csv")
p4=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  left_join(.,geneTypes) %>% 
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>% 
  ggplot(aes(Effect,colour=Classification))+
  stat_density(geom='line',position = 'identity')+
  facet_wrap(~cell_line_name)+
  labs(title="Chronos")+
  theme_bw(
    base_size = 15
  )
ggsave('results/fig4.png',p4,width = 10,height = 5)


p5=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  left_join(.,geneTypes) %>% 
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>% 
  ggplot(aes(x=Classification,y=Effect,fill=cell_line_name))+
  geom_boxplot()+
  theme_bw(
    base_size = 15
  )+
  theme(
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1)
  )
ggsave('results/fig5.png',p5,width = 10,height = 5)


#rank gene-effect
library(ggrepel)
tmp=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  left_join(.,geneTypes) %>% 
  mutate(geneType=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) 
tmp %>% group_by(cell_line_name) %>%
  arrange(desc(Effect)) %>% 
  mutate(rank_ID=1:n()) %>% 
  slice_min(.,Effect,n=5) -> top5_essential_genes
p10_A=tmp %>% group_by(cell_line_name) %>% 
  arrange(desc(Effect)) %>% 
  mutate(rank_ID=1:n()) %>% 
  ggplot(aes(x=rank_ID,y=Effect))+
  geom_point(aes(colour=geneType))+
  geom_text_repel(aes(x=rank_ID,y= Effect,label=Gene),data=top5_essential_genes)+
  facet_wrap(~cell_line_name)+  
  labs(title="Rank Effect")+
  theme_bw(
    base_size = 15
  )+
  theme(
    legend.position = 'bottom'
  )+
  guides(colour = guide_legend(nrow = 2))

tmp %>% group_by(cell_line_name) %>%
  filter(geneType=="other") %>% 
  arrange(desc(Effect)) %>% 
  mutate(rank_ID=1:n()) %>% 
  slice_min(.,Effect,n=5) -> top5_essential_genes
p10_B=tmp %>% group_by(cell_line_name) %>%  
  filter(geneType=="other") %>% 
  arrange(desc(Effect)) %>% 
  mutate(rank_ID=1:n()) %>% 
  ggplot(aes(x=rank_ID,y=Effect))+
  geom_point(aes(colour=geneType))+
  geom_text_repel(aes(x=rank_ID,y= Effect,label=Gene),data=top5_essential_genes)+
  facet_wrap(~cell_line_name)+  
  labs(title="Rank Effect")+
  theme_bw(
    base_size = 15
  )+
  theme(
    legend.position = 'bottom'
  )+
  guides(colour = guide_legend(nrow = 2))

library(ggpubr)
ps=ggpubr::ggarrange(p10_A,p10_B,nrow = 1)
ggsave("results//fig6_rank_geneEffect.png",ps,width = 10,height = 6)


#交叉比较分析
#同一化合物，不同时间点比较
tmp %>% filter(geneType!="Essential") %>% 
  select(c(cell_line_name,Gene,Effect)) %>% 
  tidyr::pivot_wider(names_from = cell_line_name, values_from = c(Effect)) %>% 
  setnames(c("Gene","y","x")) -> regression_data
#计算回归方程的系数
alphas=seq(0.1,10,by=0.1)
absolute_erros = lapply(
  alphas,
  FUN = function(alpha) {
    sum_absolute_error = regression_data %>%
      #计算垂直距离
      mutate(y_dist = abs(y - x * alpha)*(1/sqrt(1+alpha^2))) %>%
      select(y_dist) %>% sum()
  }
)
candidate_alpha=alphas[which.min(absolute_erros)]

#结果可视化，采用6倍SD，可以根据项目改动
sigma_sd = regression_data %>% mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) ) %>% pull(y_dist) %>% sd()*6
library(ggrepel)
p1=regression_data %>% 
  ggplot(aes(x=x,y=y))+
  geom_point()+
  geom_point(aes(x=x,y=y),colour='grey',
             data=regression_data %>% 
               mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
               filter(abs(y_dist) < sigma_sd)
  )+
  geom_text_repel(aes(x=x,y=y,label=Gene),max.overlaps = 15,
                  data=regression_data %>% 
                    mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
                    filter(abs(y_dist) > sigma_sd)
  )+
  geom_abline(slope = candidate_alpha,
              intercept = 0,colour='red',lwd=1.2)+
  #计算上部截距
  geom_abline(slope = candidate_alpha,
              intercept = c(sigma_sd/(1/sqrt(1+candidate_alpha^2))),
              colour='grey',lty=2,lwd=1.2)+
  #计算下部截距
  geom_abline(slope = candidate_alpha,
              intercept = c(- sigma_sd/(1/sqrt(1+candidate_alpha^2))),
              colour='grey',lty=2,lwd=1.2)+
  labs(x="DMSO",y="CompB")+
  theme_bw(
    base_size = 15
  )+xlim(c(-3,1))+ylim(c(-3,1))
ggsave("results/fig6_cross_comparisions_withoutessenGenes.jpg",p1,width = 8,height = 6)

regression_data %>% mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) ) %>% 
  ggplot(aes(x=y_dist))+
  stat_density(geom='line',position = 'identity')

#一些筛选结果
regression_data %>% 
  mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
  filter(abs(y_dist) > sigma_sd) %>% arrange(desc(abs(y_dist))) %>% write.csv(.,"results/dependency_genes.csv")
regression_data %>% 
  mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
  filter(abs(y_dist) > sigma_sd) %>% 
  filter(abs(x)<0.3 & abs(y)>0.8)

#对这些dependency genes进行功能注释，看一看在什么通路上，单个基因也可以
#富集结果没什么参考意义


paper_report=c("RRAGA",'MLST8','WDR24','LAMTOR2')
paper_report=c("RRAGA",'MLST8','WDR24','LAMTOR2',
               "RREB1","GFI1","AKT1")
regression_data %>% 
  mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
  filter(Gene %in% paper_report)
