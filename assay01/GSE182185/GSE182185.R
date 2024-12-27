#复现GSE182185中的分析
#对sgRNA进行QC
setwd("D:\\writing\\20240301_CRISPR\\GSE182185")
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

library(pheatmap)
library('factoextra')
#library("FactoMineR")
library(ggrepel)

#3BT01, p0?
#3BT1, DMSO
#3BAbe, Abemaciclib
#3BPal, Palbociclib
#3BLenva, Lenvatinib
#3BReg, Regorafenib
#3BSora, Sorafenib
#BEra, Erastin


sgRNA_count<-fread("data/GSE182185_SNU-398.txt") %>% 
  select(-c(1,2))

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
  )+xlim(c(0,500))
ggsave('results/fig2.png',p2,width = 10,height = 5)#zero counts太多了，不适合chronos，细胞死亡的多


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

#zero counts占比太高了，细胞死亡比例高，不适用chronos，尝试MAGECK或Drug-Z
p4=data.frame(zero=apply(sgRNA_count,2,FUN=function(x){round(table(x==0)[2] / length(x) ,4)*100})) %>% 
  mutate(samples=rownames(.)) %>% 
  mutate(group=gsub("(.*)\\-.*","\\1",samples)) %>% 
  ggplot(aes(x=samples,y=zero))+
  geom_bar(stat = 'identity',aes(fill=group))+
  theme_bw(base_size = 15)+
  theme(
    axis.text.x  = element_text(angle = 45),
    legend.position = 'none'
  )+
  labs(x='',y='Zero Count Ratio (%)')
ggsave('results/fig4.png',p4,width = 10,height = 5)

#样本相关性分析
mx=cor(sgRNA_count)
pheatmap(
  mx,
  cluster_rows = T,
  cluster_cols = T,
  treeheight_row = 0,
  show_rownames = T,
  show_colnames = F,
#  display_numbers = T,
 # number_format = "%.2f",
  number_color = 'white',
  fontsize_number = 15,
  filename = 'results/fig7.png'
)

#PCA分析
sample_info=data.frame(Samples=colnames(sgRNA_count)) %>% 
  mutate(Group=gsub("(.*)\\-.*","\\1",Samples))
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


#RRA，基于RRA分析结果进行可视化
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
hep38_abemaciclib_rra=fread("results/Hep3B/rra/Hep3B_Abemaciclib_vs_DMSO.gene_summary.txt") %>% 
  select(c(1,3,5,8)) %>% setnames(c('Gene','score','fdr','lfc'))
hep38_abemaciclib_drugz=fread("results/Hep3B/drugz_results/Hep3B_Abemaciclib_vs_DMSO.txt") %>% 
  select(1,4,6,7) %>% setnames(c('Gene','normZ','rank','fdr_synth'))

paper_report=c('SOD2','PSTK')
hep38_abemaciclib_rra %>% 
  left_join(.,hep38_abemaciclib_drugz) %>% 
  left_join(.,geneTypes) %>% 
  filter(Gene %in% paper_report)  %>% 
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) -> candidate1

p6=hep38_abemaciclib_rra %>% 
  left_join(.,hep38_abemaciclib_drugz) %>% 
  left_join(.,geneTypes) %>% 
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>%
  ggplot(aes(x=lfc,y=normZ))+
  geom_point(aes(colour=Classification))+
  geom_text_repel(aes(x=lfc,y= normZ,label=Gene),data=candidate1)+
  geom_point(aes(x=lfc,y= normZ,size=4),shape=21,data=candidate1)+
  facet_wrap(~Classification)+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'none')+
  labs(x='RRA LogFC',y='DrugZ score')
ggsave('results/fig6.png',p6,width = 10,height = 4)

hep38_abemaciclib_rra %>% 
  left_join(.,hep38_abemaciclib_drugz) %>% 
  left_join(.,geneTypes) %>% 
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>% write.csv(.,"results/hep3b_rra_drugz_results.csv")
hep38_abemaciclib_rra %>% 
  left_join(.,hep38_abemaciclib_drugz) %>% 
  left_join(.,geneTypes) %>% 
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>% 
  filter(normZ< -2.5 & lfc < -2) -> candidate1
p9=hep38_abemaciclib_rra %>% 
  left_join(.,hep38_abemaciclib_drugz) %>% 
  left_join(.,geneTypes) %>% 
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>%
  filter(Classification=='other') %>% ggplot(aes(x=lfc,y=normZ))+
  geom_point()+
  geom_text_repel(aes(x=lfc,y= normZ,label=Gene),data=candidate1)+
  geom_point(aes(x=lfc,y= normZ),color='darkred',data=candidate1)+
  geom_rect(aes(xmin = -7.5, xmax = -2, ymin = -5.5, ymax = -2.5), color='black', 
            linetype='dashed', alpha=0) 
ggsave('results/fig9.png',p9,width = 8,height = 6)


#联合RRA 和 drugz结果，按细胞寻找overlap target genes

# #essential QC，需要和p0进行比较，然后再可视化
geneTypes=fread("results/human.CRISPR_gene_type.depmap22Q2.csv")%>% setnames(c("Gene","geneType"))
tmp1=fread("results/AchillesCommonEssentialControls.csv") %>% mutate(geneType="Essential")
tmp2=fread("results/AchillesNonessentialControls.csv") %>% mutate(geneType="Nonessential")
geneTypes=rbind(tmp1,tmp2)
drug_report=c('Sorafenib','Erastin','Abemaciclib',
              'Lenvatinib','Palbociclib','Regorafenib')
# hep38_abemaciclib_drugz %>% 
#   filter(Gene %in% paper_report) -> candidate1
# tmp=fread("results/Hep3B/drugz_results/Hep3B_Abemaciclib_vs_DMSO.txt") %>% 
#   filter(pval_synth<0.05) %>% select(1,4,6,7) %>% setnames(c('Gene','normZ','rank','fdr_synth')) 
tmps=lapply(
  drug_report,
  FUN=function(x){
    tmp=fread(sprintf("results/Hep3B/drugz_results_p0/Hep3B_%s_vs_p0.txt",x))%>% 
      #filter(pval_synth<0.05) %>% 
      select(1,4,6,7) %>% 
      mutate(drugs=x) %>% 
      setnames(c('Gene','normZ','rank','fdr_synth','drugs')) 
  }
)
tmpdf=rbindlist(tmps)
p4.1=tmpdf %>% select(c(Gene,normZ,drugs)) %>%
  left_join(.,geneTypes) %>%
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>%
  ggplot(aes(x=Classification,y=normZ))+
  geom_boxplot(aes(fill=Classification))+
  facet_wrap(~drugs)+
  labs(title="DrugZ")+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30))
ggsave('results/fig4_1.png',p4.1,width = 7,height = 5)
p4.3=tmpdf %>% select(c(Gene,normZ,drugs)) %>%
  left_join(.,geneTypes) %>%
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>%
  #ggplot(aes(x=Classification,y=lfc))+
  #geom_boxplot(aes(fill=Classification))+
  ggplot(aes(normZ,colour=Classification))+
  stat_density(geom='line',position = 'identity')+
  facet_wrap(~drugs)+
  labs(title="DrugZ")+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'bottom')
ggsave('results/fig4_3.png',p4.3,width = 7,height = 5)


#获取RRA分析结果
drug_report=c('Sorafenib','Erastin','Abemaciclib',
              'Lenvatinib','Palbociclib','Regorafenib')
tmps=lapply(
  drug_report,
  FUN=function(x){
    tmp=fread(sprintf("results/Hep3B/rra_p0/Hep3B_%s_vs_p0.gene_summary.txt",x))%>% 
      select(1,3,4,8) %>% setnames(c('Gene','score','pvalue','lfc')) %>% 
      #filter(pvalue<0.05) %>% 
      mutate(drugs=x) 
  }
)
tmpdf=rbindlist(tmps)
p4.2=tmpdf %>% select(c(Gene,lfc,drugs)) %>%
  left_join(.,geneTypes) %>%
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>%
  ggplot(aes(x=Classification,y=lfc))+
  geom_boxplot(aes(fill=Classification))+
  facet_wrap(~drugs)+
  labs(title="RRA")+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30))
ggsave('results/fig4_2.png',p4.2,width = 7,height = 5)

p4.4=tmpdf %>% select(c(Gene,lfc,drugs)) %>%
  left_join(.,geneTypes) %>%
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>%
  #ggplot(aes(x=Classification,y=lfc))+
  #geom_boxplot(aes(fill=Classification))+
  ggplot(aes(lfc,colour=Classification))+
  stat_density(geom='line',position = 'identity')+
  facet_wrap(~drugs)+
  labs(title="RRA")+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'bottom')
ggsave('results/fig4_4.png',p4.4,width = 7,height = 5)
#stat_density(geom='line',position = 'identity')





paper_report=c('SOD2','PSTK','SPR')
drug_report=c('Sorafenib','Erastin','Abemaciclib')
hep38_abemaciclib_drugz %>% 
  filter(Gene %in% paper_report) -> candidate1
# tmp=fread("results/Hep3B/drugz_results/Hep3B_Abemaciclib_vs_DMSO.txt") %>% 
#   filter(pval_synth<0.05) %>% select(1,4,6,7) %>% setnames(c('Gene','normZ','rank','fdr_synth')) 
tmps=lapply(
  drug_report,
  FUN=function(x){
    tmp=fread(sprintf("results/Hep3B/drugz_results/Hep3B_%s_vs_DMSO.txt",x))%>% 
      filter(pval_synth<0.05) %>% select(1,4,6,7) %>% 
      mutate(drugs=x) %>% 
      setnames(c('Gene','normZ','rank','fdr_synth','drugs')) 
  }
)
tmpdf=rbindlist(tmps)
tmpdf %>% 
  filter(Gene %in% paper_report) -> candidate1
drug_report=c('Sorafenib','Erastin','Abemaciclib')
p5=tmpdf %>% 
  #mutate(drugs = drugs %>% forcats::fct_relevel(drug_report)) %>%
  ggplot(aes(x=rank,y=normZ))+
  geom_point()+
  geom_text_repel(aes(x=rank,y= normZ,label=Gene),color='red',data=candidate1)+
  geom_point(aes(x=rank,y= normZ,size=4),shape=21,data=candidate1)+
  facet_wrap(~factor(drugs,drug_report),scales = 'free_y')+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'none')+
  labs(x='',y='DrugZ score')
ggsave('results/fig5.png',p5,width = 10,height = 4)
hep3B_drugz=tmpdf
#获取RRA分析结果
drug_report=c('Sorafenib','Erastin','Abemaciclib')
# tmp=fread("results/Hep3B/rra/Hep3B_Abemaciclib_vs_DMSO.gene_summary.txt") %>% 
#   select(1,3,4,8) %>% setnames(c('Gene','score','pvalue','lfc')) %>% filter(pvalue<0.05)
tmps=lapply(
  drug_report,
  FUN=function(x){
    tmp=fread(sprintf("results/Hep3B/rra/Hep3B_%s_vs_DMSO.gene_summary.txt",x))%>% 
      select(1,3,4,8) %>% setnames(c('Gene','score','pvalue','lfc')) %>% filter(pvalue<0.05) %>% 
      mutate(drugs=x) 
  }
)
tmpdf=rbindlist(tmps)
tmpdf %>% 
  group_by(drugs) %>%
  arrange(desc(lfc)) %>% 
  mutate(rank_ID=1:n()) %>%
  filter(Gene %in% paper_report) -> candidate1
drug_report=c('Sorafenib','Erastin','Abemaciclib')
p5.1=tmpdf %>% 
  group_by(drugs) %>%
  arrange(desc(lfc)) %>% 
  mutate(rank_ID=1:n()) %>%
  #mutate(drugs = drugs %>% forcats::fct_relevel(drug_report)) %>%
  ggplot(aes(x=rank_ID,y=lfc))+
  geom_point()+
  geom_text_repel(aes(x=rank_ID,y= lfc,label=Gene),color='red',data=candidate1)+
  geom_point(aes(x=rank_ID,y= lfc,size=4),shape=21,data=candidate1)+
  facet_wrap(~factor(drugs,drug_report),scales = 'free_y')+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'none')+
  labs(x='',y='RRA LogFC')
ggsave('results/fig5_1.png',p5.1,width = 10,height = 4)

hep3B_rra=tmpdf
#联合RRA与DrugZ的分析结果
hep3B_rra#过滤条件：pvalue<0.05
hep3B_drugz#过滤条件：pvalue<0.05
hep3B_rra %>% mutate(ID=paste0(Gene,"_",drugs)) %>% 
  left_join(
    .,hep3B_drugz %>% mutate(ID=paste0(Gene,"_",drugs))
  ) %>% filter(
    normZ< -2.5
  ) %>% 
  filter(
    lfc < -1
  ) %>% count(Gene) %>% filter(n>2)


x1=hep3B_rra %>% mutate(ID=paste0(Gene,"_",drugs)) %>% 
  left_join(
    .,hep3B_drugz %>% mutate(ID=paste0(Gene,"_",drugs))
  ) %>% 
  left_join(.,geneTypes) %>% 
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>% 
  filter(Classification=='other')
candidate1=x1 %>% filter(normZ< -2.5) %>%filter(lfc < -1 )
candidate2=x1 %>% filter(Gene %in% c('PSTK','SOD2','UPF3A'))
p10=x1 %>% 
  ggplot(aes(x=lfc,y=normZ))+
  geom_point()+
  #geom_text_repel(aes(x=lfc,y= normZ,label=Gene),data=candidate1)+
  geom_point(aes(x=lfc,y= normZ),color='darkred',data=candidate1)+
  geom_text_repel(aes(x=lfc,y= normZ,label=Gene),color='blue',data=candidate2)+
  geom_point(aes(x=lfc,y= normZ,size=4),shape=21,data=candidate2)+
  geom_rect(aes(xmin = -7.5, xmax = -1, ymin = -5.5, ymax = -2.5), color='black', 
            linetype='dashed', alpha=0) +
  facet_wrap(~drugs,scales = 'free_y')+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'none')+
  labs(x='RRA LogFC',y='DrugZ score')
ggsave('results/fig10.png',p10,width = 10,height = 5)  
  
  

# #rank gene-effect
# library(ggrepel)
# tmp=chronos_gene_effect %>% 
#   tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
#   left_join(.,geneTypes) %>% 
#   mutate(geneType=ifelse(is.na(geneType),
#                                ifelse(grepl("Control",Gene),
#                                       'Non-target','other'),
#                                geneType)) 
# tmp %>% group_by(cell_line_name) %>%
#   arrange(desc(Effect)) %>% 
#   mutate(rank_ID=1:n()) %>% 
#   slice_min(.,Effect,n=5) -> top5_essential_genes
# p10_A=tmp %>% group_by(cell_line_name) %>% 
#   arrange(desc(Effect)) %>% 
#   mutate(rank_ID=1:n()) %>% 
#   ggplot(aes(x=rank_ID,y=Effect))+
#   geom_point(aes(colour=geneType))+
#   geom_text_repel(aes(x=rank_ID,y= Effect,label=Gene),data=top5_essential_genes)+
#   facet_wrap(~cell_line_name)+  
#   labs(title="Rank Effect")+
#   theme_bw(
#     base_size = 15
#   )+
#   theme(
#     legend.position = 'bottom'
#   )+
#   guides(colour = guide_legend(nrow = 2))
# 
# tmp %>% group_by(cell_line_name) %>%
#   filter(geneType=="other") %>% 
#   arrange(desc(Effect)) %>% 
#   mutate(rank_ID=1:n()) %>% 
#   slice_min(.,Effect,n=5) -> top5_essential_genes
# p10_B=tmp %>% group_by(cell_line_name) %>%  
#   filter(geneType=="other") %>% 
#   arrange(desc(Effect)) %>% 
#   mutate(rank_ID=1:n()) %>% 
#   ggplot(aes(x=rank_ID,y=Effect))+
#   geom_point(aes(colour=geneType))+
#   geom_text_repel(aes(x=rank_ID,y= Effect,label=Gene),data=top5_essential_genes)+
#   facet_wrap(~cell_line_name)+  
#   labs(title="Rank Effect")+
#   theme_bw(
#     base_size = 15
#   )+
#   theme(
#     legend.position = 'bottom'
#   )+
#   guides(colour = guide_legend(nrow = 2))
# 
# library(ggpubr)
# ps=ggpubr::ggarrange(p10_A,p10_B,nrow = 1)
# ggsave("results//fig6_rank_geneEffect.png",ps,width = 10,height = 6)
# 
# 
# #交叉比较分析
# #同一化合物，不同时间点比较
# tmp %>% filter(geneType!="Essential") %>% 
#   select(c(cell_line_name,Gene,Effect)) %>% 
#   tidyr::pivot_wider(names_from = cell_line_name, values_from = c(Effect)) %>% 
#   setnames(c("Gene","y","x")) -> regression_data
# #计算回归方程的系数
# alphas=seq(0.1,10,by=0.1)
# absolute_erros = lapply(
#   alphas,
#   FUN = function(alpha) {
#     sum_absolute_error = regression_data %>%
#       #计算垂直距离
#       mutate(y_dist = abs(y - x * alpha)*(1/sqrt(1+alpha^2))) %>%
#       select(y_dist) %>% sum()
#   }
# )
# candidate_alpha=alphas[which.min(absolute_erros)]
# 
# #结果可视化，采用6倍SD，可以根据项目改动
# sigma_sd = regression_data %>% mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) ) %>% pull(y_dist) %>% sd()*6
# library(ggrepel)
# p1=regression_data %>% 
#   ggplot(aes(x=x,y=y))+
#   geom_point()+
#   geom_point(aes(x=x,y=y),colour='grey',
#              data=regression_data %>% 
#                mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
#                filter(abs(y_dist) < sigma_sd)
#   )+
#   geom_text_repel(aes(x=x,y=y,label=Gene),max.overlaps = 15,
#                   data=regression_data %>% 
#                     mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
#                     filter(abs(y_dist) > sigma_sd)
#   )+
#   geom_abline(slope = candidate_alpha,
#               intercept = 0,colour='red',lwd=1.2)+
#   #计算上部截距
#   geom_abline(slope = candidate_alpha,
#               intercept = c(sigma_sd/(1/sqrt(1+candidate_alpha^2))),
#               colour='grey',lty=2,lwd=1.2)+
#   #计算下部截距
#   geom_abline(slope = candidate_alpha,
#               intercept = c(- sigma_sd/(1/sqrt(1+candidate_alpha^2))),
#               colour='grey',lty=2,lwd=1.2)+
#   labs(x="DMSO",y="CompB")+
#   theme_bw(
#     base_size = 15
#   )+xlim(c(-3,1))+ylim(c(-3,1))
# ggsave("results/fig6_cross_comparisions_withoutessenGenes.jpg",p1,width = 8,height = 6)
# 
# regression_data %>% mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) ) %>% 
#   ggplot(aes(x=y_dist))+
#   stat_density(geom='line',position = 'identity')
# 
# #一些筛选结果
# regression_data %>% 
#   mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
#   filter(abs(y_dist) > sigma_sd) %>% arrange(desc(abs(y_dist))) %>% write.csv(.,"results/dependency_genes.csv")
# regression_data %>% 
#   mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
#   filter(abs(y_dist) > sigma_sd) %>% 
#   filter(abs(x)<0.3 & abs(y)>0.8)
# 
# #对这些dependency genes进行功能注释，看一看在什么通路上，单个基因也可以
# #富集结果没什么参考意义

