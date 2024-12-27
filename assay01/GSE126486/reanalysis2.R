
#reanalysis 2
#直接将compound 与 DMSO 进行比较, 可能会与RRA的结果更接近
setwd("D:\\writing\\20240301_CRISPR\\GSE126486")
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

library(pheatmap)
library('factoextra')
#library("FactoMineR")
library(ggrepel)

geneTypes=fread("results/human.CRISPR_gene_type.depmap22Q2.csv")%>% setnames(c("Gene","geneType"))
tmp1=fread("results/AchillesCommonEssentialControls.csv") %>% mutate(geneType="Essential")
tmp2=fread("results/AchillesNonessentialControls.csv") %>% mutate(geneType="Nonessential")
geneTypes=rbind(tmp1,tmp2)

chronos_gene_effect<-fread("results/screen/run2/gene_effects.csv")
p1=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  left_join(.,geneTypes) %>% 
  mutate(Classification=ifelse(is.na(geneType),
                               ifelse(grepl("Control",Gene),
                                      'Non-target','other'),
                               geneType)) %>% 
  ggplot(aes(Effect,colour=Classification))+
  stat_density(geom='line',position = 'identity')+
  labs(title="Chronos")+
  theme_bw(
    base_size = 15
  )

ggsave('results/run2_fig1.png',p1,width = 6,height = 5)



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
ggsave("results//run2_fig2_rank_geneEffect1.png",ps,width = 8,height = 5)

tmp %>% group_by(cell_line_name) %>%
  filter(geneType=="other") %>% 
  arrange(desc(Effect)) %>% 
  mutate(rank_ID=1:n()) %>% 
  filter(Gene %in% paper_report) -> top5_essential_genes
p10_C=tmp %>% group_by(cell_line_name) %>%  
  filter(geneType=="other") %>% 
  arrange(desc(Effect)) %>% 
  mutate(rank_ID=1:n()) %>% 
  ggplot(aes(x=rank_ID,y=Effect))+
  geom_point(aes(colour=geneType))+
  geom_text_repel(aes(x=rank_ID,y= Effect,label=Gene),data=top5_essential_genes)+
  geom_point(aes(x=rank_ID,y= Effect),colour='grey',data=top5_essential_genes)+
  facet_wrap(~cell_line_name)+  
  labs(title="Rank Effect")+
  theme_bw(
    base_size = 15
  )+
  theme(
    legend.position = 'bottom'
  )+
  guides(colour = guide_legend(nrow = 2))

tmp %>% inner_join(.,rra_gene_summary) %>% 
  group_by(cell_line_name) %>%  
  filter(geneType=="other") %>% 
  arrange(desc(logfc)) %>% 
  mutate(rank_ID=1:n()) %>%
  filter(Gene %in% paper_report) -> top5_essential_genes
p10_D=tmp %>% inner_join(.,rra_gene_summary) %>% 
  group_by(cell_line_name) %>%  
  filter(geneType=="other") %>% 
  arrange(desc(logfc)) %>% 
  mutate(rank_ID=1:n()) %>% 
  ggplot(aes(x=rank_ID,y=logfc))+
  geom_point(aes(colour=geneType))+
  geom_text_repel(aes(x=rank_ID,y= logfc,label=Gene),data=top5_essential_genes)+
  geom_point(aes(x=rank_ID,y= logfc),colour='grey',data=top5_essential_genes)+
  facet_wrap(~cell_line_name)+  
  labs(title="RRA logfc")+
  theme_bw(
    base_size = 15
  )+
  theme(
    legend.position = 'bottom'
  )+
  guides(colour = guide_legend(nrow = 2))
ps=ggpubr::ggarrange(p10_C,p10_D,nrow = 1)
ggsave("results//run2_fig2_rank_geneEffect2.png",ps,width = 8,height = 5)


#RRA的分析结果
rra_gene_summary=readxl::read_xlsx('RRA/D18_CompB_vs_D18_DMSO.RRA.gene_summary.xlsx') %>% 
  select(c(1,7)) %>% setnames(c('Gene','logfc'))
p1=tmp %>% inner_join(.,rra_gene_summary) %>% 
  ggplot(aes(x=logfc,y=Effect))+
  geom_point()+
  theme_bw(
    base_size = 15
  )+
  labs(x='RRA LogFC',y='Chronos Effect')
tmp %>% inner_join(.,rra_gene_summary) %>% 
  ggplot(aes(x=logfc,y=Effect))+
  geom_point(aes(colour=geneType))+
  facet_wrap(~geneType)+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'none')+
  labs(x='RRA LogFC',y='Chronos Effect')
  
tmp %>% inner_join(.,rra_gene_summary) %>%
  filter(Gene %in% paper_report) -> top5_essential_genes
candidate1=tmp %>% inner_join(.,rra_gene_summary) %>%
  filter(logfc< -1 & Effect< -1.5)
p2=tmp %>% inner_join(.,rra_gene_summary) %>% 
  ggplot(aes(x=logfc,y=Effect))+
  geom_point(aes(colour=geneType))+
  geom_text_repel(aes(x=logfc,y= Effect,label=Gene,size=3.3),data=top5_essential_genes)+
  geom_point(aes(x=logfc,y= Effect,size=3),data=top5_essential_genes)+
  geom_text_repel(aes(x=logfc,y= Effect,label=Gene,size=3.3),data=candidate1)+
  geom_point(aes(x=logfc,y= Effect,size=4),shape=21,data=candidate1)+
  facet_wrap(~geneType)+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'none')+
  labs(x='RRA LogFC',y='Chronos Effect')

#输出Chronos可能的essential genes
tmp %>% inner_join(.,rra_gene_summary) %>%
  filter(Effect< -1.5) %>% write.csv(.,'results/rra_chronos.csv')
tmp %>% inner_join(.,rra_gene_summary) %>%saveRDS(.,'results/rra_chronos.rds')



ps=ggpubr::ggarrange(p1,p2,nrow = 1)
ggsave("results//run2_fig3_rra_chronos.png",ps,width = 12,height = 6)

tmpx=tmp %>% inner_join(.,rra_gene_summary) %>% 
  filter(geneType=='other') %>% 
  filter(Effect< -1.5 & abs(logfc)<0.5)
p1=tmp %>% inner_join(.,rra_gene_summary) %>% 
  filter(geneType=='other') %>% 
  ggplot(aes(x=logfc,y=Effect))+
  geom_point(colour='purple')+
  geom_text_repel(aes(x=logfc,y= Effect,label=Gene),data=tmpx)+
  geom_rect(aes(xmin = -0.5, xmax = 0.5, ymin = -2.6, ymax = -1.5), color='black', 
            linetype='dashed', alpha=0) +
  theme_bw(
    base_size = 15
  )
ggsave("results//run2_fig3_rra_chronos_s1.png",p1,width = 6,height = 6)
# install.packages("remotes")
#remotes::install_github("hughjonesd/ggmagnify")
# px1=tmp %>% inner_join(.,rra_gene_summary) %>% 
#   filter(geneType=='other') %>% 
#   filter(Effect< -1.5 & abs(logfc)<0.5) %>% 
#   ggplot(aes(x=logfc,y=Effect))+
#   geom_point(colour='purple')+
#   geom_text_repel(aes(x=logfc,y= Effect,label=Gene))+
#   theme_pubr( base_size = 8,border = TRUE)


#RRA 和 MLE结果 rank 可视化，与fig1B类似
rra=readxl::read_excel('./RRA/D18_CompB_vs_D18_DMSO.RRA.gene_summary.xlsx') %>% 
  select(!contains('pos')) %>% 
  setnames(c('id','num','pvalue','fdr','rank','goodsgRNA','lfc','lab'))
paper_report=c("RRAGA",'MLST8','WDR24','LAMTOR2',
               "RREB1","GFI1","AKT1")
rra_select=rra %>% filter(id %in% paper_report)
#文章是使用RRA实现的比较分析
rra %>% mutate(
  log10_score=-log10(pvalue)
) %>% 
  filter(lab=='Essential_Gene' | lab=='Others') %>% 
  ggplot(
  aes(
  x=rank,
  y=log10_score
  )
)+geom_point()+
  geom_point(aes(x=rank,y=-log10(pvalue)),color='red',data=rra_select)+
  ggrepel::geom_text_repel(
    aes(x=rank,y=-log10(pvalue),
        label=id),
    data=rra_select
  )


mle=readxl::read_excel("MLE/feedback-mle/ET.MLE.gene_summary.withLabel.xlsx")
mle=fread("MLE/test.gene_summary.txt") %>% 
  select(c(1:5)) %>% 
  setnames(c('gene','sgRNA','beta','z','pvalue'))
mle %>% filter(gene %in% paper_report)


paper_report=c("RRAGA",'MLST8','WDR24','LAMTOR2',
               "RREB1","GFI1","AKT1")
tmp %>% filter(Gene %in% paper_report)
tmp %>% filter(geneType=='other') %>% filter(Effect< -0.5) 

#评估和RRA的相关性（effect vs log2FC），之前的结果表明，它两个相关性很高：D:\CRISPR\chronos_vs_mageck
