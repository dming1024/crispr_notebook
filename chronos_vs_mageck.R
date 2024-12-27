
#比较mageck rra  和 chronos的分析结果
#mageck rra的分析结果在这里，一共四组：
#D:\CRISPR\QC_crispr\example_data\RRA

#chronos的分析结果在这里，也是四组：
#D:\CRISPR\Project_02\Achilles_run

library(data.table)
library(dplyr)
library(ggplot2)


chronos_gene_effect<-fread("Project_02/Achilles_run/gene_effects.csv")
#wide2long
p1=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  ggplot(aes(Effect,colour=cell_line_name))+
  stat_density(geom='line',position = 'identity')+
  labs(title="Chronos")+
  theme_bw(
    base_size = 15
  )+
  theme(
    legend.position =  c(0.2,0.8),
    legend.background = element_blank()
  )

geneTypes=readxl::read_xlsx("QC_crispr/example_data/4_types_of_genes.xlsx") %>% setnames(c("Gene","geneType"))
p3=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  inner_join(.,geneTypes) %>% 
  ggplot(aes(Effect,colour=geneType))+
  stat_density(geom='line',position = 'identity')+
  facet_wrap(~cell_line_name)+
  labs(title="Chronos")+
  theme_bw(
    base_size = 15
  )

p5=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  inner_join(.,geneTypes) %>% 
  ggplot(aes(x=geneType,y=Effect,fill=cell_line_name))+
  geom_boxplot()+
  theme_bw(
    base_size = 15
  )+
  theme(
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1)
  )
  

#计算下SSMD，每个样本一个SSMD（essential vs nonessential）
chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  inner_join(.,geneTypes) %>%
  group_by(cell_line_name,geneType) %>% 
  summarise(
    u=mean(Effect),
    sd=sd(Effect)
  ) %>% 
  filter(grepl("Gene",geneType)) %>% 
  tidyr::pivot_wider(names_from = geneType, values_from = c(u,sd)) %>% 
  mutate(SSMD=(u_Essential_Gene-u_Nonessential_Gene)/sqrt(sd_Essential_Gene+sd_Nonessential_Gene))

#mageck rra的分析结果
#D:\CRISPR\QC_crispr\example_data\RRA\A375-D14-KAT7_vs_A375-Base-KAT7
tmps=lapply(
  list.files("QC_crispr/example_data/RRA/"),
  FUN = function(x){
    tmp=readxl::read_xlsx(paste0("QC_crispr/example_data/RRA/",x,"/",x,".RRA.gene_summary.xlsx")) %>%
      mutate(group=gsub("(.*)\\_vs\\_.*","\\1",x)) %>% 
      select(c(1,7,group,label)) %>% 
      setnames(c("gene",'logfc','groups','geneType'))
  }
)
tmps_df=rbindlist(tmps)
p2=tmps_df %>% 
  ggplot(aes(logfc,colour=groups))+
  stat_density(geom='line',position = 'identity')+
  labs(title="Mageck RRA")+
  theme_bw(
    base_size = 15
  )+
  theme(
    legend.position =  c(0.2,0.8),
    legend.background = element_blank()
  )

p4=tmps_df %>% 
  ggplot(aes(logfc,colour=geneType))+
  stat_density(geom='line',position = 'identity')+
  facet_wrap(~groups)+
  labs(title="Mageck RRA")+
  theme_bw(
    base_size = 15
  )

p6=tmps_df %>% ggplot(aes(x=geneType,y=logfc,fill=groups))+
  geom_boxplot()+
  theme_bw(
    base_size = 15
  )+
  theme(
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1)
  )

#计算SSMD
ssmd1=tmps_df %>% group_by(groups,geneType) %>% 
  summarise(
    u=mean(logfc),
    sd=sd(logfc)
  ) %>% 
  filter(grepl("Gene",geneType)) %>% 
  tidyr::pivot_wider(names_from = geneType, values_from = c(u,sd)) %>% 
  mutate(SSMD=(u_Essential_Gene-u_Nonessential_Gene)/sqrt(sd_Essential_Gene+sd_Nonessential_Gene)) %>% 
  select(c(1,6)) %>% mutate(Groups=c("Mageck"))%>% setnames(c("cell_line_name","SSMD",'Groups')) 

NNMD1=tmps_df %>% group_by(groups,geneType) %>% 
  summarise(
    u=median(logfc),
    mad=mad(logfc)
  ) %>% 
  filter(grepl("Gene",geneType)) %>% 
  tidyr::pivot_wider(names_from = geneType, values_from = c(u,mad)) %>% 
  mutate(NNMD=(u_Essential_Gene-u_Nonessential_Gene)/mad_Nonessential_Gene) %>% 
  select(c(1,6)) %>% mutate(Groups=c("Mageck"))%>% setnames(c("cell_line_name","NNMD",'Groups')) 


ssmd2=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  inner_join(.,geneTypes) %>%
  group_by(cell_line_name,geneType) %>% 
  summarise(
    u=mean(Effect),
    sd=sd(Effect)
  ) %>% 
  filter(grepl("Gene",geneType)) %>% 
  tidyr::pivot_wider(names_from = geneType, values_from = c(u,sd)) %>% 
  mutate(SSMD=(u_Essential_Gene-u_Nonessential_Gene)/sqrt(sd_Essential_Gene+sd_Nonessential_Gene)) %>% 
  mutate(Groups=c("Chronos")) %>% 
  select(c('cell_line_name','SSMD','Groups'))

NNMD2=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  inner_join(.,geneTypes) %>%
  group_by(cell_line_name,geneType) %>% 
  summarise(
    u=mean(Effect),
    mad=mad(Effect)
  ) %>% 
  filter(grepl("Gene",geneType)) %>% 
  tidyr::pivot_wider(names_from = geneType, values_from = c(u,mad)) %>% 
  mutate(NNMD=(u_Essential_Gene-u_Nonessential_Gene)/mad_Nonessential_Gene) %>% 
  mutate(Groups=c("Chronos")) %>% 
  select(c('cell_line_name','NNMD','Groups'))



p7=rbind(ssmd1,ssmd2) %>% as.data.frame() %>% ggplot(aes(x=cell_line_name,y=SSMD,fill=Groups))+
  geom_bar(stat = 'identity', position='dodge') +
  scale_fill_manual(values=c("Chronos"="#1874CD","Mageck"="#7EC0EE"))+
  labs(title="SSMD")+
  theme_bw(
    base_size = 15
  )+
  theme(
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1)
  )


p8=rbind(NNMD1,NNMD2) %>% as.data.frame() %>% ggplot(aes(x=cell_line_name,y=NNMD,fill=Groups))+
  geom_bar(stat = 'identity', position='dodge') +
  scale_fill_manual(values=c("Chronos"="#1874CD","Mageck"="#7EC0EE"))+
  labs(title="NNMD")+
  theme_bw(
    base_size = 15
  )+
  theme(
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1)
  )


#log2fc vs gene_effect 相关性系数
tmp=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  mutate(id=paste0(cell_line_name,"_",Gene))
p9=tmps_df %>% mutate(id=paste0(groups,"_",gene)) %>% 
  inner_join(.,tmp) %>% 
  ggplot(aes(x=logfc,y=Effect))+
  geom_point(aes(colour=geneType))+
  facet_wrap(~groups)+
  labs(title="Correlation")+
  theme_bw(
    base_size = 15
  )

library(ggrepel)
tmp=chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  inner_join(.,geneTypes) %>% 
  mutate(geneType=as.factor(geneType))
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
  filter(geneType!="Essential_Gene") %>% 
  arrange(desc(Effect)) %>% 
  mutate(rank_ID=1:n()) %>% 
  slice_min(.,Effect,n=5) -> top5_essential_genes
p10_B=tmp %>% group_by(cell_line_name) %>%  
  filter(geneType!="Essential_Gene") %>% 
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




#不同时间点的比较D14 VS D21，同一化合物不同时间点的一致性还是挺高的
p11=tmp %>% 
  mutate(Days=gsub(".*(D[0-9]{2})\\-.*","\\1",cell_line_name)) %>% 
  mutate(Compound=gsub(".*\\-(\\w*)","\\1",cell_line_name)) %>%
  select(c(Days,Gene,Effect,Compound)) %>% 
  tidyr::pivot_wider(names_from = Days, values_from = c(Effect)) %>% 
  ggplot(aes(x=D14,y=D21))+
  geom_point()+
  facet_wrap(~Compound)+
  theme_bw(
    base_size = 15
  )

#同一时间点，不同化合物处理后，会有不同的gene effect
p12=tmp %>% 
  mutate(Days=gsub(".*(D[0-9]{2})\\-.*","\\1",cell_line_name)) %>% 
  mutate(Compound=gsub(".*\\-(\\w*)","\\1",cell_line_name)) %>%
  select(c(Days,Gene,Effect,Compound)) %>% 
  tidyr::pivot_wider(names_from = Compound, values_from = c(Effect)) %>% 
  ggplot(aes(x=NTC,y=KAT7))+
  geom_point()+
  facet_wrap(~Days)+
  theme_bw(
    base_size = 15
  )






library(ggpubr)
ps=ggarrange(p1,p2,nrow = 1)
ggsave("fig1_overall_comparision.jpg",ps,width = 10,height = 5)

ps=ggarrange(p3,p4,nrow = 1,common.legend = T)
ggsave("fig2_geneTypes_comparision.jpg",ps,width = 10,height = 5)

ps=ggarrange(p5,p6,nrow = 1,common.legend = T)
ggsave("fig2_geneTypes_comparision_boxplot.jpg",ps,width = 10,height = 5)

ps=ggarrange(p7,p8,nrow = 1,common.legend = T)
ggsave("fig3_SSMD_NNMD.jpg",ps,width = 10,height = 6)
ggsave("fig4_correlation.jpg",p9,width = 7,height = 5)

ps=ggarrange(p10_A,p10_B,nrow = 1)
ggsave("fig5_rank_gene_effect.jpg",ps,width = 12,height = 8)


ps=ggarrange(p11,p12,nrow = 2)
ggsave("fig6_cross_comparisions.jpg",ps,width = 8,height = 8)
