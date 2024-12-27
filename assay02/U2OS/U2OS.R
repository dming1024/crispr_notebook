
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
setwd("D:\\CRISPR\\20240511\\U2OS")

tmp=fread("s252a12t008_U2OS_count_all.count.txt") %>% 
  select(c(3:7)) %>% t() %>% as.data.frame()
colnames(tmp)=fread("s252a12t008_U2OS_count_all.count.txt") %>% pull(1)

tmp %>% write.csv(.,'counts.csv',quote = F)


#比较本地chrons分析结果与DepMap 中gene effect一致性
#ACH-000364, U2OS
Depmap_res=fread("D://AllDatabase//depmap/CRISPRGeneEffect.csv") %>% 
  filter(ModelID=='ACH-000364') %>% 
  select(-c(ModelID)) %>% 
  t(.) %>% 
  as.data.frame() %>% 
  setnames('Depmap') %>% 
  mutate(
    Gene=gsub("\\s.*\\(.*","",rownames(.))
  )

#local results
chronos_res=fread("run/gene_effects.csv") %>% 
  select(-c(cell_line_name)) %>% 
  t(.) %>% 
  as.data.frame() %>% 
  setnames('WuXi') %>% 
  mutate(
    Gene=rownames(.)
  ) %>% 
  inner_join(.,
             Depmap_res)

cor.test(
  chronos_res$WuXi,
  chronos_res$Depmap
)

#选几个共有的dependency 进行标记

chronos_res_labels=chronos_res %>% filter(WuXi< -1 & Depmap < -1) %>% 
  mutate(c1=WuXi+Depmap) %>% slice_min(.,c1,n=5)
  # mutate(
  #   c1=apply(.,
  #            1,
  #            FUN=function(x){
  #              (as.numeric(x[1])+as.numeric(x[3]))
  #            }),
  #   c2=apply(.,
  #            1,
  #            FUN=function(x){
  #              sd(c(as.numeric(x[1]),as.numeric(x[3])))
  #            })
  # )



p=chronos_res %>% 
  mutate(
    g1=ifelse(WuXi< -1 & Depmap < -1, 'Y','N')
  ) %>% filter(!is.na(g1)) %>% 
  ggplot(aes(x=WuXi,y=Depmap))+
  geom_point(aes(col=g1))+
  ggrepel::geom_text_repel(
    aes(x=WuXi,y=Depmap,label=Gene),data=chronos_res_labels
  )+annotate(geom="text",x=-3,y= 0,label='PCC= 0.71')+
  theme_bw(base_size = 20)+
  theme(legend.position = 'none')
ggsave('geneEffect_comparisions.png',p,width = 6,height = 6)

wuxi=chronos_res %>% filter(WuXi< -1) %>% pull(Gene)
depmap=chronos_res %>% filter(Depmap< -1) %>% pull(Gene)
length(intersect(wuxi,depmap))/length(wuxi)
