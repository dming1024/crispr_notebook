
#sgRNA count qc
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)


count_summary<-fread("QC_crispr/G1/group_G1.countsummary.txt")
p1=count_summary %>% 
  mutate(Unmapped=Reads-Mapped) %>% 
  select(c(Label,Unmapped,Mapped)) %>%
  tidyr::pivot_longer(c(Mapped,Unmapped),names_to = 'Group',values_to = 'Reads') %>%
  group_by(Label) %>% 
  mutate(Rate=round(Reads/sum(Reads),3)) %>% 
  ggplot(aes(x=Label,y=Reads))+
  geom_col(aes(fill=Group),position = position_stack(reverse = T),width = 0.8)+
  geom_text(aes(x=Label,y=Reads,label=ifelse(Group=="Mapped",Rate,'')),size=4) +
  labs(x="")+
  scale_fill_manual(values=c("Mapped"="#1874CD","Unmapped"="#7EC0EE"))+
  theme_bw(
    base_size = 15
  )+
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    legend.title=element_blank(),
    legend.position =  c(0.1,0.95),
    legend.background = element_blank()
  )
ggsave('fig1.png',p1,width = 10,height = 5)
  
p2=count_summary %>% 
  ggplot(aes(x=Label,y=GiniIndex))+
  geom_bar(stat = "identity",fill="#1874CD")+
  labs(x="")+
  theme_bw(
    base_size = 15
  )+
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
  )
ggsave('fig2.png',p2,width = 10,height = 5)


p3=count_summary %>% 
  ggplot(aes(x=Label,y=Zerocounts))+
  geom_bar(stat = "identity",fill="#1874CD")+
  labs(x="")+
  theme_bw(
    base_size = 15
  )+
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
  )
ggsave('fig3.png',p3,width = 10,height = 5)


sgRNA_count=fread("QC_crispr/G1/group_G1.count.txt")
p4=sgRNA_count %>% select(-c(1,2)) %>% 
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
ggsave('fig4.png',p4,width = 10,height = 5)


p5=sgRNA_count %>% select(-c(1,2)) %>% 
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
ggsave('fig5.png',p5,width = 10,height = 5)

p6=sgRNA_count %>% select(-c(1,2)) %>% 
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
  xlim(c(1,20000))

ggsave('fig6.png',p6,width = 10,height = 5)

gene_type=fread("QC_crispr/example_data/human.CRISPR_gene_type.depmap22Q2.csv")
experiments=c(
  "A375-D14-KAT7_vs_A375-Base-KAT7",
  "A375-D21-KAT7_vs_A375-Base-KAT7"
)
getRRA_gene_summary<-function(comparision){
  path_to_gene_summary=sprintf("QC_crispr/example_data/RRA/%s/%s.gene_summary.txt",comparision,comparision)
  rra_results=ReadRRA(path_to_gene_summary) %>% 
    select(c(1,2)) %>% 
    setnames(c("Gene","LFC")) %>% 
    left_join(.,gene_type) %>% 
    mutate(Classification=ifelse(is.na(label),
                                 ifelse(grepl("CTRL",Gene),
                                        'Non-target','other'),
                                 label)) %>% 
    mutate(Group=comparision) %>% 
    select(-label)
  return(rra_results)
}

tmps=lapply(
  experiments,
  FUN=function(x){
    tmp=getRRA_gene_summary(x)
  }
)
tmps_df=rbindlist(tmps)

p7=tmps_df %>% ggplot(aes(x=Classification,y=LFC,fill=Group))+
  geom_boxplot()+
  labs(x="Gene",y='Log2 Fold Change')+
  scale_fill_manual(values=c("#1874CD","#7EC0EE"))+
  theme_bw(
    base_size = 15
  )+
  theme(
    legend.title=element_blank(),
    legend.background = element_blank()
  )
ggsave('fig7.png',p7,width = 10,height = 5)
