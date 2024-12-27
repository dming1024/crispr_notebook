
#KEGG GO GSEA 等功能分析

library(stringr)
library(stringi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)

degs=fread("results/dependency_genes.csv")
#设置低一些的cut-off，选择更多dependency genes
sigma_sd = regression_data %>% mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) ) %>% pull(y_dist) %>% sd()*6
regression_data %>% 
  mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
  filter(abs(y_dist) > sigma_sd) %>% arrange(desc(abs(y_dist))) %>% pull(Gene) -> temp_deg

#功能分析
readRDS('results/rra_chronos.rds') %>% 
  filter(Effect< -0.5 & logfc< -0.5) %>% pull(Gene) -> temp_deg
readRDS('results/rra_chronos.rds') %>% 
  filter(Effect< -1) %>% pull(Gene) -> temp_deg

#temp_deg=degs$Gene
ego1 <- enrichGO(gene         = temp_deg,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)
p1=dotplot(ego1,showCategory=15)+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35))+ggtitle("Cell Component")

ego2 <- enrichGO(gene         = temp_deg,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)
p2=dotplot(ego2,showCategory=15)+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35))+ggtitle("Biological Process")

ego3 <- enrichGO(gene         = temp_deg,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)
p3=dotplot(ego3,showCategory=15)+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35))+ggtitle("Molecular Function")

library(KEGG.db)
library(easyConvert)
a2=easyConvert::easyConvert(
  species = "HUMAN",
  queryList = temp_deg,
  queryType = "SYMBOL"
)
ego4<-enrichKEGG(gene = a2$ENTREZID,
                 use_internal_data=T,
                 pvalueCutoff = 1,
                 qvalueCutoff  = 1) 
p4=dotplot(ego4,showCategory=10)+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))+ggtitle("KEGG Pathway")
kegg_df <- setReadable(ego4, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


ps=ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)
ggsave("results/kegg_GO.png", ps,
       width=15,height = 12)


#GSVA 或 GSEA分析
#加入mTOR抑制剂后，对各组样本进行表达检测
readRDS('results/rra_chronos.rds') %>% 
  filter(Effect< -1.5) %>% pull(Gene) -> temp_deg
gene_fpkm<-fread("GSE126486_TS84_RNA_Seq_FPKM.csv/GSE126486_TS84_RNA_Seq_FPKM.csv") %>% 
  filter(symbol %in% temp_deg)


tmp=gene_fpkm %>% dplyr::select(-c(1,2,3)) %>% as.data.frame()
rownames(tmp)=gene_fpkm$symbol
pheatmap(
  tmp[rowSums(tmp)>0.1,],
  scale = 'row',
  cutree_rows  = 2,
  filename = 'results/fig9.png'
)
dev.off()
p10=tmp[c('LYST','SNX10'),] %>% mutate(Gene=rownames(.)) %>% 
  pivot_longer(contains("_"),names_to = 'Samples',values_to = 'FPKM') %>% 
  mutate(Group=gsub("(.*)\\_[0-9]","\\1",Samples)) %>% 
  ggplot(aes(x=Group,y=FPKM))+
  geom_boxplot(aes(fill=Group))+
  theme_bw(
    base_size = 15
  )+
  theme(legend.position = 'bottom',
        #axis.text.x = element_text(angle = 60,hjust = 0.9,vjust = 0.9)
        axis.text.x = element_blank()
        )+
  facet_wrap(~Gene,scales = 'free')
ggsave('results/fig10.png',p10,width = 6,height = 5)
