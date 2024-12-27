

#计算回归曲线的，以及回归曲线中散点图的分布
chronos_gene_effect<-fread("Project_02/Achilles_run/gene_effects.csv")
geneTypes=readxl::read_xlsx("QC_crispr/example_data/4_types_of_genes.xlsx") %>% setnames(c("Gene","geneType"))

chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  inner_join(.,geneTypes) %>% count(geneType)

#同一化合物，不同时间点比较
chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  inner_join(.,geneTypes) %>% filter(geneType!="Essential_Gene") %>% 
  mutate(Days=gsub(".*(D[0-9]{2})\\-.*","\\1",cell_line_name)) %>% 
  mutate(Compound=gsub(".*\\-(\\w*)","\\1",cell_line_name)) %>%
  select(c(Days,Gene,Effect,Compound)) %>% 
  tidyr::pivot_wider(names_from = Days, values_from = c(Effect)) %>% 
  setnames(c("Gene","Compound","y","x")) -> regression_data
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
  geom_text_repel(aes(x=x,y=y,label=Gene),
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
  facet_wrap(~Compound)+
  labs(x="D21",y="D14")+
  theme_bw(
    base_size = 15
  )



#不同化合物，同一时间点比较
chronos_gene_effect %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect') %>% 
  inner_join(.,geneTypes) %>% filter(geneType!="Essential_Gene") %>% 
  mutate(Days=gsub(".*(D[0-9]{2})\\-.*","\\1",cell_line_name)) %>% 
  mutate(Compound=gsub(".*\\-(\\w*)","\\1",cell_line_name)) %>%
  select(c(Days,Gene,Effect,Compound)) %>% 
  tidyr::pivot_wider(names_from = Compound, values_from = c(Effect))%>% 
  setnames(c("Days","Gene","y","x")) -> regression_data
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
p2=regression_data %>% 
  ggplot(aes(x=x,y=y))+
  geom_point()+
  geom_point(aes(x=x,y=y),colour='grey',
             data=regression_data %>% 
               mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
               filter(abs(y_dist) < sigma_sd)
  )+
  geom_text_repel(aes(x=x,y=y,label=Gene),
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
  facet_wrap(~Days)+
  labs(x="NTC",y="KAT7")+
  theme_bw(
    base_size = 15
  )

ps=ggarrange(p1,p2,nrow = 1)
ggsave("fig6_cross_comparisions_withoutessenGenes.jpg",ps,width = 12,height = 4)


#输出outliers
regression_data %>% 
  mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
  filter(abs(y_dist) > sigma_sd)
# regression_data %>% mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )  %>% 
#   ggplot(aes(x=y_dist))+
#   geom_histogram()
# tan_=candidate_alpha=y/x
# cos_=1/sqrt(1+tan_^2)


#对mageck RRA的结果也尝试这样分析
#同一化合物，不同时间点比较
tmps_df %>% filter(geneType!="Essential_Gene") %>% 
  mutate(Days=gsub(".*(D[0-9]{2})\\-.*","\\1",groups)) %>% 
  mutate(Compound=gsub(".*\\-(\\w*)","\\1",groups)) %>%
  select(c(Days,gene,logfc,Compound)) %>% 
  tidyr::pivot_wider(names_from = Days, values_from = c(logfc)) %>%
  setnames(c("Gene","Compound","y","x")) -> regression_data
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
  geom_text_repel(aes(x=x,y=y,label=Gene),
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
  facet_wrap(~Compound)+
  labs(x="D21",y="D14")+
  theme_bw(
    base_size = 15
  )

#不同化合物，同一时间点比较
tmps_df %>% filter(geneType!="Essential_Gene") %>% 
  mutate(Days=gsub(".*(D[0-9]{2})\\-.*","\\1",groups)) %>% 
  mutate(Compound=gsub(".*\\-(\\w*)","\\1",groups)) %>%
  select(c(Days,gene,logfc,Compound)) %>% 
  tidyr::pivot_wider(names_from = Compound, values_from = c(logfc)) %>%
  setnames(c("Days","Gene","y","x")) -> regression_data
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
p2=regression_data %>% 
  ggplot(aes(x=x,y=y))+
  geom_point()+
  geom_point(aes(x=x,y=y),colour='grey',
             data=regression_data %>% 
               mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
               filter(abs(y_dist) < sigma_sd)
  )+
  geom_text_repel(aes(x=x,y=y,label=Gene),
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
  facet_wrap(~Days)+
  labs(x="NTC",y="KAT7")+
  theme_bw(
    base_size = 15
  )

ps=ggarrange(p1,p2,nrow = 1)
ggsave("fig6_cross_comparisions_withoutessenGenes_mageck.jpg",ps,width = 12,height = 4)
