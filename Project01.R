
#CRISPR screeing data

#在beta score水平，对结果进行可视化与分析
#MLE结果分析
library(MAGeCKFlute)
gdata=ReadBeta("Project_01/s252a09010-ZWR_mle.gene_summary.txt")
#对不同样本进行normalization
gdata_c=NormalizeBeta(gdata,samples = c('PD7_DMSO','PD7_EP'))
DensityView(gdata,samples = c('PD7_DMSO','PD7_EP'))
DensityView(gdata_c,samples = c('PD7_DMSO','PD7_EP'))

gdata_c$Control = rowMeans(gdata_c[,'PD7_DMSO', drop = FALSE])
gdata_c$Treatment = rowMeans(gdata_c[,'PD7_EP', drop = FALSE])

#correlation plot
p1 = ScatterView(gdata_c, 
                 "Control", 
                 "Treatment", 
                 groups = c("top", "bottom"), 
                 auto_cut_diag = TRUE, 
                 display_cut = TRUE
                 )
p1

#rank plot
gdata_c$Diff = gdata_c$Treatment - gdata_c$Control
gdata_c$Rank = rank(gdata_c$Diff)
rankdata = gdata_c$Treatment - gdata_c$Control
names(rankdata) = gdata_c$Gene
RankView(rankdata)

#square plot
p2=SquareView(gdata_c, ctrlname = "Control", treatname = "Treatment",
           label = "Gene",groups=c('topcenter','midright','bottomcenter'),
           x_cutoff = c(-1,0),
           y_cutoff = c(-1,0))#打上essential_genes标签
#会找到一些sgRNA：C2CD4A, CHGA等positive selection
#negative selection: SMG7,POLB

p1 = ScatterView(gdata_c, x = "Control", y = "Treatment", label = "Gene", 
                 model = "ninesquare", top = 5, display_cut = TRUE, force = 2)
print(p1)

#可以结合：rra/PD7_EP_DMSO.gene_summary文件使用
#获取logfc，pvalue
rra_gdata=ReadRRA("Project_01/rra/PD7_EP_DMSO.gene_summary.txt")
head(rra_gdata)

rra_sdata=ReadsgRRA("Project_01/rra/PD7_EP_DMSO.sgrna_summary.txt")
head(rra_sdata)
#volcano plot
rra_gdata$LogFDR = -log10(rra_gdata$FDR)
VolcanoView(rra_gdata,x='LFC',y='FDR',
            Label = 'Official',
            y_cutoff = 0.05)

SquareView(gdata_c,ctrlname = "Control", treatname = "Treatment",
           label = 'Gene',label.top = T,
           groups = c("topcenter"),
           groupnames = c("group1"))

#https://rdrr.io/bioc/MAGeCKFlute/src/R/SquareView.R

file3 = file.path(system.file("extdata", package = "MAGeCKFlute"),
                 "testdata/mle.gene_summary.txt")
dd = ReadBeta(file3)
groups=c('midleft')
groupnames=c('group1')
names(groupnames)=groups
SquareView(dd, ctrlname = "dmso", treatname = "plx", label = "Gene",
           groups = groups,
           groupnames = groupnames)
ddx=dd
names(ddx)=c('Symbol','Control','Treatment')
p1 = ScatterView(ddx, x = "Control", y = "Treatment", label = "Symbol",
                 groups = c("midleft", "topcenter", "midright", "bottomcenter"),
                 groupnames = c("Group1", "Group2", "Group3", "Group4"),
                 display_cut = TRUE,x_cut=c(-1,0),y_cut=c(-1,0))



#选择基因进行功能富集分析
Square9 = p2$data
idx=Square9$group=="topcenter"
geneList = Square9$Diff
names(geneList) = Square9$Gene[idx]
universe = Square9$Gene
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 6, bottom = 0)

#自定义square plot
library(ggplot2)
gdata_c %>% ggplot(aes(x=Control,y=Treatment))+
  geom_point(color="grey",size=0.5)+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))+
  xlab("Control beta score")+ylab("Treatment beta score")+
  geom_smooth(method = "lm", se=FALSE, color="grey") 

with(gdata_c,
     lm(Control~Treatment))
plot(gdata_c$Control,gdata_c$Treatment)
abline(b=0.894,a=-0.009768)
