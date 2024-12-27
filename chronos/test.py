

#基于chronos技术，对CRISPR数据进行分析，获取gene effect结果
import numpy as np
import pandas as pd
import chronos
import os
from matplotlib import pyplot as plt
import seaborn as sns


from matplotlib import rcParams
rcParams['axes.titlesize'] = 14
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False
rcParams['savefig.dpi'] = 200
rcParams['savefig.transparent'] = False
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = '11'
rcParams['figure.dpi'] = 200
rcParams["savefig.facecolor"] = (1, 1, 1.0, 0.2)

rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 7


#读入数据：
#AvanaSequenceMap: sgRNA与基因
#AvanaGuideMap: 实验设计
#AvanaReadcounts: sgRNA counts值
sequence_map = pd.read_csv("Data/SampleData/AvanaSequenceMap.csv")
guide_map = pd.read_csv("Data/SampleData/AvanaGuideMap.csv")
readcounts = chronos.read_hdf5("Data/SampleData/AvanaReadcounts.hdf5")


#获取essential和non-essential genes
common_essentials = pd.read_csv("Data/SampleData/AchillesCommonEssentialControls.csv")["Gene"]
nonessentials = pd.read_csv("Data/SampleData/AchillesNonessentialControls.csv")["Gene"]
#生成positive 和 negative 对照
positive_controls = guide_map.sgrna[guide_map.gene.isin(common_essentials)]
negative_controls = guide_map.sgrna[guide_map.gene.isin(nonessentials)]


#options：过滤outlier，通过log FC值
chronos.nan_outgrowths(readcounts=readcounts, guide_gene_map=guide_map,
                                   sequence_map=sequence_map)

#QC data
reportdir = "./Data/reports"
if not os.path.isdir(reportdir):
    os.mkdir(reportdir)
from chronos import reports
metrics = reports.qc_initial_data("Initial QC", readcounts, sequence_map,guide_map, 
        negative_controls, positive_controls,
        directory=reportdir
       )

#creating and training model
logdir = "./Data/logs"
if not os.path.isdir(logdir):
    os.mkdir(logdir)
model = chronos.Chronos(
    sequence_map={"avana": sequence_map},
    guide_gene_map={"avana": guide_map},
    readcounts={"avana": readcounts},
    negative_control_sgrnas={"avana": negative_controls},
    log_dir=logdir
)
#training
model.train(301, report_freq=50, burn_in_period=50, ge_only=0)

#save results and model
savedir = "Data/Achilles_run"
if not os.path.isdir(savedir):
    os.mkdir(savedir)

model.save(savedir, overwrite=True)
print("Saved files:\n\n" + '\n'.join(['\t' + s for s in os.listdir(savedir)
                if s.endswith("csv")
                or s.endswith("hdf5")
                or s.endswith("json")
                ]))
#load model
model_restored = chronos.load_saved_model(savedir)
gene_effects = model.gene_effect
#output to CSV
gene_effects.to_csv("gene_effects.csv")

#options: CNV矫正
cn = chronos.read_hdf5("Data/SampleData/OmicsCNGene.hdf5")

#个别基因可能没有CNV数据，需要剔除
try:
    corrected, shifts = chronos.alternate_CN(gene_effects, cn)
except ValueError as e:
    print(e)
    
for col in set(gene_effects.columns) - set(cn.columns):
    cn[col] = 1
corrected, shifts = chronos.alternate_CN(gene_effects, cn)
corrected.to_csv("gene_effects_corrected.csv")

#options: 联合多组学的QC
maf = pd.read_csv("Data/SampleData/OmicsSomaticMutations.csv")
cancer_relevant = maf[
  (
      maf.Driver | maf.LikelyDriver  
  ) & (
      maf.LikelyGoF
  )
]
cancer_relevant = cancer_relevant[~cancer_relevant.duplicated(subset=["ModelID", "Gene"])]
cancer_relevant['truecol'] = True
gof_matrix_base = pd.pivot(cancer_relevant, index="ModelID", columns="Gene", values="truecol")

expression_addictions = pd.read_csv("Data/SampleData/RNAiExpressionAddictions.csv")['Gene']
addiction_expressions = chronos.read_hdf5("Data/SampleData/OmicsExpressionProteinCodingGenesTPMLogp1.hdf5")[
    expression_addictions
]

metrics = reports.dataset_qc_report("ChronosAvana", savedir, 
                          common_essentials, nonessentials,
                          gof_matrix_base, addiction_expressions,
                          cn, directory="Data/reports",
                          gene_effect_file="gene_effect_corrected.hdf5"
                         )






