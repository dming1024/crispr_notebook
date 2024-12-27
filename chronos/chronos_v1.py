

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

import yaml
with open('data.yaml', 'r') as file:
    configs = yaml.safe_load(file)

#读入数据：
sequence_map = pd.read_csv(configs['sequence_map'])
guide_map = pd.read_csv(configs['guide_map'])
readcounts = chronos.read_hdf5(configs['readcounts'])


#获取essential和non-essential genes
common_essentials = pd.read_csv(configs['common_essentials'])["Gene"]
nonessentials = pd.read_csv(configs['nonessentials'])["Gene"]
#生成positive 和 negative 对照
positive_controls = guide_map.sgrna[guide_map.gene.isin(common_essentials)]
negative_controls = guide_map.sgrna[guide_map.gene.isin(nonessentials)]


#options：过滤outlier，通过log FC值
chronos.nan_outgrowths(readcounts=readcounts, guide_gene_map=guide_map,
                                   sequence_map=sequence_map)

#QC data
reportdir = configs['output_qc']
if not os.path.isdir(reportdir):
    os.makedirs(reportdir)
from chronos import reports
metrics = reports.qc_initial_data("Initial QC", readcounts, sequence_map,guide_map, 
        negative_controls, positive_controls,
        directory=reportdir
       )

#creating and training model
logdir = configs['log_dir']
if not os.path.isdir(logdir):
    os.makedirs(logdir)
model = chronos.Chronos(
    sequence_map={configs['project_name']: sequence_map},
    guide_gene_map={configs['project_name']: guide_map},
    readcounts={configs['project_name']: readcounts},
    negative_control_sgrnas={configs['project_name']: negative_controls},
    log_dir=logdir
)
#training
model.train(301, report_freq=50, burn_in_period=50, ge_only=0)

#save results and model
savedir = configs['output_results']
if not os.path.isdir(savedir):
    os.makedirs(savedir)

model.save(savedir, overwrite=True)
print("Saved files:\n\n" + '\n'.join(['\t' + s for s in os.listdir(savedir)
                if s.endswith("csv")
                or s.endswith("hdf5")
                or s.endswith("json")
                ]))
#load model
#model_restored = chronos.load_saved_model(savedir)
gene_effects = model.gene_effect
#output to CSV
gene_effects.to_csv(savedir+"/gene_effects.csv")



