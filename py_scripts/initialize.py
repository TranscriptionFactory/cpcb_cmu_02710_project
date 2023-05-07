import sys
sys.path.insert(0, './CSHMM-for-time-series-scRNA-Seq/')
import pandas as pd
import os
import subprocess

import scdiff_init

data_file = 'data/tcell_data2_kc6.E'

tf_dna_path = '/net/dali/home/mscbio/aar126/genomics/jubilant-barnacle/scdiff/tf_dna/Mouse_TF_targets.txt'

#initial clusters - number of time points? 
kc = 6

scdiff_init.run_scdiff_init(data_file, tf_dna_path, kc)