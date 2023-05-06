import sys
sys.path.insert(0, './CSHMM-for-time-series-scRNA-Seq/')
import pandas as pd
import os
import subprocess

import scdiff_init

data_file = 'data/tcell_data_copy.E'

tf_dna_path = '/net/dali/home/mscbio/aar126/genomics/jubilant-barnacle/scdiff/tf_dna/Mouse_TF_targets.txt'

#initial clusters - number of time points? 
kc = 7

scdiff_init.run_scdiff_init(data_file, tf_dna_path, kc)