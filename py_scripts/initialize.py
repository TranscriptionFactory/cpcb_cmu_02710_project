import sys
sys.path.insert(0, './CSHMM-for-time-series-scRNA-Seq/')
import pandas as pd
import os
import subprocess

import scdiff_init

data_file = 'data/tcell_data_copy.E'

scdiff_init.run_scdiff_init(data_file)