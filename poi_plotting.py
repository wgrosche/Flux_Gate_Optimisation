# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 11:06:33 2020

@author: wilke
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

configuration = ["BabySFC","Generated","airss_end_conds", "5"]



cond_list = pd.DataFrame(np.array(pd.read_csv("Output"+configuration[0]+"/"+configuration[1]+"_Fields/num_sens_vals/"+configuration[2]+"_50_"+configuration[3]+"_finalised_run.csv")), columns = ["after","before"])

sns.relplot(x=cond_list.index, y="before", data=cond_list);