# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 21:39:27 2020

@author: wilke
"""
import numpy as np # For numerical fast numerical calculations
import matplotlib.pyplot as plt # For making plots
import pandas as pd # Deals with data
import seaborn as sns # Makes beautiful plots
sns.set_theme(style="darkgrid")


# configs = ["BabySFC", "PSI"]
# fields = ["Generated", "COMSOL"]
# starting= ["Corners", "Tim", "fine_tuned"]
# scalings = ["50","5","20"]

# for config in configs:
#     for field in fields:
#         for start in starting:
#             for scale in scalings:
#                 configuration = [config,field,start, scale]
#                 cond_list_corners = pd.DataFrame(np.array(pd.read_csv("Output"+configuration[0]+"/"+configuration[1]+"_Fields/cond/"+configuration[2]+"_50_"+configuration[3]+"_finalised_run.csv")), columns = ["cond"])
#                 plot1 = sns.relplot(x=cond_list_corners.index, y="cond", data=cond_list_corners,kind="line")
#                 plot1.set_axis_labels("Iterations [au]", "Condition Number [au]")
#                 plt.title("Condtion Number for the Descent: ["+configuration[0]+","+configuration[1]+", "+configuration[2]+", 50:"+configuration[3]+"]")
#                 plot1.savefig("WriteUp/"+configuration[1]+"_Plots/"+configuration[0]+"_"+configuration[2]+"_50_"+configuration[3])
                
configuration = ["PSI","Generated","Tim_10000", "5"]
# cond_list = pd.DataFrame(np.array(pd.read_csv("OutputBabySFC/COMSOL_Fields/airss_end_conds_50_50_finalised_run.csv")), columns = ["before","after"])
#sns.relplot(x=cond_list.index, y="before", data=cond_list);

cond_list_corners = pd.DataFrame(np.array(pd.read_csv("Output"+configuration[0]+"/"+configuration[1]+"_Fields/cond/"+configuration[2]+"_50_"+configuration[3]+"_finalised_run.csv")), columns = ["cond"])




plot1 = sns.relplot(x=cond_list_corners.index, y="cond", data=cond_list_corners,kind="line")
plot1.set_axis_labels("Iterations [au]", "Condition Number [au]")
plt.title("Condtion Number for the Descent: ["+configuration[0]+","+configuration[1]+", "+configuration[2]+", 50:"+configuration[3]+"]")
plot1.savefig("WriteUp/"+configuration[1]+"_Plots/"+configuration[0]+"_"+configuration[2]+"_50_"+configuration[3])



