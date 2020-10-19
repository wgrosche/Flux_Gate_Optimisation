# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 18:04:10 2020

@author: wilke
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np # For numerical fast numerical calculations
import pandas as pd # Deals with data
import seaborn as sns # Makes beautiful plots

configuration = ["BabySFC","Generated","fine_tuned", "5"]
sns.set_theme(style="darkgrid")

poi_list = pd.DataFrame(np.array(pd.read_csv("Output"+configuration[0]+"/"+configuration[1]+"_Fields/Stability/"+configuration[2]+"_50_"+configuration[3]+"_finalised_run_all.csv")), 
                        columns = ["mag","x","y","z","-x", "-y", "-z"])

# poi_list = pd.DataFrame(np.array(pd.read_csv("Output"+configuration[0]+"/"+configuration[1]+"_Fields/Stability/"+configuration[2]+"_50_"+configuration[3]+"_finalised_run_2.csv")), 
                        # columns = ["mag"])

poi_fg_1 = np.zeros((int(len(poi_list)/8), 7))
for i in range(0, int(len(poi_list)/8)):
    poi_fg_1[i,:] = poi_list.values[8*i+1,:]
    

sns.relplot(poi_fg_1[:,0], poi_fg_1[:,3], kind = "line")
plt.title("Shift in the Condition Number as a Function of Flux Gate Perturbation")
plt.xlabel("Distance from Optimised Position [m]")
plt.ylabel("Condition Number Change [au]")


# sns.relplot(data = poi_list, kind = "line")
# plt.title("Shift in the Condition Number as a Function of Flux Gate Perturbation")
# plt.xlabel("Distance from Optimised Position [m]")
# plt.ylabel("Condition Number Change [au]")
