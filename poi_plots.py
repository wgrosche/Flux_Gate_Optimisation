# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 12:25:38 2020

@author: wilke
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np # For numerical fast numerical calculations
import pandas as pd # Deals with data
import seaborn as sns # Makes beautiful plots

configuration = ["BabySFC","Generated","Corners_10000", "5"]
sns.set_theme(style="darkgrid")

poi_list = pd.DataFrame(np.array(pd.read_csv("Output"+configuration[0]+"/"+configuration[1]+"_Fields/poi/"+configuration[2]+"_50_"+configuration[3]+"_finalised_run.csv")), columns = ["0","1","2"])


viewing_angle = [0,2]
directions= ["x","y","z"]

fig2 = plt.figure()
ax2 = fig2.add_subplot(111, aspect='equal')

coils_dim = [5*0.262, 9*0.262, 5*0.262] 
coils_centre = [0,0,0]
msr_dim = [0.97, 0.97, 0.87] 
msr_centre = [0,0.315,0]
coils_corner = (np.subtract((coils_centre[viewing_angle[0]],coils_centre[viewing_angle[1]]),(coils_dim[viewing_angle[0]]/2,coils_dim[viewing_angle[1]]/2)))
msr_corner = (np.subtract((msr_centre[viewing_angle[0]],msr_centre[viewing_angle[1]]),(msr_dim[viewing_angle[0]]/2,msr_dim[viewing_angle[1]]/2)))


ax2.set_xlim([coils_centre[viewing_angle[0]]-coils_dim[viewing_angle[0]]/2-0.5,coils_centre[viewing_angle[0]]+coils_dim[viewing_angle[0]]/2+0.5])
ax2.set_ylim([coils_centre[viewing_angle[1]]-coils_dim[viewing_angle[1]]/2-0.5,coils_centre[viewing_angle[1]]+coils_dim[viewing_angle[1]]/2+0.5])


ax2.add_patch(patches.Rectangle(coils_corner,coils_dim[viewing_angle[0]],coils_dim[viewing_angle[1]],fill=False, color = 'blue')) 
ax2.add_patch(patches.Rectangle(msr_corner,msr_dim[viewing_angle[0]],msr_dim[viewing_angle[1]], fill=False, color = 'red'))

end_poi, = plt.plot(np.array(poi_list.values)[-9:-1,viewing_angle[0]], np.array(poi_list.values)[-9:-1,viewing_angle[1]], 'kx', label = "End")
start_poi, = plt.plot(np.array(poi_list.values)[0:8,viewing_angle[0]], np.array(poi_list.values)[0:8,viewing_angle[1]], 'rx', Label = "Start")
plt.title("POI for the Descent: ["+configuration[0]+","+configuration[1]+", "+configuration[2]+", 50:"+configuration[3]+"]")
plt.xlabel(directions[viewing_angle[0]]+" [m]")
plt.ylabel(directions[viewing_angle[1]]+" [m]")
#plt.legend(handles=[black_patch, red_patch])
first_legend = plt.legend(handles = [start_poi, end_poi], loc='best')

fig2.savefig("WriteUp/"+configuration[1]+"_Plots/"+configuration[0]+"_"+configuration[2]+"_50_"+configuration[3]+"_view_"+str(viewing_angle[1]), bbox_inches='tight')
