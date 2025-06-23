# -*- coding: utf-8 -*-
"""
@author: pauldipjoti
"""

import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import pandas as pd

#%%
current_directory = os.getcwd()
data = pd.read_csv('mtlcad_3D_spheres_FF11_final.csv')
positions = data.iloc[:, 0:3].values
radii = data.iloc[:, 3].values
num_spheres = len(radii)
mean_radius = np.mean(radii);                          print(mean_radius) 
mean_radius2 = np.mean(radii ** 2);                    print(mean_radius2) 
mean_radius3 = np.mean(radii ** 3);                    print(mean_radius3)

sphere_total_volume = np.sum((4 / 3) * np.pi * (radii ** 3))
box_volume = 2 ** 3
FF = (sphere_total_volume / box_volume) * 100;         print(FF)  

