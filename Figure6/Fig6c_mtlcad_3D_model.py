# -*- coding: utf-8 -*-
import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
current_directory = os.getcwd()

#%%
os.chdir(r"mtlcad_3D_spheres-FF11-final")
data = pd.read_csv("mtlcad_3D_model_FF11_comsol_result.txt", sep=r"\s+", skiprows=5, header=None, names=["lambda", "Power (W)"]).apply(pd.to_numeric, errors="coerce")
lambda_F11, Power_F11 = data["lambda"].tolist(), data["Power (W)"].tolist()

data = pd.read_csv('mtlcad_3D_spheres_FF11_final.csv')
positions = data.iloc[:, 0:3].values
radii_F11 = data.iloc[:, 3].values
num_spheres_F11 = len(radii_F11);                          # print(num_spheres_F11)
mean_radius_F11 = np.mean(radii_F11);                      # print(mean_radius_F11*1e3) 
std_radii_F11 = np.std(radii_F11);                         print(std_radii_F11/mean_radius_F11) 
mean_radius2 = np.mean(radii_F11 ** 2);                    # print(mean_radius2) 
mean_radius3 = np.mean(radii_F11 ** 3);                    # print(mean_radius3)
mean_radius6 = np.mean(radii_F11 ** 6);
sphere_total_volume = np.sum((4 / 3) * np.pi * (radii_F11 ** 3))
box_volume = 2 ** 3                                        # COMSOL simulation box 
FF = (sphere_total_volume / box_volume);                   # print(FF*1e2)  
index = 1.77
P0_F11 = ((3e8/index)**2/6)*(1-(1/index**2))*(1-(1/index))*(4*np.pi*mean_radius2*1e-12/377)
Power_F11_P0 = Power_F11/(P0_F11*num_spheres_F11)  

lambda_F11_R0 = lambda_F11[25:55]
Power_F11_R0 = ((3e8/index)**2*1e-12/377) * (2*np.pi/np.array(lambda_F11_R0))**4 * ((index**2-1)/(index**2+2))**2 * (mean_radius6/(index**2*mean_radius3))
Power_F11_R0 = Power_F11_R0 * (1.8/100)                     
os.chdir(current_directory)
 
#%%
os.chdir(r"mtlcad_3D_spheres-FF37-final")
data = pd.read_csv("mtlcad_3D_model_FF37_comsol_result.txt", sep=r"\s+", skiprows=5, header=None, names=["lambda", "Power (W)"]).apply(pd.to_numeric, errors="coerce")
lambda_F37, Power_F37 = data["lambda"].tolist(), data["Power (W)"].tolist()

data = pd.read_csv('mtlcad_3D_spheres_FF37_final.csv')
positions = data.iloc[:, 0:3].values
radii_F37 = data.iloc[:, 3].values
num_spheres_F37 = len(radii_F37);                          # print(num_spheres_F37)
mean_radius_F37 = np.mean(radii_F37);                      # print(mean_radius_F37*1e3) 
std_radii_F37 = np.std(radii_F37);                         print(std_radii_F37/mean_radius_F37) 
mean_radius2 = np.mean(radii_F37 ** 2);                    # print(mean_radius2) 
mean_radius3 = np.mean(radii_F37 ** 3);                    # print(mean_radius3)
mean_radius6 = np.mean(radii_F37 ** 6);
sphere_total_volume = np.sum((4 / 3) * np.pi * (radii_F37 ** 3))
box_volume = 2 ** 3                                        # COMSOL simulation box 
FF = (sphere_total_volume / box_volume);                   # print(FF*1e2)  
index = 1.77
P0_F37 = ((3e8/index)**2/6)*(1-(1/index**2))*(1-(1/index))*(4*np.pi*mean_radius2*1e-12/377)
Power_F37_P0 = Power_F37/(P0_F37*num_spheres_F37)  

lambda_F37_R0 = lambda_F37[25:51]
Power_F37_R0 = ((3e8/index)**2*1e-12/377) * (2*np.pi/np.array(lambda_F37_R0))**4 * ((index**2-1)/(index**2+2))**2 * (mean_radius6/(index**2*mean_radius3))
Power_F37_R0 = Power_F37_R0 * (0.825/100)                     
os.chdir(current_directory)

#%%
os.chdir(r"mtlcad_3D_spheres-FF50-final")
data = pd.read_csv("mtlcad_3D_model_FF50_comsol_result.txt", sep=r"\s+", skiprows=5, header=None, names=["lambda", "Power (W)"]).apply(pd.to_numeric, errors="coerce")
lambda_F50, Power_F50 = data["lambda"].tolist(), data["Power (W)"].tolist()

data = pd.read_csv('mtlcad_3D_spheres_FF50_final.csv')
positions = data.iloc[:, 0:3].values
radii_F50 = data.iloc[:, 3].values
num_spheres_F50 = len(radii_F50);                          # print(num_spheres_F50)
mean_radius_F50 = np.mean(radii_F50);                      # print(mean_radius_F50*1e3)
std_radii_F50 = np.std(radii_F50);                         print(std_radii_F50/mean_radius_F50)  
mean_radius2 = np.mean(radii_F50 ** 2);                    # print(mean_radius2) 
mean_radius3 = np.mean(radii_F50 ** 3);                    # print(mean_radius3)
mean_radius6 = np.mean(radii_F50 ** 6);
sphere_total_volume = np.sum((4 / 3) * np.pi * (radii_F50 ** 3))
box_volume = 2 ** 3                                        # COMSOL simulation box 
FF = (sphere_total_volume / box_volume);                   # print(FF*1e2)  
index = 1.77
P0_F50 = ((3e8/index)**2/6)*(1-(1/index**2))*(1-(1/index))*(4*np.pi*mean_radius2*1e-12/377)
Power_F50_P0 = Power_F50/(P0_F50*num_spheres_F50)  

lambda_F50_R0 = lambda_F50[25:45]
Power_F50_R0 = ((3e8/index)**2*1e-12/377) * (2*np.pi/np.array(lambda_F50_R0))**4 * ((index**2-1)/(index**2+2))**2 * (mean_radius6/(index**2*mean_radius3))
Power_F50_R0 = Power_F50_R0 * (0.5/100)                     
os.chdir(current_directory)

#%%
plt.figure(1)
plt.plot(lambda_F11/mean_radius_F11,Power_F11_P0, "ob-", label="Config 1: $f$ = 11.6%", linewidth=1.75, markersize=4)
plt.plot(lambda_F37/mean_radius_F37,Power_F37_P0, "or-", label="Config 2: $f$ = 38.2%", linewidth=1.75, markersize=4)
plt.plot(lambda_F50/mean_radius_F50,Power_F50_P0, "ok-", label="Config 3: $f$ = 50.3%", linewidth=1.75, markersize=4)

plt.plot(lambda_F11_R0/mean_radius_F11,Power_F11_R0, "b-.", linewidth=1.15)
plt.plot(lambda_F37_R0/mean_radius_F37,Power_F37_R0, "r-.", linewidth=1.15)
plt.plot(lambda_F50_R0/mean_radius_F50,Power_F50_R0, "k-.", linewidth=1.15)

plt.xscale("log");      plt.yscale("log");      plt.xlim(0.1, 1e2);     plt.ylim(1e-3, 3)
plt.xlabel(r"$\mathrm{\lambda} / \langle R \rangle$", fontsize=16)
plt.ylabel(r"$\mathrm{P}/(N \times \mathbb{P}_{0})$", fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=14)  
plt.tick_params(axis='both', which='minor', labelsize=14)  
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend(fontsize=14, loc="best");            plt.tight_layout()
plt.show()
# svg_path = os.path.join(current_directory, 'Fig6c_mtlcad_3D_model.svg')
# plt.savefig(svg_path, format='svg', bbox_inches="tight", pad_inches=0)

