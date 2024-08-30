import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import pandas as pd 

plt.rcParams.update({'font.size': 12})
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)

#thresholds 
kk=0.19396
kkk=0.29232
kkpi=0.26302
kketa=0.29760
kksigma=0.32556
kkpipi=0.33208
kkomega=0.34937
kketaprime=0.35806
kkpieta=0.36667
kkphi=0.37345
kkkk=0.38792

import os 
import subprocess 
import math 
import sys 

threebody_path_ubuntu = '/home/digonto/Codes/Practical_Lattice_v2/3body_quantization_v2/'
threebody_path_macos = '../../../'
macos_path2 = '/Users/digonto/GitHub/jackknife_codes/'
ubuntu_path2 = '/home/digonto/Codes/Practical_Lattice_v2/jackknife_codes/'

from sys import platform 
#print(platform)
if platform=="linux" or platform=="linux2":
    print("platform = ",platform)
    jackpath = ubuntu_path2
    threebody_path = threebody_path_ubuntu
elif platform=="darwin":
    print("platform = ",platform)
    jackpath = macos_path2
    threebody_path = threebody_path_macos

sys.path.insert(1, jackpath)

def random_color_generator():
    #color = np.random.randint(0, 256, size=3)
    color = np.random.uniform(0.0, 1.0, size=3)
    return tuple(color)
 
random_color = random_color_generator()
print(random_color)

def param_plotter():
    plt.rcParams.update({'font.size': 22,'legend.fontsize': 12})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    #plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    fig, ax = plt.subplots(2,1, figsize = (12,10))
    ax[0].set_title("$\mathcal{K}_{3,iso,0}$",fontsize=35)
    ax[1].set_title("$\mathcal{K}_{3,iso,1}$",fontsize=35)
    
    #ax.set_xlim([non_int_list[0].min(),non_int_list[0].max()])
    #ax[0].set_ylim([-1E6,1E6])
    #ax[1].set_ylim([0,2])
    ax[0].set_xlim([0.0,2.0])
    ax[1].set_xlim([0,2.0])

    #ax[0].set_ylabel("$\\frac{|E_{lat} - E_{QC}|}{\\sigma_{lat}}$")
    #ax[1].set_ylabel("$\\frac{|E_{lat} - E_{QC}|}{\\sigma_{lat}}$")
    ax[1].set_xlabel("$\chi^2/ndof$")

    filename = 'K3iso_vs_chisquare.dat'

    (chisquare,  K3iso0, K3iso0_err, K3iso1, K3iso1_err, fit_levels, comments) = np.genfromtxt(filename,dtype=(float,float,float,float,float,float,"O"),skip_header=1,unpack=True)

    print(comments)

    for i in range(0,len(chisquare),1):
        thisColor = random_color_generator()
        K3iso1_val = K3iso1[i]
        print(K3iso1_val)
        comment_string = comments[i]
        comment_string1 = str(comment_string,'utf-8')
        
        ax[0].errorbar(chisquare[i], K3iso0[i] , xerr=None, yerr=K3iso0_err[i],
                        marker='o',markerfacecolor="None", markersize=20, color=thisColor,
                        linestyle='none',markeredgewidth=2, capsize=20,zorder=4)
        
        ax[0].plot(chisquare[i], K3iso0[i] , 
                        marker='o',markerfacecolor="white", markersize=20, color=thisColor,
                        linestyle='none', markeredgewidth=2, zorder=5,label=comment_string1)
        
        if(K3iso1_val!=0.0):
            ax[1].errorbar(chisquare[i], K3iso1[i] , xerr=None, yerr=K3iso1_err[i],
                            marker='o',markerfacecolor="None", markersize=20, color=thisColor,
                            linestyle='none',markeredgewidth=2, capsize=20,zorder=4)
        
            ax[1].plot(chisquare[i], K3iso1[i] , 
                            marker='o',markerfacecolor="white", markersize=20, color=thisColor,
                            linestyle='none', markeredgewidth=2, zorder=5,label=comment_string1)

    
        ax[0].legend()
        ax[1].legend()
    fig.tight_layout() 
    outfile = 'K3iso_parameters_plot.pdf'
    outfile1 = 'K3iso_parameters_plot.png'   
    plt.savefig(outfile)
    plt.savefig(outfile1)
    plt.show()
    plt.close()

param_plotter()