import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy.linalg import block_diag
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PyPDF2 import PdfMerger

import os 
import subprocess 
import math 
import sys 

from timeit import default_timer as timer

from iminuit import minimize 

import random 

threebody_path_ubuntu = '/home/digonto/Codes/Practical_Lattice_v2/3body_quantization/'
threebody_path_macos = '/Users/digonto/GitHub/3body_quantization/'
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

import jackknife 
from jackknife import jackknife_resampling, jackknife_average, jackknife_error 
from lattice_data_covariance import covariance_between_states_szscl21_based

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

def plotter_with_1K3iso():
    plt.rcParams.update({'font.size': 22,'legend.fontsize': 12})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    #plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    K3iso0 = 2788251.35396196
    K3iso0_err = 3065027.045413032

    L_list = ["20","24"]

    irrep_list = ['000_A1m','100_A2','110_A2','111_A2','200_A2']

    mom_list = ['000','100','110','111','200']

    marker_list = ['o','*','s','v','^']

    num_states = 20 

    fig, ax = plt.subplots(figsize = (12,5))
    ax.set_ylabel("$\mathcal{K}_{3,iso}$")
    ax.set_xlabel("$a_t E_{cm}$")

    ax.set_xlim([0.26,0.45])
    ax.set_ylim([-1E7,1E7])
    
    for L in L_list:
        for i in range(0,len(irrep_list),1):
            for j in range(0,num_states,1):
                state_file = "K3iso_jackavg_" + L + "_P" + mom_list[i] + "_state_" + str(j) + ".dat"
                if(os.path.exists(state_file)):
                    (En, Ecm, norm, F3, K3iso, F2, G, K2inv, Hinv) = np.genfromtxt(state_file, unpack=True)
                    ax.plot(Ecm, K3iso, linestyle='solid', linewidth=1, color='darkred',zorder=4)
            
            central_file = "K3iso_jackavg_centralval_" + L + "_P" + mom_list[i] + ".dat"

            (central_Ecm,central_K3iso) = np.genfromtxt(central_file, unpack=True)

            if(L=="20"):
                markerface = "white"
            elif(L=="24"):
                markerface = "darkred"

            ax.plot(central_Ecm, central_K3iso , 
                    marker=marker_list[i],markerfacecolor=markerface, markersize=10, color="darkred",
                    linestyle='none', markeredgewidth=1, zorder=5,label="[" + mom_list[i] + "]")
    
    
    
    K3iso0_arr = np.full(1000,K3iso0)
    K3iso0_ini_arr = np.full(1000,K3iso0-K3iso0_err)
    K3iso0_fin_arr = np.full(1000,K3iso0+K3iso0_err)
    Ecm_for_K3iso = np.linspace(0.26,0.5,1000)
    ax.plot(Ecm_for_K3iso,K3iso0_arr,linestyle='solid',zorder=2,color='teal')
    ax.fill_between(Ecm_for_K3iso,K3iso0_ini_arr,K3iso0_fin_arr,color='teal',alpha=0.2)
    ax.legend()
    fig.tight_layout()
    outputfile='K3iso_with_latticedata_1.pdf'
    outputfile1='K3iso_with_latticedata_1.png'
    plt.savefig(outputfile)
    plt.savefig(outputfile1)
    plt.show()

def plotter_with_3K3iso():
    plt.rcParams.update({'font.size': 22,'legend.fontsize': 12})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    #plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    K3iso0 = 2788251.35396196
    K3iso0_err = 3065027.045413032

    K3iso0_oneparam = 387605.32299431
    K3iso0_err_oneparam = 273455.7377215678

    K3iso0_twoparam = 179534.12529502
    K3iso0_err_twoparam = 134555.88380875264
    K3iso1_twoparam = -14373406.842944479
    K3iso1_err_twoparam = 3243769.5369793973


    L_list = ["20","24"]

    irrep_list = ['000_A1m','100_A2','110_A2','111_A2','200_A2']

    mom_list = ['000','100','110','111','200']

    marker_list = ['o','*','s','v','^']

    num_states = 20 

    fig, ax = plt.subplots(figsize = (12,5))
    ax.set_ylabel("$\mathcal{K}_{3,iso}$")
    ax.set_xlabel("$a_t E_{cm}$")

    ax.set_xlim([0.26,0.45])
    ax.set_ylim([-1E7,1E7])
    
    for L in L_list:
        for i in range(0,len(irrep_list),1):
            for j in range(0,num_states,1):
                state_file = "K3iso_jackavg_" + L + "_P" + mom_list[i] + "_state_" + str(j) + ".dat"
                if(os.path.exists(state_file)):
                    (En, Ecm, norm, F3, K3iso, F2, G, K2inv, Hinv) = np.genfromtxt(state_file, unpack=True)
                    ax.plot(Ecm, K3iso, linestyle='solid', linewidth=0.5, color='darkred',zorder=4)
            
            central_file = "K3iso_jackavg_centralval_" + L + "_P" + mom_list[i] + ".dat"

            (central_Ecm,central_K3iso) = np.genfromtxt(central_file, unpack=True)

            if(L=="20"):
                markerface = "white"
            elif(L=="24"):
                markerface = "darkred"

            ax.plot(central_Ecm, central_K3iso , 
                    marker=marker_list[i],markerfacecolor=markerface, markersize=7, color="darkred",
                    linestyle='none', markeredgewidth=0.5, zorder=5,label= L + ", [" + mom_list[i] + "]")
    
    
    
    K3iso0_arr = np.full(1000,K3iso0)
    K3iso0_one_param_arr = np.full(1000,K3iso0_oneparam)
    K3iso0_ini_arr = np.full(1000,K3iso0-K3iso0_err)
    K3iso0_fin_arr = np.full(1000,K3iso0+K3iso0_err)

    K3iso0_one_param_ini_arr = np.full(1000,K3iso0_oneparam-K3iso0_err_oneparam)
    K3iso0_one_param_fin_arr = np.full(1000,K3iso0_err_oneparam+K3iso0_err_oneparam)

    K3iso_twoparam = np.zeros(1000)
    K3iso_ini_twoparam = np.zeros(1000)
    K3iso_fin_twoparam = np.zeros(1000)
   

    Ecm_for_K3iso = np.linspace(0.26,0.5,1000)

    for i in range(0,len(Ecm_for_K3iso),1):
        sval = Ecm_for_K3iso[i]**2
        K3iso_twoparam_val = K3iso0_twoparam + (sval - kkpi*kkpi)*K3iso1_twoparam
        K3iso_ini_twoparam_val = (K3iso0_twoparam - K3iso0_err_twoparam) + (sval - kkpi*kkpi)*(K3iso1_twoparam - K3iso1_err_twoparam)
        K3iso_fin_twoparam_val = (K3iso0_twoparam + K3iso0_err_twoparam) + (sval - kkpi*kkpi)*(K3iso1_twoparam + K3iso1_err_twoparam)
        K3iso_twoparam[i] = K3iso_twoparam_val
        K3iso_ini_twoparam[i] = K3iso_ini_twoparam_val
        K3iso_fin_twoparam[i] = K3iso_fin_twoparam_val

    ax.plot(Ecm_for_K3iso,K3iso0_arr,linestyle='solid',zorder=2,color='teal',label="000_A1m_state0")
    ax.fill_between(Ecm_for_K3iso,K3iso0_ini_arr,K3iso0_fin_arr,color='teal',alpha=0.2)
    
    ax.plot(Ecm_for_K3iso,K3iso0_one_param_arr,linestyle='solid',zorder=2,color='green',label="one_param_cutoff=0.29")
    ax.fill_between(Ecm_for_K3iso,K3iso0_one_param_ini_arr,K3iso0_one_param_fin_arr,color='green',alpha=0.2)
    
    ax.plot(Ecm_for_K3iso,K3iso_twoparam,linestyle='solid',zorder=2,color='blue',label="two_param_cutoff=0.31")
    ax.fill_between(Ecm_for_K3iso,K3iso_ini_twoparam, K3iso_fin_twoparam ,color='blue',alpha=0.2)
    
    
    ax.legend()
    fig.tight_layout()
    outputfile='K3iso_with_latticedata_3.pdf'
    outputfile1='K3iso_with_latticedata_3.png'
    plt.savefig(outputfile)
    plt.savefig(outputfile1)
    plt.show()


#plotter_with_1K3iso()

plotter_with_3K3iso()
