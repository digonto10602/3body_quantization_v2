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

threebody_path_ubuntu = '/home/digonto/Codes/Practical_Lattice_v2/3body_quantization/'
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


def config_maker( max_nsq ):
    nx = []
    ny = []
    nz = []

    count = 1
    for i in range(-max_nsq,max_nsq+1,1):
        for j in range(-max_nsq,max_nsq+1,1):
            for k in range(-max_nsq,max_nsq+1,1):
                config_sq = i*i + j*j + k*k 
                if(config_sq<=max_nsq):
                    print("config ",count,": ",i,j,k)
                    count = count + 1 
                    nx.append(i)
                    ny.append(j)
                    nz.append(k)
                else:
                    continue
        
    return nx,ny,nz 


def energy( atm, 
            xi, 
            Lbyas,
            nx,
            ny,
            nz  ):
    onebyxi = 1.0/xi 
    pi = np.pi
    twopibyLbyas = 2.0*pi/Lbyas 
    nsquare = nx*nx + ny*ny + nz*nz 

    return np.sqrt(atm*atm + onebyxi*onebyxi*twopibyLbyas*twopibyLbyas*nsquare)


def irrep_list_maker(energy_file):
    df = pd.read_csv(energy_file,delim_whitespace=True, header=None)
    dflist = df.values.tolist()
    irrep_list = []
    for i in range(len(dflist)):
        irrep = dflist[i][1]
    
        if len(irrep_list)==0:
            irrep_list.append(irrep)
        else:
            check_irrep_list_flag = 0
            for j in range(len(irrep_list)):
                temp_irrep = irrep_list[j]
            
                if(temp_irrep == irrep):
                    check_irrep_list_flag = 1
                    break 
            if(check_irrep_list_flag==0):
                irrep_list.append(irrep)

    return dflist, irrep_list 

def irrep_energy_list_maker(full_energy_list, fixed_irrep):
    Ecm_list = []
    Elat_list = []

    for i in range(len(full_energy_list)):
        if(full_energy_list[i][1]==fixed_irrep):
            Ecm_list.append(float(full_energy_list[i][2]))
            Elat_list.append(float(full_energy_list[i][3]))

    return Ecm_list, Elat_list 

def canonical_mom_maker(nx, ny, nz):
    temp_list = []
    temp_list.append(abs(nx))
    temp_list.append(abs(ny))
    temp_list.append(abs(nz))
    temp_list.sort(reverse=True)
    
    return temp_list



def non_int_spectrum_plotter(irrep):
    plt.rcParams.update({'font.size': 12})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    #plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    #spectrum_file1 = "spectrum/spectrum.L20." + str(irrep)
    #spectrum_file2 = "spectrum/spectrum.L24." + str(irrep)
    non_int_file = "non_int." + str(irrep)
    title = non_int_file.split(".")

    non_int_list = np.genfromtxt(non_int_file,unpack=True)
    #L1, E1, Eerr1 = np.genfromtxt(spectrum_file1,unpack=True)
    #L2, E2, Eerr2 = np.genfromtxt(spectrum_file2,unpack=True)

    #print(len(non_int_list))

    en_list_for_200_A2_cm = np.array([0.2814,0.2923,0.3249,0.3436,0.3468,0.3612,0.3660])
    L_array_np_for_cm = [20 for i in range(len(en_list_for_200_A2_cm))]
    en_list_for_200_A2_lat = np.array([0.3353,0.3445,0.3727,0.3890,0.3918,0.4047,0.4089])
    L_array_np_for_lat = [20 for i in range(len(en_list_for_200_A2_lat))]
    en_list_for_200_A2_24_cm = np.array([0.2765,0.2849,0.3112,0.3248,0.3274,0.3368,0.3381])
    L_array_np_for_24_cm = [24 for i in range(len(en_list_for_200_A2_24_cm))]
    fig, ax = plt.subplots(figsize = (5,10))
    ax.set_title(title[1],fontsize=22)
    #ax.set_xlim([non_int_list[0].min(),non_int_list[0].max()])
    ax.set_ylim([0.24,0.38])
    #ax.set_ylim([5,8])
    ax.set_xticks([20, 22, 24], ['20', ' ', '24'])
    #ax.set_yticks(np.arange(0.15, 0.31, step=0.09))
    ax.set_yticks([0.24,kkpi, 0.31, kkpipi,0.38], ['0.24','$KK\pi$', '0.31', '$KK\pi\pi$','0.38'])
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    #ax.axhline(y = kk, color = 'grey', linestyle = '--',linewidth=1)
    ax.axhline(y = kkpi, color = 'grey', linestyle = '--',linewidth=1)
    ax.axhline(y = kkpipi, color = 'grey', linestyle = '--',linewidth=1)
    #ax.axvline(x = 0.0, color = 'black', linestyle = '-',linewidth=1)
    for i in range(len(non_int_list)-1):
        xval = non_int_list[0]
        yval = non_int_list[i+1]
        ax.plot(xval,yval,linestyle='-',color='black',zorder=2)        
    
    #ax.errorbar(L1, E1 , xerr=None, yerr=Eerr1,
    #            marker='o',markerfacecolor="None", markersize=10, color="black",
    #            linestyle='none',capsize=5)
    #ax.errorbar(L2, E2 , xerr=None, yerr=Eerr2,
    #            marker='o',markerfacecolor="None", markersize=10, color="black",
    #            linestyle='none',capsize=5)
    ax.scatter(L_array_np_for_cm,en_list_for_200_A2_cm, marker="_", s=520,color="red",zorder=4)
    ax.scatter(L_array_np_for_24_cm,en_list_for_200_A2_24_cm, marker="_", s=520,color="red",zorder=4)
    #ax.scatter(L_array_np_for_lat,en_list_for_200_A2_lat, marker="_", s=520,color="green",zorder=4)

    fig.tight_layout()
    
    output_file = "nonint_spectrum_plot." + str(irrep) + ".pdf"
    
    plt.savefig(output_file)
    plt.draw()
    plt.close()


def spectrum_plotter_with_L20andL24():
    plt.rcParams.update({'font.size': 22,'legend.fontsize': 12})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    #plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    irrep_list = ['000_A1m','100_A2','110_A2','111_A2','200_A2']


    fig, ax = plt.subplots(1,5,figsize = (30,10))

    for i in range(0,len(irrep_list),1):

        
        irrep = irrep_list[i]

        spectrum_file1 = "V20_" + str(irrep)
        spectrum_file2 = "V24_" + str(irrep)
        #spectrum_file2 = "spectrum/spectrum.L24." + str(irrep)
        non_int_file = "non_int." + str(irrep)

        

        non_int_list = np.genfromtxt(non_int_file,unpack=True)
        L1, E1, Eerr1, EerrA = np.genfromtxt(spectrum_file1,unpack=True)
        L2, E2, Eerr2, EerrA2 = np.genfromtxt(spectrum_file2,unpack=True) 

        title = non_int_file.split(".")
        #ax[i].set_title(title[1],fontsize=35)
        #ax.set_xlim([non_int_list[0].min(),non_int_list[0].max()])
        ax[i].set_ylim([0.25,0.35])
        ax[i].set_xlim(18,26)
        ax[i].set_xticks([20, 24], ['20 ', '24'])
        #ax.set_yticks(np.arange(0.15, 0.31, step=0.09))
        if(i==0):
            ax[i].set_yticks([0.25,kkpi, 0.28, 0.31, kkpipi,0.34], ['0.25','$KK\pi$', '0.28', '0.31', '$KK\pi\pi$','0.34'])
        else:
            ax[i].set_yticks([])

        ax[i].tick_params(axis='both', which='major', labelsize=35)
        ax[i].tick_params(axis='both', which='minor', labelsize=35)
        #ax.axhline(y = kk, color = 'grey', linestyle = '--',linewidth=1)
        ax[i].axhline(y = kkpi, color = 'grey', linestyle = '--',linewidth=1)
        ax[i].axhline(y = kkpipi, color = 'grey', linestyle = '--',linewidth=1)


        for j in range(0,len(non_int_list)-1,1):
            xval = non_int_list[0]
            yval = non_int_list[j+1]
            ax[i].plot(xval,yval,linestyle='-',color='black',zorder=2)        
    
        ax[i].errorbar(L1, E1 , xerr=None, yerr=EerrA,
                    marker='o',markerfacecolor="None", markersize=20, color="black",
                    linestyle='none',markeredgewidth=2, capsize=20,zorder=4)
        ax[i].plot(L1, E1 , 
                    marker='o',markerfacecolor="white", markersize=20, color="black",
                    linestyle='none', markeredgewidth=2, zorder=5)
    
        ax[i].errorbar(L2, E2 , xerr=None, yerr=EerrA2,
                    marker='o',markerfacecolor="None", markersize=20, color="black",
                    linestyle='none',markeredgewidth=2, capsize=20,zorder=4)
        ax[i].plot(L2, E2 , 
                    marker='o',markerfacecolor="white", markersize=20, color="black",
                    linestyle='none', markeredgewidth=2, zorder=5)
    
        
        ax[i].axhline(y=0.35,linestyle='dashed',color='darkred')
        if(irrep=='000_A1m'):
            atEcm_cutoff=0.4092
            ax[i].axhline(y=atEcm_cutoff,linestyle='dashed',color='darkred')
        elif(irrep=='100_A2'):
            atEcm_cutoff= 0.3408
            ax[i].axhline(y=atEcm_cutoff,linestyle='dashed',color='darkred')
        elif(irrep=='110_A2'):
            atEcm_cutoff= 0.3548
            ax[i].axhline(y=atEcm_cutoff,linestyle='dashed',color='darkred')
        elif(irrep=='111_A2'):
            atEcm_cutoff= 0.3665
            ax[i].axhline(y=atEcm_cutoff,linestyle='dashed',color='darkred')
        elif(irrep=='200_A2'):
            atEcm_cutoff= 0.3368
            ax[i].axhline(y=atEcm_cutoff,linestyle='dashed',color='darkred')
        
        non_int_list = []
        #ax[i].legend() 
    fig.tight_layout()
    
    output_file = "KKpi_spectrum_plot_with_L20_and_L24.pdf"
    
    plt.savefig(output_file)
    #plt.draw()
    plt.close()

def spectrum_plotter_with_L20andL24_with_QCstates(cutoff_val):
    plt.rcParams.update({'font.size': 22,'legend.fontsize': 12})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    #plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    irrep_list = ['000_A1m','100_A2','110_A2','111_A2','200_A2']

    mom_list = ['000','100','110','111','200']

    fig, ax = plt.subplots(1,5,figsize = (30,10))

    for i in range(0,len(irrep_list),1):

        
        irrep = irrep_list[i]

        moms = mom_list[i]

        spec_drive = threebody_path + "/lattice_data/KKpi_interacting_spectrum/twoptvar_analysis/masses/"
        spectrum_file1 = spec_drive + "V20_" + str(irrep)
        spectrum_file2 = spec_drive + "V24_" + str(irrep)
        #spectrum_file2 = "spectrum/spectrum.L24." + str(irrep)
        non_int_file = spec_drive + "non_int." + str(irrep)

        QC_drive = threebody_path + "/test_files/QC_states_twoparams_fitting_L_20_24/cutoff_" + str(cutoff_val) + "/"
        QC_file1 = QC_drive + "QC_states_jackknifed_L_20_nP_" + moms + "_energycutoff_" + str(cutoff_val) + "_two_params.dat"
        QC_file2 = QC_drive + "QC_states_jackknifed_L_24_nP_" + moms + "_energycutoff_" + str(cutoff_val) + "_two_params.dat"

        non_int_list = np.genfromtxt(non_int_file,unpack=True)
        L1, E1, Eerr1, EerrA = np.genfromtxt(spectrum_file1,unpack=True)
        L2, E2, Eerr2, EerrA2 = np.genfromtxt(spectrum_file2,unpack=True) 

        QCL1, QCspec1, QCspecErr1, diff1, diff1err = np.genfromtxt(QC_file1,unpack=True)
        QCL2, QCspec2, QCspecErr2, diff2, diff2err = np.genfromtxt(QC_file2,unpack=True) 

        title = non_int_file.split(".")
        ax[i].set_title(irrep,fontsize=35)
        #ax.set_xlim([non_int_list[0].min(),non_int_list[0].max()])
        ax[i].set_ylim([0.25,0.35])
        ax[i].set_xlim(17,27)
        ax[i].set_xticks([20, 24], ['20 ', '24'])
        #ax.set_yticks(np.arange(0.15, 0.31, step=0.09))
        if(i==0):
            ax[i].set_yticks([0.25,kkpi, 0.28, 0.31, kkpipi,0.34], ['0.25','$KK\pi$', '0.28', '0.31', '$KK\pi\pi$','0.34'])
        else:
            ax[i].set_yticks([])

        ax[i].tick_params(axis='both', which='major', labelsize=35)
        ax[i].tick_params(axis='both', which='minor', labelsize=35)
        #ax.axhline(y = kk, color = 'grey', linestyle = '--',linewidth=1)
        ax[i].axhline(y = kkpi, color = 'grey', linestyle = '--',linewidth=1)
        ax[i].axhline(y = kkpipi, color = 'grey', linestyle = '--',linewidth=1)


        for j in range(0,len(non_int_list)-1,1):
            xval = non_int_list[0]
            yval = non_int_list[j+1]
            ax[i].plot(xval,yval,linestyle='-',color='black',zorder=2)        
    
        ax[i].errorbar(L1, E1 , xerr=None, yerr=EerrA,
                    marker='o',markerfacecolor="None", markersize=20, color="black",
                    linestyle='none',markeredgewidth=2, capsize=20,zorder=4)
        ax[i].plot(L1, E1 , 
                    marker='o',markerfacecolor="white", markersize=20, color="black",
                    linestyle='none', markeredgewidth=2, zorder=5)
    
        ax[i].errorbar(L2, E2 , xerr=None, yerr=EerrA2,
                    marker='o',markerfacecolor="None", markersize=20, color="black",
                    linestyle='none',markeredgewidth=2, capsize=20,zorder=4)
        ax[i].plot(L2, E2 , 
                    marker='o',markerfacecolor="white", markersize=20, color="black",
                    linestyle='none', markeredgewidth=2, zorder=5)
    
        
        ax[i].axhline(y=0.35,linestyle='dashed',color='darkred')
        if(irrep=='000_A1m'):
            atEcm_cutoff=0.4092
            ax[i].axhline(y=atEcm_cutoff,linestyle='dashed',color='darkred')
        elif(irrep=='100_A2'):
            atEcm_cutoff= 0.3408
            ax[i].axhline(y=atEcm_cutoff,linestyle='dashed',color='darkred')
        elif(irrep=='110_A2'):
            atEcm_cutoff= 0.3548
            ax[i].axhline(y=atEcm_cutoff,linestyle='dashed',color='darkred')
        elif(irrep=='111_A2'):
            atEcm_cutoff= 0.3665
            ax[i].axhline(y=atEcm_cutoff,linestyle='dashed',color='darkred')
        elif(irrep=='200_A2'):
            atEcm_cutoff= 0.3368
            ax[i].axhline(y=atEcm_cutoff,linestyle='dashed',color='darkred')
        


        ax[i].errorbar(QCL1 + 1.0, QCspec1 , xerr=None, yerr=QCspecErr1,
                    marker='o',markerfacecolor="None", markersize=10, color="red",
                    linestyle='none',markeredgewidth=1, capsize=20,zorder=4)
        ax[i].plot(QCL1 + 1.0, QCspec1 , 
                    marker='o',markerfacecolor="white", markersize=10, color="red",
                    linestyle='none', markeredgewidth=1, zorder=5)
    
        ax[i].errorbar(QCL2 + 1.0, QCspec2 , xerr=None, yerr=QCspecErr2,
                    marker='o',markerfacecolor="None", markersize=10, color="red",
                    linestyle='none',markeredgewidth=1, capsize=20,zorder=4)
        ax[i].plot(QCL2 + 1.0, QCspec2 , 
                    marker='o',markerfacecolor="white", markersize=10, color="red",
                    linestyle='none', markeredgewidth=1, zorder=5)
    
        
        #ax[i].fill_between(QCL1 + 1.0,QCspec1-QCspecErr1,QCspec1+QCspecErr1, alpha=0.5, color='red')
        #ax[i].fill_between(QCL2 + 1.0,QCspec2-QCspecErr2,QCspec2+QCspecErr2, alpha=0.5, color='red')
        
        non_int_list = []
        #ax[i].legend() 
    fig.tight_layout()
    
    output_file1 = "KKpi_spectrum_plot_with_L20_and_L24_with_QC_states.pdf"
    output_file2 = "KKpi_spectrum_plot_with_L20_and_L24_with_QC_states.png"
    
    plt.savefig(output_file1)
    plt.savefig(output_file2)
    #plt.draw()
    plt.close()

def diff_plotter(cutoff_val):
    plt.rcParams.update({'font.size': 22,'legend.fontsize': 12})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    #plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    irrep_list = ['000_A1m','100_A2','110_A2','111_A2','200_A2']

    mom_list = ['000','100','110','111','200']

    marker_list = ['o','v','<','s','D']

    fig, ax = plt.subplots(2,1, figsize = (12,10))

    ax[0].set_title("L=20",fontsize=35)
    ax[1].set_title("L=24",fontsize=35)
    
    #ax.set_xlim([non_int_list[0].min(),non_int_list[0].max()])
    ax[0].set_ylim([0,3])
    ax[1].set_ylim([0,2])
    ax[0].set_xlim([0.26,0.385])
    ax[1].set_xlim([0.26,0.385])

    ax[0].set_ylabel("$\\frac{|E_{lat} - E_{QC}|}{\\sigma_{lat}}$")
    ax[1].set_ylabel("$\\frac{|E_{lat} - E_{QC}|}{\\sigma_{lat}}$")
    ax[1].set_xlabel("$E_{cm}$")


    for i in range(0,len(irrep_list),1):
        print(i)
        spec_drive = threebody_path + "/lattice_data/KKpi_interacting_spectrum/twoptvar_analysis/masses/"
        spectrum_file1 = spec_drive + "V20_" + str(irrep_list[i])
        spectrum_file2 = spec_drive + "V24_" + str(irrep_list[i])

        L1, E1, Eerr1, EerrA = np.genfromtxt(spectrum_file1,unpack=True)
        L2, E2, Eerr2, EerrA2 = np.genfromtxt(spectrum_file2,unpack=True) 

        QC_drive = threebody_path + "/test_files/QC_states_twoparams_fitting_L_20_24/cutoff_" + str(cutoff_val) + "/"
        QC_file1 = QC_drive + "QC_states_jackknifed_L_20_nP_" + mom_list[i] + "_energycutoff_" + str(cutoff_val) + "_two_params.dat"
        QC_file2 = QC_drive + "QC_states_jackknifed_L_24_nP_" + mom_list[i] + "_energycutoff_" + str(cutoff_val) + "_two_params.dat"

        energycutoff_string = QC_file1.split('_') 
        for j in range(0,len(energycutoff_string),1):
            if(energycutoff_string[j]=='energycutoff'):
                cutoff = float(energycutoff_string[j+1])

        
        QCL1, QCspec1, QCspecErr1, diff1, diff1err = np.genfromtxt(QC_file1,unpack=True)
        QCL2, QCspec2, QCspecErr2, diff2, diff2err = np.genfromtxt(QC_file2,unpack=True) 

        fitted_spec1 = []
        fitted_spec1_diff = []
        fitted_spec1_diff_err = []
        nonfitted_spec1 = []
        nonfitted_spec1_diff = []
        nonfitted_spec1_diff_err = []
        fitted_spec2 = []
        fitted_spec2_diff = []
        fitted_spec2_diff_err = []
        nonfitted_spec2 = []
        nonfitted_spec2_diff = []
        nonfitted_spec2_diff_err = []                                 

        for j in range(0,len(E1),1):
            if(E1[j]<=cutoff):
                cutoff_ind1 = j
                print(j,E1[j],L1[j])
                #fitted_spec1.append(QCspec1[j])
                fitted_spec1.append(E1[j])
                fitted_spec1_diff.append(diff1[j])
                fitted_spec1_diff_err.append(diff1err[j])
            else:
                #nonfitted_spec1.append(QCspec1[j])
                nonfitted_spec1.append(E1[j])
                nonfitted_spec1_diff.append(diff1[j])
                nonfitted_spec1_diff_err.append(diff1err[j])
        for j in range(0,len(E2),1):
            if(E2[j]<=cutoff):
                cutoff_ind2 = j
                print(j,E2[j],L2[j])
                #fitted_spec2.append(QCspec2[j])
                fitted_spec2.append(E2[j])
                fitted_spec2_diff.append(diff2[j])
                fitted_spec2_diff_err.append(diff2err[j])
            else:
                #nonfitted_spec2.append(QCspec2[j])
                nonfitted_spec2.append(E2[j])
                nonfitted_spec2_diff.append(diff2[j])
                nonfitted_spec2_diff_err.append(diff2err[j]) 

        ax[0].errorbar(fitted_spec1, fitted_spec1_diff , xerr=None, yerr=fitted_spec1_diff_err,
                    marker=marker_list[i],markerfacecolor="None", markersize=20, color="red",
                    linestyle='none',markeredgewidth=1, capsize=20,zorder=4)
        
        ax[0].plot(fitted_spec1, fitted_spec1_diff , 
                    marker=marker_list[i],markerfacecolor="white", markersize=20, color="red",
                    linestyle='none', markeredgewidth=1, zorder=5,label=mom_list[i])

        ax[1].errorbar(fitted_spec2, fitted_spec2_diff, xerr=None, yerr=fitted_spec2_diff_err,
                    marker=marker_list[i],markerfacecolor="None", markersize=20, color="red",
                    linestyle='none',markeredgewidth=1, capsize=20,zorder=4)
        ax[1].plot(fitted_spec2, fitted_spec2_diff, 
                    marker=marker_list[i],markerfacecolor="white", markersize=20, color="red",
                    linestyle='none', markeredgewidth=1, zorder=5,label=mom_list[i])
        
        ax[0].errorbar(nonfitted_spec1, nonfitted_spec1_diff , xerr=None, yerr=nonfitted_spec1_diff_err,
                    marker=marker_list[i],markerfacecolor="None", markersize=20, color="grey",
                    linestyle='none',markeredgewidth=1, capsize=20,zorder=4)
        
        ax[0].plot(nonfitted_spec1, nonfitted_spec1_diff , 
                    marker=marker_list[i],markerfacecolor="white", markersize=20, color="grey",
                    linestyle='none', markeredgewidth=1, zorder=5)

        ax[1].errorbar(nonfitted_spec2, nonfitted_spec2_diff , xerr=None, yerr=nonfitted_spec2_diff_err,
                    marker=marker_list[i],markerfacecolor="None", markersize=20, color="grey",
                    linestyle='none',markeredgewidth=1, capsize=20,zorder=4)
        ax[1].plot(nonfitted_spec2, nonfitted_spec2_diff , 
                    marker=marker_list[i],markerfacecolor="white", markersize=20, color="grey",
                    linestyle='none', markeredgewidth=1, zorder=5)
        
        ax[0].axvline(x=cutoff,linestyle='dashed',zorder=1)
        ax[1].axvline(x=cutoff,linestyle='dashed',zorder=1)
        ax[0].legend()
        ax[1].legend()
    
    fig.tight_layout()
    outfile1 = "diff_plot.pdf"
    plt.savefig(outfile1)
    outfile2 = "diff_plot.png"
    plt.savefig(outfile2)
    #plt.close()
    plt.show()


#for L=20 and L=24
energyfile = "S2I2_energies"
#dflist, irrep_list = irrep_list_maker(energyfile)
#print(irrep_list)
#irrep = "200_A2"
#non_int_spectrum_plotter(irrep)
irrep_list = ["000_A1m","100_A2","110_A2","111_A2","200_A2"]
Pmom_list = ["000","100","110","111","200"]

for i in range(len(irrep_list)):
    irrep = irrep_list[i]
    #print(irrep)
    Pmom = Pmom_list[i]
    #spectrum_plotter_with_3_KDF(irrep, Pmom)

cutoff_val = 0.33
spectrum_plotter_with_L20andL24_with_QCstates(cutoff_val)

diff_plotter(cutoff_val)

def test():
    filename = "QC_states_jackknifed_L_20_nP_000_energycutoff_0.29_one_param.dat"

    newname = filename.split("_")

    for i in range(0,len(newname),1):
        if(newname[i]=='energycutoff'):
            val1 = newname[i+1]
            print(i+1, newname[i+1])
            num = 5
            print(num + float(val1))

#test()