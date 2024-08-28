# Here we plot the F3 function for each boost frame and compare them with 
# free spectrum generated from threebody_non_int_spectrum() code 


#This plots F3_iso for different total momentum P
#This is for KKpi system

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import pandas as pd
import plotly.validators.scatter.marker


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

def energy_list_maker_fixedmom(full_energy_list, fixed_mom):
    Ecm_list = []
    Elat_list = []


    for i in range(len(full_energy_list)):
        if(full_energy_list[i][1][0]==fixed_mom[0] and full_energy_list[i][1][1]==fixed_mom[1] and full_energy_list[i][1][2]==fixed_mom[2]):
            Ecm_list.append(float(full_energy_list[i][2]))
            Elat_list.append(float(full_energy_list[i][3]))

    return Ecm_list, Elat_list 


def plotter():

    plt.rcParams.update({'font.size': 18})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    Lbyas = 24



    drive = "/home/digonto/Codes/Practical_Lattice_v2/3body_quantization/test_files/F3_for_pole_KKpi_L24/"

    testfile_drive = "/home/digonto/Codes/Practical_Lattice_v2/3body_quantization/test_files/"
    filename1 = drive + "ultraHQ_F3_for_pole_KKpi_L" + str(int(Lbyas)) + "_nP_000.dat"
    filename2 = drive + "ultraHQ_F3_for_pole_KKpi_L" + str(int(Lbyas)) + "_nP_100.dat"
    filename3 = drive + "ultraHQ_F3_for_pole_KKpi_L" + str(int(Lbyas)) + "_nP_110.dat"
    filename4 = drive + "ultraHQ_F3_for_pole_KKpi_L" + str(int(Lbyas)) + "_nP_111.dat"
    filename5 = drive + "ultraHQ_F3_for_pole_KKpi_L" + str(int(Lbyas)) + "_nP_200.dat"

    spectrum1 = drive + "non_int_L" + str(int(Lbyas)) + "_points_000_A1m.dat"
    spectrum2 = drive + "non_int_L" + str(int(Lbyas)) + "_points_100_A2.dat"
    spectrum3 = drive + "non_int_L" + str(int(Lbyas)) + "_points_110_A2.dat"
    spectrum4 = drive + "non_int_L" + str(int(Lbyas)) + "_points_111_A2.dat"
    spectrum5 = drive + "non_int_L" + str(int(Lbyas)) + "_points_200_A2.dat"



    atmpi = 1#0.06906

    #raulF3file = "detF3inv_test_FRL_L5"
    filelist = [filename1, filename2, filename3, filename4, filename5]
    non_int_filelist = [spectrum1, spectrum2, spectrum3, spectrum4, spectrum5] 

    titlelist = ["000","100","110","111","200"]
    fig,ax = plt.subplots(5,1,figsize=(16,45))
    counter = 0
    #energy_file = testfile_drive + "S2I2.energies"
    #full_energy_list, irrep_list = irrep_list_maker(energy_file)
    selected_irrep_list = ["000_A1m","100_A2","110_A2","111_A2","200_A2"]


    for file in filelist:
        print("file loaded = ",file)
        filename = file 
        (En, EcmR, norm, ReF3sum, F2, G, K2inv, Hinv) = np.genfromtxt(filename, unpack=True)

        (L_spec, non_int_energies_ecm) = np.genfromtxt(non_int_filelist[counter], unpack=True)
    
        np_zero_y_val = np.zeros((len(non_int_energies_ecm)))#np.array(zero_y_val)

    
        ax[counter].set_xlim(0.262,0.42)
        #ax[counter].set_ylim(-1E-17,1E-17)
        ax[counter].set_ylim(-0.0000051,0.0000051)
        ax[counter].tick_params(axis='both', which='major', labelsize=25)
        ax[counter].tick_params(axis='both', which='minor', labelsize=25)

        if((counter+1)!=len(filelist)):
            ax[counter].set_xticks([])
        #ax.set_yscale('log')
        #L = 6
        HFGtitle = "$(\mathcal{K}_2^{-1} + F_2 + G)^{-1}$"
        titlestring = "$P = $" + str(titlelist[counter])
        #ax[counter].set_title(titlestring)
        ax[counter].plot(EcmR,ReF3sum, linewidth=2, zorder=4,label=titlestring)
        ax[counter].scatter(non_int_energies_ecm,np_zero_y_val, marker='o', s=100, edgecolor="teal", facecolor='white',zorder=5)
        ax[counter].axhline(y=0,linestyle='--',color='black',zorder=2)

        ax[counter].axvline(x=0.37,linestyle='--',color='black',zorder=2)
    
        ax[4].set_xlabel("$E_{cm}$", fontsize=25)
        ax[counter].set_ylabel("$F_{3,\\textrm{iso}}$",fontsize=25)
        #plt.draw()
        counter=counter+1


    outputfile_str = "F3_KKpi_L20_v3.pdf"

    plt.show()
    #plt.savefig(outputfile_str)

    #plt.close()

    #print(full_energy_list[0][1][2])


plotter() 