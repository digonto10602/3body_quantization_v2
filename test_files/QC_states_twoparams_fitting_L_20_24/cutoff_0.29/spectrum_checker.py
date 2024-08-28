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


def QC(F3, K3):
    return 1.0/F3 + K3 

def plotter():

    plt.rcParams.update({'font.size': 18})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    Lbyas = 24

    K3iso0 = 589836.99644428
    K3iso1 = -41238939.44689889
    cutoff = 0.29

    kkpi=0.26302

    drive = "/home/digonto/Codes/Practical_Lattice_v2/3body_quantization/test_files/F3_for_pole_KKpi_L24/"
    QC_spec_drive = "/home/digonto/Codes/Practical_Lattice_v2/3body_quantization/test_files/QC_states_twoparams_fitting_L_20_24/cutoff_" + str(cutoff) + "/"
    spec_drive = "/home/digonto/Codes/Practical_Lattice_v2/3body_quantization/lattice_data/KKpi_interacting_spectrum/twoptvar_analysis/masses/"

    filename1 = drive + "ultraHQ_F3_for_pole_KKpi_L" + str(int(Lbyas)) + "_nP_000.dat"
    filename2 = drive + "ultraHQ_F3_for_pole_KKpi_L" + str(int(Lbyas)) + "_nP_100.dat"
    filename3 = drive + "ultraHQ_F3_for_pole_KKpi_L" + str(int(Lbyas)) + "_nP_110.dat"
    filename4 = drive + "ultraHQ_F3_for_pole_KKpi_L" + str(int(Lbyas)) + "_nP_111.dat"
    filename5 = drive + "ultraHQ_F3_for_pole_KKpi_L" + str(int(Lbyas)) + "_nP_200.dat"

    nonintspectrum1 = drive + "non_int_L" + str(int(Lbyas)) + "_points_000_A1m.dat"
    nonintspectrum2 = drive + "non_int_L" + str(int(Lbyas)) + "_points_100_A2.dat"
    nonintspectrum3 = drive + "non_int_L" + str(int(Lbyas)) + "_points_110_A2.dat"
    nonintspectrum4 = drive + "non_int_L" + str(int(Lbyas)) + "_points_111_A2.dat"
    nonintspectrum5 = drive + "non_int_L" + str(int(Lbyas)) + "_points_200_A2.dat"

    spec_file1 = spec_drive + "V" + str(int(Lbyas)) + "_000_A1m"
    spec_file2 = spec_drive + "V" + str(int(Lbyas)) + "_100_A2"
    spec_file3 = spec_drive + "V" + str(int(Lbyas)) + "_110_A2"
    spec_file4 = spec_drive + "V" + str(int(Lbyas)) + "_111_A2"
    spec_file5 = spec_drive + "V" + str(int(Lbyas)) + "_200_A2"

    QC_state1 = QC_spec_drive + "QC_states_jackknifed_L_" + str(int(Lbyas)) + "_nP_000_energycutoff_" + str(cutoff) + "_two_params.dat"
    QC_state2 = QC_spec_drive + "QC_states_jackknifed_L_" + str(int(Lbyas)) + "_nP_100_energycutoff_" + str(cutoff) + "_two_params.dat"
    QC_state3 = QC_spec_drive + "QC_states_jackknifed_L_" + str(int(Lbyas)) + "_nP_110_energycutoff_" + str(cutoff) + "_two_params.dat"
    QC_state4 = QC_spec_drive + "QC_states_jackknifed_L_" + str(int(Lbyas)) + "_nP_111_energycutoff_" + str(cutoff) + "_two_params.dat"
    QC_state5 = QC_spec_drive + "QC_states_jackknifed_L_" + str(int(Lbyas)) + "_nP_200_energycutoff_" + str(cutoff) + "_two_params.dat"


    atmpi = 1#0.06906

    #raulF3file = "detF3inv_test_FRL_L5"
    filelist = [filename1, filename2, filename3, filename4, filename5]
    non_int_filelist = [nonintspectrum1, nonintspectrum2, nonintspectrum3, nonintspectrum4, nonintspectrum5] 
    spec_filelist = [spec_file1, spec_file2, spec_file3, spec_file4, spec_file5]
    QC_spec_filelist = [QC_state1, QC_state2, QC_state3, QC_state4, QC_state5]

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

        (nonint_L_spec, non_int_energies_ecm) = np.genfromtxt(non_int_filelist[counter], unpack=True)

        (L_QC, E_QC, Eerr_QC, diffQC, differrQC) = np.genfromtxt(QC_spec_filelist[counter], unpack=True)

        (L_spec, E_spec, E_err1_spec, E_err2_spec) = np.genfromtxt(spec_filelist[counter], unpack=True)
    
        np_zero_y_val = np.zeros((len(non_int_energies_ecm)))#np.array(zero_y_val)

        QC_val = []

        for i in range(0,len(ReF3sum),1):
            s = EcmR[i]**2 
            Msq = kkpi**2 
            K3iso = K3iso0 + (s - Msq)*K3iso1 
            QC_val.append(QC(ReF3sum[i],K3iso))

        np_QC_val = np.array(QC_val) 

    
        ax[counter].set_xlim(0.262,0.42)
        ax[counter].set_ylim(-1E7,1E7)
        #ax[counter].set_ylim(-0.0000051,0.0000051)
        ax[counter].tick_params(axis='both', which='major', labelsize=25)
        ax[counter].tick_params(axis='both', which='minor', labelsize=25)

        if((counter+1)!=len(filelist)):
            ax[counter].set_xticks([])
        #ax.set_yscale('log')
        #L = 6
        HFGtitle = "$(\mathcal{K}_2^{-1} + F_2 + G)^{-1}$"
        titlestring = "$P = $" + str(titlelist[counter])
        #ax[counter].set_title(titlestring)
        ax[counter].plot(EcmR,np_QC_val, linewidth=2, zorder=4,label=titlestring)
        ax[counter].scatter(non_int_energies_ecm,np_zero_y_val, marker='o', s=100, edgecolor="teal", facecolor='white',zorder=5)
        ax[counter].scatter(E_spec,np.zeros((len(E_spec))), marker='o', s=100, edgecolor="darkred", facecolor='orange',zorder=5)
        ax[counter].scatter(E_QC,np.zeros((len(E_QC))), marker='o', s=100, edgecolor="blue", facecolor='indigo',zorder=5)
        
        ax[counter].axhline(y=0,linestyle='--',color='black',zorder=2)

        ax[counter].axvline(x=0.37,linestyle='--',color='black',zorder=2)

        #ax[2].scatter(0.3053831258637795,0,marker='s',s=100)
        #ax[2].scatter(0.3085861914201294,0,marker='s',s=100)
    
        ax[4].set_xlabel("$E_{cm}$", fontsize=25)
        ax[2].set_ylabel("$F_{3,\\textrm{iso}}^{-1} + \\mathcal{K}_{3,\\textrm{iso}}$",fontsize=25)
        #plt.draw()
        counter=counter+1


    outputfile_str = "F3_KKpi_L20_v3.pdf"

    plt.show()
    #plt.savefig(outputfile_str)

    #plt.close()

    #print(full_energy_list[0][1][2])


plotter() 