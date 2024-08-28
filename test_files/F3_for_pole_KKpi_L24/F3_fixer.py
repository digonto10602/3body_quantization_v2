#This file fixes the F3 data files and generates one file for each irrep

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import pandas as pd 
from IPython.display import display

def F3_fixer_L24():
    Lbyas = 24

    moms_list = ["100","110","111","200"]

    highest_file_num = 5 # this number changes on what is the highest number 
                         # of files for an irrep that we will use to modify 
                         # the final F3 file 

    for moms in moms_list:
        filename_string1 = "ultraHQ_F3_for_pole_KKpi_L24_nP_" + moms + "_"
        file_list = []
        for i in range(0,highest_file_num,1):
            filename = filename_string1 + str(i) + ".dat"
            if(os.path.isfile(filename)):
                file_list.append(filename)

        En_tot = []
        Ecm_tot = []
        norm_vec_tot = []
        F3_tot = []
        F2_tot = []
        G_tot = []
        K2inv_tot = []
        Hinv_tot = []

        for i in range(0,len(file_list),1):
            file = file_list[i]
            print("using file = ",file)
            (En, Ecm, norm_vec, F3, F2, G, K2inv, Hinv) = np.genfromtxt(file, unpack=True)

            for j in range(0,len(Ecm),1):
                En_tot.append(En[j])
                Ecm_tot.append(Ecm[j]) 
                norm_vec_tot.append(norm_vec[j])
                F3_tot.append(F3[j])
                F2_tot.append(F2[j])
                G_tot.append(G[j])
                K2inv_tot.append(K2inv[j])
                Hinv_tot.append(Hinv[j])

        df = pd.DataFrame(list(zip(En_tot, Ecm_tot, norm_vec_tot, F3_tot, F2_tot, G_tot, K2inv_tot, Hinv_tot)))

        sorted_df = df.sort_values(1)
        re_sorted_df = sorted_df.reset_index(drop=True)
        print("sorting completed")
        outfile = "ultraHQ_F3_for_pole_KKpi_L24_nP_" + moms + ".dat"

        fout = open(outfile,"w")

        for i in range(0,len(En_tot),1):
            En_final = re_sorted_df[0][i]
            Ecm_final = re_sorted_df[1][i]
            norm_vec_final = re_sorted_df[2][i]
            F3_final = re_sorted_df[3][i]
            F2_final = re_sorted_df[4][i]
            G_final = re_sorted_df[5][i]
            K2inv_final = re_sorted_df[6][i]
            Hinv_final = re_sorted_df[7][i]

            fout.write(     str(En_final) + '\t'
                       +    str(Ecm_final) + '\t'
                       +    str(norm_vec_final) + '\t'
                       +    str(F3_final) + '\t'
                       +    str(F2_final) + '\t'
                       +    str(G_final) + '\t'
                       +    str(K2inv_final) + '\t'
                       +    str(Hinv_final) + '\n')

        fout.close()
        print("file = ",outfile,"generated.")    

F3_fixer_L24() 