#This file copies the mass jack files from the directories and puts it
#over here, the two sections are divided based on the L values we are using

#This also generates the Ecm_data.ini file needed for scattering devel 
#to generate the spectrum file in Ecm frame 

import numpy as np 
import sys
import subprocess
import os.path


#Input 
L20_irrep_list = ["000_A1m","100_A2","110_A2","111_A2","200_A2"]

#L20 and L24 irrep list are the same 
L24_irrep_list = L20_irrep_list 

L20_ensemble_name = "szscl21_20_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265"
L24_ensemble_name = 'szscl21_24_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265'

L20_mass_directory = "../szscl21_20_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per/fits_mhi/"  #000_A1m/out/t0_9/prin_corrs/ord*

L24_mass_directory = "../szscl21_24_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265/fits_mhi/"

L20_t0_list = ["9", #000_A1m
               "10",#100_A2
               "13",#110_A2
               "8", #111_A2
               "10" #200_A2
              ]

L24_t0_list = ["12", #000_A1m
               "12",#100_A2
               "12",#110_A2
               "12", #111_A2
               "12" #200_A2
              ]

#first copy the L20 files 

gap = " "
max_state = 30
entry_for_ecmini = []

for i in range(0,len(L20_irrep_list),1):
    for j in range(0,max_state,1):
        state_num = j
        filename =  ( L20_mass_directory + L20_irrep_list[i] 
                    + "/out/t0_" + L20_t0_list[i] 
                    + "/prin_corrs/ord" + str(state_num) 
                    + "/m.jack"  )
        searching_file = filename 
        #print("searching file = ",searching_file)
        if(os.path.exists(filename)):
            renamed_mass_file = L20_ensemble_name + "_" + L20_irrep_list[i] + "_state_" + str(state_num) 
            print("found file = ",renamed_mass_file)

            runcommand = "cp" + gap + filename + gap + renamed_mass_file

            subprocess.call([runcommand],shell=True, text=True, executable='/bin/bash')
            
            split_irrep = L20_irrep_list[i].split("_",1)
            boost = split_irrep[0]
            irrep = split_irrep[1]
            ecm_entry =  (   "20" + '\t'
                        +   boost + '\t'
                        +   irrep + '\t'
                        +   str(state_num) + '\t'
                        +   renamed_mass_file + '\t'
                        +   "[[]]" + '\t'
                        +   "{{}}" 
                        )
            entry_for_ecmini.append(ecm_entry)

for i in range(0,len(L24_irrep_list),1):
    for j in range(0,max_state,1):
        state_num = j
        filename =  ( L24_mass_directory + L24_irrep_list[i] 
                    + "/out/t0_" + L24_t0_list[i] 
                    + "/prin_corrs/ord" + str(state_num) 
                    + "/m.jack"  )
        searching_file = filename 
        #print("searching file = ",searching_file)
        if(os.path.exists(filename)):
            renamed_mass_file = L24_ensemble_name + "_" + L24_irrep_list[i] + "_state_" + str(state_num) 
            print("found file = ",renamed_mass_file)

            runcommand = "cp" + gap + filename + gap + renamed_mass_file

            subprocess.call([runcommand],shell=True, text=True, executable='/bin/bash')

            split_irrep = L24_irrep_list[i].split("_",1)
            boost = split_irrep[0]
            irrep = split_irrep[1]
            ecm_entry =  (   "24" + '\t'
                        +   boost + '\t'
                        +   irrep + '\t'
                        +   str(state_num) + '\t'
                        +   renamed_mass_file + '\t'
                        +   "[[]]" + '\t'
                        +   "{{}}" 
                        )
            print("ecm entry = ",ecm_entry)
            entry_for_ecmini.append(ecm_entry)


ecm_ini_file = "Ecm_data.ini"
f = open(ecm_ini_file,"w")

for i in range(0,len(entry_for_ecmini),1):
    string1 = entry_for_ecmini[i]

    f.write(string1 + '\n')

f.close() 




        