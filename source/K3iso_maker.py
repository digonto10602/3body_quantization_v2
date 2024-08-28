# K3iso maker and fitter

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PyPDF2 import PdfMerger

import os 
import subprocess 

def E_to_Ecm(En, P):
    return np.sqrt(En**2 - P**2)

def Esq_to_Ecmsq(En, P):
    return En**2 - P**2

def Ecmsq_to_Esq(Ecm, P):
    return Ecm**2 + P**2

def P000():
    drive = "/home/digonto/Codes/Practical_Lattice_v2/KKpi_interacting_spectrum/Lattice_data/KKpi_L20/"
    filename1 = drive + "mass_000_A1m_t011_MassJackFiles_mass_t0_11_reorder_state0.jack"
    filename2 = drive + "mass_000_A1m_t011_MassJackFiles_mass_t0_11_reorder_state1.jack"
    filename3 = drive + "mass_000_A1m_t011_MassJackFiles_mass_t0_11_reorder_state2.jack"
    filename4 = drive + "mass_000_A1m_t011_MassJackFiles_mass_t0_11_reorder_state3.jack"
    filename5 = drive + "mass_000_A1m_t011_MassJackFiles_mass_t0_11_reorder_state4.jack"

    filelist = [filename1, filename2, filename3, filename4, filename5]

    
    xi = 3.444
    pi = np.arccos(-1.0)
    L = 20
    twopibyxiL = 2.0*pi/(xi*L)
    nPx = 0
    nPy = 0
    nPz = 0

    Px = nPx*twopibyxiL
    Py = nPy*twopibyxiL
    Pz = nPz*twopibyxiL

    P = np.sqrt(Px*Px + Py*Py + Pz*Pz)

    statecounter=0

    for files in filelist:
        zerocol, energylab = np.genfromtxt(files,skip_header=1,unpack=True)

        outputfile = "K3iso_P" + str(nPx) + str(nPy) + str(nPz) + "_state" + str(statecounter) + ".dat"

        fout = open(outputfile,"w")



        for i in range(len(energylab)):
            K3iso = subprocess.check_output(['./generate_K3iso',str(nPx),str(nPy),str(nPz),str(energylab[i])],shell=False)
            result = K3iso.decode('utf-8')
            finresult = float(result)
            Ecm = E_to_Ecm(energylab[i],P)
            print(nPx,nPy,nPz,statecounter, energylab[i],E_to_Ecm(energylab[i],P),float(result))
            output = str(energylab[i]) + '\t' + str(Ecm) + '\t' + str(finresult) + '\n' 
            #output = output.rstrip('\n')
            fout.write(output)
    
        fout.close()
        statecounter = statecounter + 1

def P100():
    drive = "/home/digonto/Codes/Practical_Lattice_v2/KKpi_interacting_spectrum/Lattice_data/KKpi_L20/"
    filename1 = drive + "mass_100_A2_t012_MassJackFiles_mass_t0_12_reorder_state0.jack"
    filename2 = drive + "mass_100_A2_t012_MassJackFiles_mass_t0_12_reorder_state1.jack"
    filename3 = drive + "mass_100_A2_t012_MassJackFiles_mass_t0_12_reorder_state2.jack"
    filename4 = drive + "mass_100_A2_t012_MassJackFiles_mass_t0_12_reorder_state3.jack"
    filename5 = drive + "mass_100_A2_t012_MassJackFiles_mass_t0_12_reorder_state4.jack"
    filename6 = drive + "mass_100_A2_t012_MassJackFiles_mass_t0_12_reorder_state5.jack"

    filelist = [filename1, filename2, filename3, filename4, filename5, filename6]

    
    xi = 3.444
    pi = np.arccos(-1.0)
    L = 20
    twopibyxiL = 2.0*pi/(xi*L)
    nPx = 0
    nPy = 0
    nPz = 1

    Px = nPx*twopibyxiL
    Py = nPy*twopibyxiL
    Pz = nPz*twopibyxiL

    P = np.sqrt(Px*Px + Py*Py + Pz*Pz)

    statecounter=0

    for files in filelist:
        zerocol, energylab = np.genfromtxt(files,skip_header=1,unpack=True)

        outputfile = "K3iso_P" + str(nPx) + str(nPy) + str(nPz) + "_state" + str(statecounter) + ".dat"

        fout = open(outputfile,"w")



        for i in range(len(energylab)):
            K3iso = subprocess.check_output(['./generate_K3iso',str(nPx),str(nPy),str(nPz),str(energylab[i])],shell=False)
            result = K3iso.decode('utf-8')
            finresult = float(result)
            Ecm = E_to_Ecm(energylab[i],P)
            print(nPx,nPy,nPz,statecounter, energylab[i],E_to_Ecm(energylab[i],P),float(result))
            output = str(energylab[i]) + '\t' + str(Ecm) + '\t' + str(finresult) + '\n' 
            #output = output.rstrip('\n')
            fout.write(output)
    
        fout.close()
        statecounter = statecounter + 1

def P110():
    drive = "/home/digonto/Codes/Practical_Lattice_v2/KKpi_interacting_spectrum/Lattice_data/KKpi_L20/"
    filename1 = drive + "mass_110_A2_t013_MassJackFiles_mass_t0_13_reorder_state0.jack"
    filename2 = drive + "mass_110_A2_t013_MassJackFiles_mass_t0_13_reorder_state1.jack"
    filename3 = drive + "mass_110_A2_t013_MassJackFiles_mass_t0_13_reorder_state2.jack"
    filename4 = drive + "mass_110_A2_t013_MassJackFiles_mass_t0_13_reorder_state3.jack"
    filename5 = drive + "mass_110_A2_t013_MassJackFiles_mass_t0_13_reorder_state4.jack"
    filename6 = drive + "mass_110_A2_t013_MassJackFiles_mass_t0_13_reorder_state5.jack"

    filelist = [filename1, filename2, filename3, filename4, filename5, filename6]

    
    xi = 3.444
    pi = np.arccos(-1.0)
    L = 20
    twopibyxiL = 2.0*pi/(xi*L)
    nPx = 0
    nPy = 1
    nPz = 1

    Px = nPx*twopibyxiL
    Py = nPy*twopibyxiL
    Pz = nPz*twopibyxiL

    P = np.sqrt(Px*Px + Py*Py + Pz*Pz)

    statecounter=0

    for files in filelist:
        zerocol, energylab = np.genfromtxt(files,skip_header=1,unpack=True)

        outputfile = "K3iso_P" + str(nPx) + str(nPy) + str(nPz) + "_state" + str(statecounter) + ".dat"

        fout = open(outputfile,"w")



        for i in range(len(energylab)):
            K3iso = subprocess.check_output(['./generate_K3iso',str(nPx),str(nPy),str(nPz),str(energylab[i])],shell=False)
            result = K3iso.decode('utf-8')
            finresult = float(result)
            Ecm = E_to_Ecm(energylab[i],P)
            print(nPx,nPy,nPz,statecounter, energylab[i],E_to_Ecm(energylab[i],P),float(result))
            output = str(energylab[i]) + '\t' + str(Ecm) + '\t' + str(finresult) + '\n' 
            #output = output.rstrip('\n')
            fout.write(output)
    
        fout.close()
        statecounter = statecounter + 1

def P111():
    drive = "/home/digonto/Codes/Practical_Lattice_v2/KKpi_interacting_spectrum/Lattice_data/KKpi_L20/"
    filename1 = drive + "mass_111_A2_t012_MassJackFiles_mass_t0_12_reorder_state0.jack"
    filename2 = drive + "mass_111_A2_t012_MassJackFiles_mass_t0_12_reorder_state1.jack"
    filename3 = drive + "mass_111_A2_t012_MassJackFiles_mass_t0_12_reorder_state2.jack"
    filename4 = drive + "mass_111_A2_t012_MassJackFiles_mass_t0_12_reorder_state3.jack"
    filename5 = drive + "mass_111_A2_t012_MassJackFiles_mass_t0_12_reorder_state4.jack"
    filename6 = drive + "mass_111_A2_t012_MassJackFiles_mass_t0_12_reorder_state5.jack"

    filelist = [filename1, filename2, filename3, filename4, filename5, filename6]

    
    xi = 3.444
    pi = np.arccos(-1.0)
    L = 20
    twopibyxiL = 2.0*pi/(xi*L)
    nPx = 1
    nPy = 1
    nPz = 1

    Px = nPx*twopibyxiL
    Py = nPy*twopibyxiL
    Pz = nPz*twopibyxiL

    P = np.sqrt(Px*Px + Py*Py + Pz*Pz)

    statecounter=0

    for files in filelist:
        zerocol, energylab = np.genfromtxt(files,skip_header=1,unpack=True)

        outputfile = "K3iso_P" + str(nPx) + str(nPy) + str(nPz) + "_state" + str(statecounter) + ".dat"

        fout = open(outputfile,"w")



        for i in range(len(energylab)):
            Ecmsq = Esq_to_Ecmsq(energylab[i],P)
            
            if(Ecmsq<0.26**2 or Ecmsq>0.43**2):
                continue 
            else:
                Ecm = np.sqrt(Ecmsq)
                K3iso = subprocess.check_output(['./generate_K3iso',str(nPx),str(nPy),str(nPz),str(energylab[i])],shell=False)
            result = K3iso.decode('utf-8')
            finresult = float(result)
            
            print(nPx,nPy,nPz,statecounter, energylab[i],E_to_Ecm(energylab[i],P),float(result))
            output = str(energylab[i]) + '\t' + str(Ecm) + '\t' + str(finresult) + '\n' 
            #output = output.rstrip('\n')
            fout.write(output)

        fout.close()
        statecounter = statecounter + 1

def P200():
    drive = "/home/digonto/Codes/Practical_Lattice_v2/KKpi_interacting_spectrum/Lattice_data/KKpi_L20/"
    filename1 = drive + "mass_200_A2_t010_MassJackFiles_mass_t0_10_reorder_state0.jack"
    filename2 = drive + "mass_200_A2_t010_MassJackFiles_mass_t0_10_reorder_state1.jack"
    filename3 = drive + "mass_200_A2_t010_MassJackFiles_mass_t0_10_reorder_state2.jack"
    filename4 = drive + "mass_200_A2_t010_MassJackFiles_mass_t0_10_reorder_state3.jack"
    filename5 = drive + "mass_200_A2_t010_MassJackFiles_mass_t0_10_reorder_state4.jack"
    filename6 = drive + "mass_200_A2_t010_MassJackFiles_mass_t0_10_reorder_state5.jack"
    filename7 = drive + "mass_200_A2_t010_MassJackFiles_mass_t0_10_reorder_state6.jack"
    filename8 = drive + "mass_200_A2_t010_MassJackFiles_mass_t0_10_reorder_state7.jack"

    filelist = [filename1, filename2, filename3, filename4, filename5, filename6, filename7, filename8]

    
    xi = 3.444
    pi = np.arccos(-1.0)
    L = 20
    twopibyxiL = 2.0*pi/(xi*L)
    nPx = 0
    nPy = 0
    nPz = 2

    Px = nPx*twopibyxiL
    Py = nPy*twopibyxiL
    Pz = nPz*twopibyxiL

    P = np.sqrt(Px*Px + Py*Py + Pz*Pz)

    statecounter=0

    for files in filelist:
        zerocol, energylab = np.genfromtxt(files,skip_header=1,unpack=True)

        outputfile = "K3iso_P" + str(nPx) + str(nPy) + str(nPz) + "_state" + str(statecounter) + ".dat"

        fout = open(outputfile,'w')

        

        for i in range(len(energylab)):
            Ecmsq = Esq_to_Ecmsq(energylab[i],P)
            
            if(Ecmsq<0.26**2 or Ecmsq>0.43**2):
                continue 
            else:
                Ecm = np.sqrt(Ecmsq)
                K3iso = subprocess.check_output(['./generate_K3iso',str(nPx),str(nPy),str(nPz),str(energylab[i])],shell=False)
            result = K3iso.decode('utf-8')
            finresult = float(result)
            #Ecm = E_to_Ecm(energylab[i],P)
            print(nPx,nPy,nPz,statecounter, energylab[i],E_to_Ecm(energylab[i],P),float(result))
            output = str(energylab[i]) + '\t' + str(Ecm) + '\t' + str(finresult) + '\n' 
            #output = output.rstrip('\n')
            fout.write(output)
    
        fout.close()
        statecounter = statecounter + 1

#this one calculates the K3iso using the averaged data file
def jackknifeavg_lattice_data():
    #this drive path has been set from macbook
    drive = "/Users/digonto/GitHub/3body_quantization/lattice_data/KKpi_interacting_spectrum/Three_body/"
    #this drive path is for ubuntu
    drive = "/home/digonto/Codes/Practical_Lattice_v2/3body_quantization/lattice_data/KKpi_interacting_spectrum/Three_body/"
    
    filename1 = drive + "KKpi_spectrum.000_A1m"
    filename2 = drive + "KKpi_spectrum.100_A2"
    filename3 = drive + "KKpi_spectrum.110_A2"
    filename4 = drive + "KKpi_spectrum.111_A2"
    filename5 = drive + "KKpi_spectrum.200_A2"

    filelist = [filename1, filename2, filename3, filename4, filename5]

    
    xi = 3.444
    pi = np.arccos(-1.0)
    L = 20
    twopibyxiL = 2.0*pi/(xi*L)
    


    nconfig = [[0,0,0],[0,0,1],[0,1,1],[1,1,1],[0,0,2]]

    
    statecounter=0
    counter = 0
    
    for files in filelist:
        Ls, Ecm, err1, err2 = np.genfromtxt(files,unpack=True)

        nPx = nconfig[counter][0]
        nPy = nconfig[counter][1]
        nPz = nconfig[counter][2]

        Px = nPx*twopibyxiL
        Py = nPy*twopibyxiL
        Pz = nPz*twopibyxiL

        P = np.sqrt(Px*Px + Py*Py + Pz*Pz)


        

        for i in range(len(Ecm)):
            Ecm_val_ini = Ecm[i] - err2[i]
            Ecm_val_fin = Ecm[i] + err2[i]
            Elab_ini = np.sqrt(Ecmsq_to_Esq(Ecm_val_ini,P))
            Elab_fin = np.sqrt(Ecmsq_to_Esq(Ecm_val_fin,P))
            Elab = np.linspace(Elab_ini, Elab_fin, 500)

            outputfile = "K3iso_jackavg_P" + str(nPx) + str(nPy) + str(nPz) + "_state_" + str(i) + ".dat"

            fout = open(outputfile,"w")

            for j in range(len(Elab)):
                Elab_val = Elab[j]
                K3iso = subprocess.check_output(['./generate_K3iso',str(nPx),str(nPy),str(nPz),str(Elab_val)],shell=False)
                result = K3iso.decode('utf-8')
                finresult = float(result)
                calcEcmsq = Esq_to_Ecmsq(Elab_val,P)
                calcEcm = np.sqrt(calcEcmsq)
                print(nPx,nPy,nPz,i, Elab_val,calcEcm,float(result))
                output = str(Elab_val) + '\t' + str(calcEcm) + '\t' + str(finresult) + '\n' 
                #output = output.rstrip('\n')
                fout.write(output)
    
            fout.close()
        statecounter = statecounter + 1
        counter = counter + 1 

def jackknifeavg_centralvalue_lattice_data():
    #this drive path has been set from macbook
    drive = "/Users/digonto/GitHub/3body_quantization/lattice_data/KKpi_interacting_spectrum/Three_body/"
    #this drive path is for ubuntu
    drive = "/home/digonto/Codes/Practical_Lattice_v2/3body_quantization/lattice_data/KKpi_interacting_spectrum/Three_body/"
    
    filename1 = drive + "KKpi_spectrum.000_A1m"
    filename2 = drive + "KKpi_spectrum.100_A2"
    filename3 = drive + "KKpi_spectrum.110_A2"
    filename4 = drive + "KKpi_spectrum.111_A2"
    filename5 = drive + "KKpi_spectrum.200_A2"

    filelist = [filename1, filename2, filename3, filename4, filename5]

    
    xi = 3.444
    pi = np.arccos(-1.0)
    L = 20
    twopibyxiL = 2.0*pi/(xi*L)
    


    nconfig = [[0,0,0],[0,0,1],[0,1,1],[1,1,1],[0,0,2]]

    
    statecounter=0
    counter = 0
    
    for files in filelist:
        Ls, Ecm, err1, err2 = np.genfromtxt(files,unpack=True)

        nPx = nconfig[counter][0]
        nPy = nconfig[counter][1]
        nPz = nconfig[counter][2]

        Px = nPx*twopibyxiL
        Py = nPy*twopibyxiL
        Pz = nPz*twopibyxiL

        P = np.sqrt(Px*Px + Py*Py + Pz*Pz)

        outputfile = "K3iso_jackavg_centralval_P" + str(nPx) + str(nPy) + str(nPz) + ".dat"

        fout = open(outputfile,"w")

        

        for i in range(len(Ecm)):
            
            Ecm_val = Ecm[i]
            Elab_val = np.sqrt(Ecmsq_to_Esq(Ecm_val,P))
            K3iso = subprocess.check_output(['./generate_K3iso',str(nPx),str(nPy),str(nPz),str(Elab_val)],shell=False)
            result = K3iso.decode('utf-8')
            finresult = float(result)
            calcEcmsq = Esq_to_Ecmsq(Elab_val,P)
            calcEcm = np.sqrt(calcEcmsq)
            print(nPx,nPy,nPz,i, Elab_val,calcEcm,float(result))
            output = str(Elab_val) + '\t' + str(calcEcm) + '\t' + str(finresult) + '\n' 
            #output = output.rstrip('\n')
            fout.write(output)                  
            
            
    
        fout.close()
        statecounter = statecounter + 1
        counter = counter + 1 


#P000()
#P100()
#P110()
#P111()
#P200()

jackknifeavg_lattice_data()
jackknifeavg_centralvalue_lattice_data()




############################### TEST ############################################

#################################################################################
#file = open("temp.temp",'w')

#for i in range(10):
#    file.write(str(i) + '\n')

#file.close()
#command = ['./generate_K3iso','0','0','0','0.235373']    
#args = [float(arg) if arg.isdigit() else float(arg) for arg in command[1:]]
#checkres = subprocess.check_output(command,shell=False, stderr=subprocess.STDOUT)
#print(checkres.decode('ascii'))
##################################################################################

######################################################################################################################################
#hello = "hiaaa"
#output = subprocess.call(["echo " + str(hello) + " |awk '{print(substr($0,1,3)}'"]) #subprocess.check_output("echo ")

#output = subprocess.check_output(['/usr/bin/echo ', hello])

#output = subprocess.run(["echo",hello,"|awk '{print(substr($0,1,1)}'"],capture_output=True, shell=True)

#output1 = output.decode('ascii')

#print(output)
#print(output.stdout.decode())
#print(float(output1) + 2)
######################################################################################################################################
