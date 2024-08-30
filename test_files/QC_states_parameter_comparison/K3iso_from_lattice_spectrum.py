#Here we plot the K3iso found from the quantization condition
#using the lattice spectrum, we will use the interpolator to 
#draw the -F3inv that encompasses the lattice energy value within 
#their uncertainties. Along with that we will also plot the K3iso 
#value we found from the fitting. 

#This function fits the
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

def E_to_Ecm(En, P):
    return np.sqrt(En**2 - P**2)

def Esq_to_Ecmsq(En, P):
    return En**2 - P**2

def Ecmsq_to_Esq(Ecm, P):
    return Ecm**2 + P**2

def QC3(K3iso, F3inv):
    return F3inv + K3iso 

def sign_func(val):
    if(val>0.0):
        return 1.0; 
    else:
        return -1.0; 