#!/usr/env/bin python

# this script is used to evaluate dS/dlambda for a gaussian basis set (lambda = H)
# and update H based on the REM update rule:
# lambda_{i+1} = lambda_{i} - alpha * dS/dlambda / d2S/dL2 (newton raphson)
# lambda_{i+1} = lambda_{i} - alpha * dS/dlambda (steep descent)

import numpy as np
import math as math
import sys as sys

learning_rate = float(sys.argv[1])

ecap = 0.1 #max update kcal/mol
nframes = 4000.0
nmol = 18.0 #number of gag
nmol2 = 7.0 #number of ip6
rcut = 22.0
DEBUG = 0
t1DEBUG = "25"
t2DEBUG = "25"

ff_in=open("../gausswall.ff", "r")
ff_out=open("gausswall.ff", "w")
update_out=open("updates.dat", "w")

# E = H / (sigma * sqrt(2pi)) * exp( -(r-rh)^2 / (2*sigma^2) )
def eval_dEdH(r, sigma, rh) :
    #prefactor = 1.0 / (sigma * math.sqrt(2.0*math.pi))
    return math.exp( -(r - rh)**2.0 / (2.0 * sigma**2.0)) / (sigma * math.sqrt(2.0*math.pi))


all_lines = ff_in.readlines()
ff_in.close()

dSdlambda=0.0

for line in all_lines :

    strtok = line.split()
    type1 = strtok[1]
    type2 = strtok[2]
    H = float(strtok[4])
    rh = float(strtok[5])
    sigma = float(strtok[6])

    dSdL_temp_CG = 0.0
    dSdL_sq_temp_CG = 0.0
    count_temp_CG = 0.0
    
    rstats_file = "%s_%s.hist" % (type1, type2)
    hist_data = np.genfromtxt(rstats_file)

    if(len(hist_data) > 10) :
        for ii, (r, counts) in enumerate(hist_data) :
            dEdH = eval_dEdH(r, sigma, rh)
            if(r < rcut and counts > 0) :
                if (type1 == t1DEBUG and type2 == t2DEBUG and DEBUG) :
                    print "25 25 data: %f, %f, %f\n" % (r, counts, dEdH)
                dSdL_temp_CG += (dEdH * float(counts))
                dSdL_sq_temp_CG += (dEdH*dEdH * float(counts))
                count_temp_CG += float(counts)
        if (type1 == t1DEBUG and type2 == t2DEBUG and DEBUG) :
            print "25 25 total CG data: %f, %f, %f\n" % (dSdL_temp_CG, dSdL_sq_temp_CG, count_temp_CG)

    if(type2 != 36) :
        dSdL_temp_CG /= (nframes * nmol * (nmol-1))
        dSdL_sq_temp_CG /= (nframes * nmol * (nmol-1))
    else :
        dSdL_temp_CG /= (nframes * nmol * (nmol2))
        dSdL_sq_temp_CG /= (nframes * nmol * (nmol2))

    dSdL_temp_AA = 0.0
    count_temp_AA = 0.0

    rstats_file_AA = "../../RF/%s_%s.hist" % (type1, type2)
    hist_data = np.genfromtxt(rstats_file_AA)

    for ii, (r, counts) in enumerate(hist_data) :
        dEdH = eval_dEdH(r, sigma, rh)
        if(r < rcut and counts > 0) :
            dSdL_temp_AA += (dEdH * float(counts))
            count_temp_AA += float(counts)
    if (type1 == t1DEBUG and type2 == t2DEBUG and DEBUG) :
        print "25 25 total AA data: %f, %f\n" % (dSdL_temp_AA, count_temp_AA)
    if(type2 != 36) :
        dSdL_temp_AA /= (nframes * nmol * (nmol-1))
    else :
        dSdL_temp_AA /= (nframes * nmol * (nmol2))

    #dSdL = beta * (dUdL_{AA} - dUdL_{CG})
    #d2SdL2 = beta*(d2U/dL2_{AA} - d2U/dL2_{CG}) + beta^2 * ( ( (dU/dL)^2 )_{CG} - (dU/dL)_{CG}^2 )
    # lambda_{i+1} = lambda_{i} - alpha * dS/dlambda / d2S/dL2
    alpha = learning_rate
    beta = 0.60 #kcal/mol

    dSdL = beta * (dSdL_temp_AA - dSdL_temp_CG)
    d2SdL2 = beta * (0.0 - 0.0) + beta * beta * ( dSdL_sq_temp_CG - dSdL_temp_CG**2.0 )

    Hnew = 0.0
    update = 0.0
    steep_flag = 0
    if (steep_flag) :
        update = alpha * dSdL
    else :
        if(d2SdL2 > 0.0) :
            update = alpha * dSdL / d2SdL2
        elif(d2SdL2 <= 0.0) :
            update = 0.0 
        
        if(abs(update) > ecap) :
            update = np.sign(update) * ecap

    Hnew = H - update

    strtok[4] = "%3.4f" % Hnew
    ff_out.write(" ".join(strtok)+"\n")
    update_out.write("%3.4f %3.4f %3.4f\n" % (dSdL, d2SdL2, update))


ff_out.close()
update_out.close()
