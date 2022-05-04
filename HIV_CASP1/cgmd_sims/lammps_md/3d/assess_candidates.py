import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def lognorm(x, mu, sig) :
    return (1/(x*sig*(2*np.pi)**0.5) * np.exp( -0.5 * ((np.log(x)-mu)/sig)**2 ))

cdata = np.zeros(10)
for i in range(10) :
    fopen = open("trial_%d/log.1" % (i+1),"r")
    all_lines = fopen.readlines()
    fopen.close()

    for l in all_lines :
        strtok = l.split()
        if len(strtok) < 21 :
            continue
        if strtok[0] == "50000000" :
            cdata[i] = float(strtok[20])
            break

histo, bin_edges = np.histogram(cdata, bins=20, density=True)
bin_centers = 0.5*(bin_edges[-2] + bin_edges[1:])

# filter out zero values
nz_inds = histo != 0
histo = histo[nz_inds]
bin_centers = bin_centers[nz_inds]

cmax = np.max(cdata)
cmin = np.min(cdata)
m0 = np.mean(np.log(cdata))
s0 = np.std(np.log(cdata))

p0 = [ m0, s0 ] 
popt, pcov = curve_fit(lognorm, 
                       bin_centers, 
                       histo, 
                       p0=p0,
                       method='lm')

# calculate mean
mean_c = np.exp( popt[0] + 0.5*popt[1]**2.0 )

# find trial closest to mean
err_c = 10000
trial_mean = 1
for i, c in enumerate(cdata) :
    cdiff = np.abs( c - mean_c )
    if cdiff < err_c :
        err_c = cdiff
        trial_mean = i+1
print(trial_mean)
