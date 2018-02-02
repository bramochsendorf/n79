import numpy as np
import itertools
from astropy.coordinates import SkyCoord
from astropy import units as u


####################################
# FUNCTIONS
####################################


def Kroupa(m, p1=0.3, p2=1.3, p3=2.3, break1=0.08, break2=0.5, integral_form=False):
    """ 
    Kroupa 2001 IMF 
    """
    m = np.array(m)
    binv = ((break1**(-(p1-1)) - 0.03**(-(p1-1)))/(1-p1) +
            (break2**(-(p2-1)) - break1**(-(p2-1))) * (break1**(p2-p1))/(1-p2) +
            (- break2**(-(p3-1))) * (break1**(p2-p1)) * (break2**(p3-p2))/(1-p3))
    b = 1./binv
    c = b * break1**(p2-p1)
    d = c * break2**(p3-p2)
    zeta = (b*(m**(-(p1))) * (m<break1) +
            c*(m**(-(p2))) * (m>=break1) * (m<break2) +
            d*(m**(-(p3))) * (m>=break2))

    return zeta*m

def integrate(massfunc, mlow, mhigh, **kwargs):
    """
    Integrate the mass function over some range
    """
    import scipy.integrate

    return scipy.integrate.quad(massfunc, mlow, mhigh, **kwargs)


####################################
# DATA INPUT
####################################


# This file contains the photometry for YSOs (source ID 0 - 3524) and dust clumps (source ID 3525 - 4698), formatted for the Robitaille 2006 SED fitter. Only the first 3524 are fitted with the SED fitter (because of limited photometry for the dust clumps). PACS100 is left out of fitting due to high level of artifacts. 
data_sage = np.genfromtxt('image_data_files/yso_list/catalogue_merged_yso_clumps_uplimit.txt', delimiter=' ',usecols=[0,1,2,30],skip_header=0,filling_values=-999)

# This file contains the ouput of the SED fitter for all YSO models. Some SEDs did not converge with the fitting and have been cleaned from the data.
# with the cleaned data file that doesn't contain '0 fit' lines, the lines are alternating and we can use itertools to slice the datafile and read in all model parameters:
with open('image_data_files/yso_list/merged_ysos_clumps_nopacs100_uplimit.txt') as f_in:
    source_pre = np.loadtxt(itertools.islice(f_in, 3, 9967, 2),skiprows=3,dtype=str)

with open('image_data_files/yso_list/merged_ysos_clumps_nopacs100_uplimit.txt') as f_in:
    source_para_pre = np.loadtxt(itertools.islice(f_in, 4, 9968, 2),skiprows=3,dtype=str)


####################################
# DATA MANIPULATION
####################################


# separate good and bad fits.
# threshold of maximum chi-squared value, completeness in solar masses per year
thres = 8.0 
comp = 8.

# grab the model fits below this goodness-of-fit threshold
source_raw = source_pre[np.where(np.float_(source_para_pre[:,2])/(np.float_(source_pre[:,1])) <= thres)]
source_para = source_para_pre[np.where(np.float_(source_para_pre[:,2])/np.float_(source_pre[:,1]) <= thres)]

# grab the model parameters from good fits
mass = np.float_(source_para[:,6])
mdot = np.float_(source_para[:,9])
mdisk = np.float_(source_para[:,13])
chi_sq = np.float_(source_para[:,2])
lum = np.float_(source_para[:,29])

# now limit only to >8 solar mass objects (completeness limit (Ochsendorf et al. 2016))
complete = np.where(mass >= comp)

source = source_raw[complete]
mass_com = mass[complete]
mdot_com = mdot[complete]
mdisk_com = mdisk[complete]

# calculate total mass using Kroupa IMF and YSO mass, assuming the YSO dominates the luminosity.
mass_imf = []
for yso_mass in mass_com:
    mass_kroupa = np.logspace(-2, 2.5,500) # mass to integrate
    massive = np.argmin((np.abs(mass_kroupa-yso_mass))) # couple mass of source to index 
    kroupa = Kroupa(mass_kroupa) # build kroupa function 
    factor = 1/kroupa[massive] # scale kroupa such that there is exactly 1 source at source mass
    totalmass = integrate(Kroupa,0.08,yso_mass) 
    sf = totalmass[0]*factor
    mass_imf.append(sf)


correct_index = np.zeros(0)
for i in np.arange(0,len(complete[0]-1)):
    correct_index = np.hstack([correct_index,source[i][0]])

# separate into stage 1/2/3 following Robitaille et al. 2006
stage_1 = np.where(mdot_com/mass_com > 10**(-6))        
stage_2 = np.where(np.logical_and(mdot_com/mass_com < 10**(-6), mdisk_com/mass_com > 10**(-6)))
stage_3 = np.where(np.logical_and(mdot_com/mass_com < 10**(-6), mdisk_com/mass_com < 10**(-6)))



####################################
# DATA OUTPUT
####################################


# output results in text file
# first grab RA and DEC from original file
fit_ra = data_sage[:,1]
fit_dec = data_sage[:,2] 

# three columns: RA, DEC, and mass
ra_dec_mass = np.zeros(len(fit_dec[complete][stage_1]), dtype=[('RA',float),('DEC',float),('cl_mass',float)])
ra_dec_mass['RA'], ra_dec_mass['DEC'], ra_dec_mass['cl_mass'] = fit_ra[np.int_(correct_index)][stage_1],fit_dec[np.int_(correct_index)][stage_1],np.float_(mass_imf)[stage_1]

# we now adjust mass of H72.97-69.39 as it is severely underestimated (did not include Far-IR). Total luminosity from SED fit (aperture_photometry) is single ~150 Msol star. For massive stars, the evolutionary track is nearly horizontal, so the current luminosity should be close zo ZAMS mass. Then with an IMF we find totall mass ~ 3112 Msun
ra_dec_mass['cl_mass'][60] = 3117

np.savetxt('image_data_files/sfr_radec_myso.txt', ra_dec_mass, fmt="%s")
