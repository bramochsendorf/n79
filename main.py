import pandas as pd
import numpy as np
import scipy.stats as stats
from astropy.io import fits
from sys import exit
from astropy.wcs import WCS
from plots import plot_1,plot_2,plot_s1,plot_s2

####################################
# CONSTANTS
####################################


syst_unc = 0.3 # systematic uncertainty in dex of x-factor (MAGMA) (recommonded by Bolatto 2013) and molecular mass (Jameson et. al 2016)
m_h = 1.6738232e-24 # mass of proton
pccm = 3.086e18 # pc to cm
sol_mass = 2e33 # solar mass in gram
grav = 6.67e-8 # gravitational constant in dyn cm^-2 g^-2
myr = 365.25*24*3600*1e6 # Myr in seconds
hz = 1.249e13 # 24 micron in Hz
aspx = 3136 # the amount of arcseconds per pix (56"x56")
srpx = 2.3504e-11 * 3136 # the amount of steradian per pix


####################################
# FUNCTIONS
####################################


def radius(xpix,ypix,xcenter,ycenter):
    """
    calculate the radius from (xpix, ypix) to (xcenter, ycenter)
    """
    rad = np.sqrt((xpix-xcenter)**2+(ypix-ycenter)**2)
    
    return rad
    
def get_par(xcenter_cube,ycenter_cube,radius,array,wcs_mol,pix_size,tracer):
    """
    calculate SFR_myso and ISM parameters as a function of aperture size
    """
    
    xcube0,ycube0,zcube0 = w.all_pix2world(xcenter_cube,ycenter_cube,0,0)
    radius2 = radius*(20./pix_size) # everything normalized to the HI pixel size (20")
    a,b = wcs_mol.all_world2pix(xcube0,ycube0,0)
    nx,ny = array.shape
    y,x = np.ogrid[-a:nx-a,-b:ny-b]
    
    circle = np.where(x*x + y*y <= radius2*radius2)
    
    if tracer == 'ha' or tracer == 'um24':
        mol_mass = np.nansum(array.data[circle[1],circle[0]])
    elif tracer == 'disp':
        array.data[array.data == 0] = np.nan    
        mol_mass = np.nanmean(array.data[circle[1],circle[0]]) # return average dispersion value.
    else:
        mol_mass = np.nansum(array.data[circle[1],circle[0]])*(pix_size/4.)**2
        
    return mol_mass
    
def calc_par(xcen,ycen,xpix,ypix,apertures,data_hi,data_co,data_j16,data_hii,data_ha,data_um24,wcs_hi,wcs_co,wcs_j16,wcs_hii,wcs_ha,wcs_um24):
    """
    calculate SFR_myso and ISM parameters as a function of aperture size
    """    
    
    distance = radius(xpix,ypix,xcen,ycen)
    
    # pixel sizes from input maps (in arcseconds)
    pixsize_hi = 20.
    pixsize_co = 10.
    pixsize_j16 = 56.
    pixsize_hii = 15.
    pixsize_ha = 56.
    pixsize_um24 = 56.
    
    d = []
    
    for apsize in apertures:
        # sfr_myso is calculated by summing up masses and 0.5 myr formation timescale
        sfr_myso = np.sum(sfr_rate[np.where(distance < apsize)])/0.5
        hi = get_par(xcen,ycen,apsize,data_hi,w_flat,pixsize_hi,'hi')
        co = get_par(xcen,ycen,apsize,data_co,wcs_co,pixsize_co,'co')
        j16 = get_par(xcen,ycen,apsize,data_j16,wcs_j16,pixsize_j16,'j16')
        hii = get_par(xcen,ycen,apsize,data_hii,wcs_hii,pixsize_hii,'hii')
        ha = get_par(xcen,ycen,apsize,data_ha,w_ha,pixsize_ha,'ha')
        um24 = get_par(xcen,ycen,apsize,data_um24,w_ha,pixsize_um24,'um24')   
        disp = (get_par(xcen,ycen,apsize,data_co,wcs_co,pixsize_co,'co')*get_par(xcen,ycen,apsize,co_mom2,w_co,pixsize_co,'disp'))/get_par(xcen,ycen,apsize,data_co,wcs_co,pixsize_co,'co')
        # append everything through a list comprehension
        d.append({'sfr_myso':sfr_myso, 'hi':hi, 'co':co, 'j16':j16, 'hii':hii, 'ha':ha, 'um24':um24, 'disp':disp })
    
    # save to dataframe
    result = pd.DataFrame(d)
     
    return result

def sfr(ha,um24):
    """
    calculate SFR from Ha and 24 micron (Calzetti et al. 2007)
    """
    l_ir = 4.*np.pi*(5e4*3.086e18)**2*(um24*srpx*1e-23*1e6*hz) # IR (input map in MJy/sr)
    l_ha = 4.*np.pi*(5e4*3.086e18)**2*(ha*5.661e-18*aspx*0.1) # halpha (0.1 comes from decirayleigh)        
    em_final = 5.3e-42*(l_ha + 0.032*l_ir)*1e6
    
    return em_final 
    
def vir(par,apertures):
    """
    calculate virial parameter. Factor 5 comes from 20" ~ 5 pc @ the distance of the LMC
    """
    disp = par['disp']
    co = par['co']
    vir_par =(5*(disp*1e5)**2*(aperture*5*pccm))/(grav*(co)*sol_mass)
    
    return vir_par

def sfr_sfe(par,apertures,rad):
    """
    calculate star formation efficiencies (only molecular and total gas)
    adjust SFR by factor of 2 (Pellegrini et al. 2011; Ochsendorf et al. 2017)
    factor 5 comes from 20" per pixel ~ 5 pc at the distance of the LMC
    first extract relevant columns from dataframe
    """
    sfr_ha = par['sfr_ha']
    j16 = par['j16']
    hi = par['hi']
    hii = par['hii']
    
    apertures = apertures*5.
    index = np.argmin(np.abs(apertures - rad))
    sfe = sfr_ha[index]*2./(j16[index]+hi[index]+hii[index])
    sfe_mol = sfr_ha[index]*2./(j16[index])
    
    return index,sfe,sfe_mol


####################################
# DATA IMPORT
####################################


# Sum of layers 42 to 50 in the original HI cube to highlight LMC's spiral arms
hi_cut = fits.open('image_data_files/lmc_spiralarms.fits')[0]

# HI column density image: reshape such that aplpy can read it, remove third axis (including keywords), and transform to Msun/pc^-2 (factor 0.8e-20)
hi_flat = fits.open('image_data_files/lmc_HI.fits')[0]
hi_flat.data = np.reshape(hi_flat.data,[hi_flat.shape[1],hi_flat.shape[2]]) * 0.8e-20
hi_flat.header['NAXIS'] = 2 
del hi_flat.header['CDELT3'], hi_flat.header['CRPIX3'], hi_flat.header['CRVAL3'], hi_flat.header['CTYPE3']

# CO-based molecular gas: multiply by alpha_co factor (8.6; Bollato et al. 2013)
co = fits.open('image_data_files/co_mom1_10as.fits')[0]
co.data = co.data * 8.6

# Dust-based molecular gas (courtesy Katie Jameson); mask negative values
j16 = fits.open('image_data_files/LMC_H2_bembb_gdr500it2_beam_clean_MF3pix.fits')[0]
bad = np.where(j16 < 0)
j16.data[bad] = 0 

# Ionized gas surface density 
hii = fits.open('image_data_files/lmc_ionized_halpha_surf_dens.fits')[0]

# CO velocity dispersion
co_mom2 = fits.open('image_data_files/co_mom2_10as.fits')[0]

# Halpha (SHASSA) + 24micron (SAGE)
ha = fits.open('image_data_files/ha.fits')[0]
um24 = fits.open('image_data_files/smoothed_24um.fits')[0]

# subtract global background from HI and HII map (average of four corners away from the LMC)
bg = 3.5e20*0.8e-20
hi_flat.data = hi_flat.data - bg

bg_hii = 0.2
hii.data = hii.data - bg_hii

# read in YSO RA, DEC, and SFR rate (Ochsendorf et al. 2016)
sfr_ysolist = np.genfromtxt('image_data_files/sfr_radec_myso.txt', delimiter=' ',usecols=[0,1,2],skip_header=0,filling_values=-999)
x, y,sfr_rate = sfr_ysolist[:,0],sfr_ysolist[:,1],sfr_ysolist[:,2]

# read in WCS data from all images
w = WCS(hi_cut.header)
w_flat = WCS(hi_flat.header)
w_co = WCS(co.header)
w_j16 = WCS(j16.header)
w_hii = WCS(hii.header)
w_ha = WCS(ha.header)
w_um24 = WCS(um24.header)

xpix,ypix,zpix = w.all_world2pix(x,y,0,0)

####################################
# CALCULATIONS 
####################################

# define incremental aperture sizes to measure quantities in
aperture = np.logspace(0,2,50)

# Calculate the probability density function of YSOs based on location
xy = np.vstack([x,y])
z = stats.gaussian_kde(xy,bw_method=0.15)(xy)

# Center points for star forming regions
x_n79,y_n79 = 714,507 # N79 
x_dor30,y_dor30 = 343,547 # 30 Dor
x_n11,y_n11 = 703,778 # N11

# calculate SFR_myso and ISM masses for N79, 30Dor, and N11. Returns a Pandas dataframe.
n79_par = calc_par(x_n79,y_n79,xpix,ypix,aperture,hi_flat,co,j16,hii,ha,um24,w_flat,w_co,w_j16,w_hii,w_ha,w_um24)
dor30_par = calc_par(x_dor30,y_dor30,xpix,ypix,aperture,hi_flat,co,j16,hii,ha,um24,w_flat,w_co,w_j16,w_hii,w_ha,w_um24)
n11_par = calc_par(x_n11,y_n11,xpix,ypix,aperture,hi_flat,co,j16,hii,ha,um24,w_flat,w_co,w_j16,w_hii,w_ha,w_um24)

# calculate SFR_halpha + 24micron. Add to dataframe.
n79_par['sfr_ha'] = sfr(n79_par['ha'],n79_par['um24'])
dor30_par['sfr_ha'] = sfr(dor30_par['ha'],dor30_par['um24'])
n11_par['sfr_ha'] = sfr(n11_par['ha'],n11_par['um24'])

# calculate virial parameter. Add to dataframe.
n79_par['vir'] = vir(n79_par,aperture)
dor30_par['vir'] = vir(dor30_par,aperture)
n11_par['vir'] = vir(n11_par,aperture)

# calculate SFE for total gas (j16+hi+hii) and molecular gas only (j16)   
# HII region radii from Lopez et al. 2010; should contain 90% of all flux. N79 = 66 pc, 30 Dor = 102 pc, N11 = 150 pc 
n79_ind,n79_sfe,n79_sfe_mol = sfr_sfe(n79_par,aperture,66)
dor30_ind,dor30_sfe,dor30_sfe_mol = sfr_sfe(dor30_par,aperture,102)
n11_ind,n11_sfe,n11_sfe_mol = sfr_sfe(n11_par,aperture,150)


####################################
# MAKE PLOTS
####################################

output = '/data1/python/output/paper3/github/'

plot_1(hi_flat,hi_cut,output)
plot_2(aperture,n79_par,n79_ind,n79_sfe,n79_sfe_mol,dor30_par,dor30_ind,dor30_sfe,dor30_sfe_mol,n11_par,n11_ind,n11_sfe,n11_sfe_mol,output)
plot_s1(aperture,n79_par['vir'],dor30_par['vir'],n11_par['vir'],output)
plot_s2(aperture,n79_par,dor30_par,n11_par,output)
