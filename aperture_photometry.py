import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import CircularAnnulus
from sys import exit 
from astropy.wcs import WCS

plt.close('all')


####################################
# FUNCTIONS
####################################


def magtojy(mag,zero):
    """
    function to transform magnitude to flux (Jy)
    """
    jy = 10.**(-0.4*mag)*zero
    return jy

def func(nu, scale1,scale2,beta1,beta2, T1, T2):
    """
    fit a two-temperature blackbody
    """
    return mod_blackbody_lam(nu, scale1, beta1, T1) + mod_blackbody_lam(nu, scale2, beta2, T2)

def mod_blackbody_lam(nu, scale, beta, T):
    """ 
    blackbody as a function of frequency (nu) and temperature (K).
    returns units of erg/s/cm^2/cm/Steradian
    """
    from scipy.constants import h,k,c
    k_0 = 4.8 # cm^2 g^-1
    nu_0 = 1.2e12 # 250 micron = 1.2 THz
    k_abs = k_0 * (nu/nu_0)**beta
    return scale * k_abs * 2*h*nu**3 / (c**2 * (np.exp(h*nu / (k*T)) - 1)) * 1e23 # 1e20 to convert to Jy
       
       
####################################
# DATA IMPORT
####################################


# image names
images = ['irac_36','irac_45','irac_58','irac_80','wise_22','mips_70','pacs_100','pacs_160','spire_250','spire_350','spire_500']

# read in first image to get dimensions
irac_3 = fits.open('/data1/Data/lmc_smc/n79/cutout/irac_36.fits')[0]

# fill up cube
data_cube = []

for im in images:
    implane = fits.open('/data1/Data/lmc_smc/n79/cutout/'+im+'.fits')[0]
    data_cube.append(implane)
  
  
####################################
# DATA MANIPULATION
####################################


# cube parameters (wavelength of images, resolution in arcseconds, and pixel sizes)
wave = np.float_([3.6,4.5,5.8,8,22,70,100,160,250,350,500]) 
fwhm = [2.,2.,2.,2.,12.,20.,7.,12.,18.2,24.9,36.3]
pix = [2.,2.,2.,2.,1.375,4.8,2.,3.,6.,10.,14]

# define poistion and aperture for photometry
positions = [72.971945, -69.391315] # S-1 (H72.97)
source = 'S-1'

# perform photometry, both the science and bg apertures.
sed = np.zeros(0)
for i in range(0,len(images)):
    # get WCS for each image
    slice = data_cube[i]
    w = WCS(slice.header) 
    # position of source in pixel
    pos = w.all_world2pix(positions[0],positions[1], 0) # position of aperture in pixel
   
    # source is very extended at shorter wavelengths, so adjust annuli for these waves.
    if i < 5:
        apertures = CircularAperture(pos, r=(5*fwhm[i])/(pix[i])) # take aperture of diameter x*FWHM
        annulus_apertures = CircularAnnulus(pos,  r_in=(7.*fwhm[i])/(pix[i]), r_out=(10.*fwhm[i])/(pix[i]))
        apers = [apertures, annulus_apertures]
    else:
        apertures = CircularAperture(pos, r=(3*fwhm[i])/(2.*pix[i])) # take aperture of diameter x*FWHM
        annulus_apertures = CircularAnnulus(pos,  r_in=(7.*fwhm[i])/(2.*pix[i]), r_out=(10.*fwhm[i])/(2.*pix[i]))
        apers = [apertures, annulus_apertures]
        
    phot_table = aperture_photometry(slice, apers)
    
    # area of apertures are different; normalize bg by (1) divide the sum in the annulus by its area
    bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area() 
    # total sum is then mean times area
    bkg_sum = bkg_mean * apertures.area() 
    # subtract bg
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    # incorporate in SED
    sed = np.hstack([sed,final_sum]) 

# multiply by 0.0046 (MJy/sr --> Jy/pix)
sed = sed * (np.float_(pix)/206265)**2.*1e6 

# Ks from ISRF, transformed to 2MASS system with Nakajima et al. 2008
mag = [14.53]

# zeropoints in Jy for 2MASS (http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html)
zero = [666.7]

# magnitude --> flux
ubvijhk = map(magtojy,mag,zero)

# transform lambda (microns) to nu
nu_wave = 3e14/wave

wave2 = np.float_([2.190])
nu_wave2 = 3e14/wave2
sed2 = np.float_([0.689,0.032])

# fit with two-component blackbody
# get the best fitting parameter values and their 1 sigma errors
param_bounds=([-np.inf,-np.inf,0.5,0.5,-np.inf,-np.inf],[np.inf,np.inf,2.5,2.5,np.inf,np.inf])
popt, pcov = curve_fit(func, nu_wave, sed, p0=(1e-4,1e-4, 1.0, 1.0,200,60), sigma=0.2*sed,bounds=param_bounds)

best = popt
sigma = np.sqrt(np.diag(pcov))

# set up xaxis and best fit models
nu = np.logspace(9,15.5,1000)
ybest1 = mod_blackbody_lam(nu, best[0],best[2],best[4]) 
ybest2 = mod_blackbody_lam(nu, best[1],best[3],best[5]) 

# integrate model flux
int_flux = np.trapz(ybest1+ybest2,x=nu) # the pdf weighted average (expectation value of the PDF)

# calculate total luminosity (4*pi*d^2)*flux, normalized by L_sol and Jy-->erg (3.85e33*1e23). Output result
dist = 5e4*3.086e18
total_lum = (4.*np.pi*dist**2*int_flux)/(3.85e33*1e23)
print 'Total luminosity (1e6 Lsun) = ', total_lum/1e6


####################################
# PLOTTING ROUTINES
####################################


fig = plt.figure(figsize=(6, 3))
font = {'family' : 'normal','weight' : 'normal', 'size'   : 12}
plt.rc('font', **font)

plt.subplot(111)
plt.subplots_adjust(left=0.2, bottom=0.2, right=0.875, top=0.85, wspace=0.35, hspace=0.3)

plt.scatter(wave,nu_wave*sed*1e-23,c='k')
plt.errorbar(wave,nu_wave*sed*1e-23, yerr=0.2*nu_wave*sed*1e-23, fmt='o',c='k')
[plt.scatter(wave2,nu_wave2*ubvijhk*1e-23, edgecolor='black',facecolor='black', marker='o', s=25, alpha=1, zorder=2) for x in ubvijhk] 

plt.plot(3e14/nu, nu*ybest1*1e-23, label='Best fitting\nmodel',ls='--',c='darkgray',linewidth=1.5)
plt.plot(3e14/nu, nu*ybest2*1e-23, label='Best fitting\nmodel',ls='--',c='lightgray',linewidth=1.5)
plt.plot(3e14/nu, nu*(ybest1+ybest2)*1e-23, label='Best fitting\nmodel',c='k',linewidth=1.5)

plt.text(1.25e2,2.5e-8,r'$L_\mathrm{IR}$ = '+'{:3.1e}'.format(total_lum)+r' M$_\odot$',fontsize=12)

axes = plt.gca()
plt.xlabel(r'$\lambda$ ($\mu$m)')
plt.ylabel(r'$\nu$$F_\mathrm{\nu}$ (erg s$^{-1}$ cm$^{-2}$)')
axes.set_xlim([0.5,5000])
axes.set_ylim([1e-12,1e-7])   
axes.set_xscale("log")
axes.set_yscale("log")
plt.savefig('/data1/python/output/paper3/github/h7297_sedfit.pdf',dpi=75)
    
exit()
    
# save spectrum for JWST ETC (micron, mJy)  
h72_sed = np.zeros(len(nu), dtype=[('var1',float),('var2',float)])
h72_sed['var1'], h72_sed['var2']= 3e14/nu[::-1], (ybest1[::-1]+ybest2[::-1])
np.savetxt('/data1/python/output/paper3/github/h72_sed.txt', h72_sed, fmt="%s")

    
    
    