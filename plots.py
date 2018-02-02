import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits
import aplpy

####################################
# PLOTTING ROUTINES
####################################

def plot_1(hi_flat,hi_cut,output):
    
    ####################################
    # LMC & YSO HEATMAP PLOTS (FIGURE 1)
    ####################################

    plt.close('all')

    fig = plt.figure(figsize=(11.25, 3.75))

    f1 = aplpy.FITSFigure(hi_flat,figure=fig, subplot=[0.1,0.1,0.275,0.8])
    f1.show_grayscale(stretch='log', vmin=1e20*0.8e-20,vmax=6e21*0.8e-20)
    f1.set_theme('publication')
    f1.tick_labels.set_xformat('ddd')
    f1.tick_labels.set_yformat('ddd')
    f1.recenter(79.0, -69.1, radius=5)

    # ALL candidate YSOs (Ochsendorf et al. 2016)
    sfr_ysolist = np.genfromtxt('image_data_files/sfr_radec_allyso.txt', delimiter=' ',usecols=[0,1],skip_header=0,filling_values=-999)
    x, y = sfr_ysolist[:,0],sfr_ysolist[:,1]
    
    # Calculate the point density for this dataset
    xy = np.vstack([x,y])
    z = stats.gaussian_kde(xy,bw_method=0.15)(xy)

    hi_cut.header['NAXIS'] = 2 
    del hi_cut.header['CDELT3'], hi_cut.header['CRPIX3'], hi_cut.header['CRVAL3'], hi_cut.header['CTYPE3']

    f2 = aplpy.FITSFigure(hi_cut,figure=fig, subplot=[0.45,0.1,0.275,0.8])
    f2.show_grayscale(stretch='log', vmin=2e-1,vmax=1.75)
    f2.set_theme('publication')
    f2.show_markers(x, y, c=z, s=10, facecolor='',edgecolor='',alpha=0.75,vmin=0.005)
    f2.tick_labels.set_xformat('ddd')
    f2.tick_labels.set_yformat('ddd')
    f2.recenter(79.0, -69.1, radius=5)
    f2.add_label(90.2, -64.3,r'Candidate YSOs',size=12,horizontalalignment='left') 

    # YSOs (Ochsendorf et al. 2016)
    sfr_ysolist = np.genfromtxt('image_data_files/sfr_radec_myso.txt', delimiter=' ',usecols=[0,1,2],skip_header=0,filling_values=-999)
    x, y,sfr_rate = sfr_ysolist[:,0],sfr_ysolist[:,1],sfr_ysolist[:,2]

    # Calculate the point density for this dataset
    xy = np.vstack([x,y])
    z = stats.gaussian_kde(xy,bw_method=0.15)(xy)

    f3 = aplpy.FITSFigure(hi_cut,figure=fig, subplot=[0.8,0.1,0.275,0.8])
    f3.show_grayscale(stretch='log', vmin=2e-1,vmax=1.75)
    f3.set_theme('publication')
    f3.show_markers(x, y, c=z, s=25, facecolor='',edgecolor='',alpha=0.7,vmin=0.01)
    f3.tick_labels.set_xformat('ddd')
    f3.tick_labels.set_yformat('ddd')
    f3.recenter(79.0, -69.1, radius=5)
    f3.add_label(90.2, -64.3,r'MYSOs ($M$ > 8 M$_\odot$)',size=12,horizontalalignment='left') 
    f3.add_label(81.2, -73.6,r'257 km/s - 270 km/s',size=12,horizontalalignment='left') 

    plt.savefig(output+'/figure1.pdf', bbox_inches='tight',dpi=225)


def plot_2(tryout,n79_par,n79_ind,n79_sfe,n79_sfe_mol,dor30_par,dor30_ind,dor30_sfe,dor30_sfe_mol,n11_par,n11_ind,n11_sfe,n11_sfe_mol,output):

    ####################################
    # SFR/SFE PLOTS (FIGURE 2)
    ####################################
    
    plt.close('all')
    
    fig = plt.figure(figsize=(10.5, 3))
    font = {'family' : 'normal','weight' : 'normal', 'size'   : 9}
    plt.rc('font', **font)
    plt.subplots_adjust(left=0.2, bottom=0.2, right=0.875, top=0.85, wspace=0.35, hspace=0.35)

    plt.subplot(131)
    plt.step(np.log10(tryout*5),np.log10(n79_par['sfr_myso']),label='N79',color='blue')
    plt.step(np.log10(tryout*5),np.log10(dor30_par['sfr_myso']),label='30Dor',color='green')
    plt.step(np.log10(tryout*5),np.log10(n11_par['sfr_myso']),label='N11',color='darkgray',linewidth=0.5,zorder=3)
    plt.errorbar(1.275,3.2,yerr=[0.1],color='black',capsize=2,linewidth=0.5)

    plt.scatter(np.log10(66),np.log10(n79_par['sfr_ha'][n79_ind]*2),marker='*',s=50,color='blue')
    plt.scatter(np.log10(102),np.log10(dor30_par['sfr_ha'][dor30_ind]*2),marker='*',s=50,color='green')
    plt.scatter(np.log10(150),np.log10(n11_par['sfr_ha'][n11_ind]*2),marker='*',s=50,facecolor='darkgray',edgecolor='darkgray')

    plt.errorbar(np.log10(66),np.log10(n79_par['sfr_ha'][n79_ind]*2),yerr=[0.24],color='blue',capsize=2,linewidth=0.5)
    plt.errorbar(np.log10(102),np.log10(dor30_par['sfr_ha'][dor30_ind]*2),yerr=[0.13],color='green',capsize=2,linewidth=0.5)
    plt.errorbar(np.log10(150),np.log10(n11_par['sfr_ha'][n11_ind]*2),yerr=[0.165],color='darkgray',capsize=2,linewidth=0.5)

    axes = plt.gca()
    axes.set_ylim([3,5])
    axes.set_xlim([1.2,2.4])
    plt.legend(loc='upper left',prop={'size':8},frameon=False)
    plt.xlabel(r'log[$R$ (pc)]')
    plt.ylabel(r'log[SFR (M$_\odot$ Myr$^{-1}$)]')
    plt.title('SFR')

    plt.subplot(132)
    plt.step(np.log10(tryout*5),np.log10(n79_par['sfr_myso']/(n79_par['j16'])),label='N79',color='blue')
    plt.step(np.log10(tryout*5),np.log10(dor30_par['sfr_myso']/(dor30_par['j16'])),label='30Dor',color='green')
    plt.step(np.log10(tryout*5),np.log10(n11_par['sfr_myso']/(n11_par['j16'])),label='N11',color='darkgray',linewidth=0.5,zorder=3)
    plt.errorbar(1.275,-2.65,yerr=[0.3],color='black',capsize=2,linewidth=0.5)

    plt.scatter(np.log10(66),np.log10(n79_sfe_mol),marker='*',s=50,color='blue')
    plt.scatter(np.log10(102),np.log10(dor30_sfe_mol),marker='*',s=50,color='green')
    plt.scatter(np.log10(150),np.log10(n11_sfe_mol),marker='*',s=50,facecolor='darkgray',edgecolor='darkgray')

    plt.errorbar(np.log10(66),np.log10(n79_sfe_mol),yerr=[np.sqrt(0.24**2+0.3**2)],color='blue',capsize=2,linewidth=0.5)
    plt.errorbar(np.log10(102),np.log10(dor30_sfe_mol),yerr=[np.sqrt(0.13**2+0.3**2)],color='green',capsize=2,linewidth=0.5)
    plt.errorbar(np.log10(150),np.log10(n11_sfe_mol),yerr=[np.sqrt(0.165**2+0.3**2)],color='darkgray',capsize=2,linewidth=0.5)

    axes = plt.gca()
    axes.set_ylim([-3.0,-0.5])
    axes.set_xlim([1.2,2.4])
    plt.legend(loc='upper right',prop={'size':8},frameon=False)
    plt.xlabel(r'log[$R$ (pc)]')
    plt.ylabel(r'log[SFR$_\mathrm{mYSO}$/$M_\mathrm{H_2}$ (Myr$^{-1}$)]')
    plt.title(r'SFE = SFR/$M_\mathrm{gas}$')

    plt.subplot(133)
    plt.step(np.log10(tryout*5),np.log10(n79_par['sfr_myso']/(n79_par['j16'] + n79_par['hi'] + n79_par['hii'])),label='N79',color='blue')
    plt.step(np.log10(tryout*5),np.log10(dor30_par['sfr_myso']/(dor30_par['j16'] + dor30_par['hi'] + dor30_par['hii'])),label='30Dor',color='green')
    plt.step(np.log10(tryout*5),np.log10(n11_par['sfr_myso']/(n11_par['j16'] + n11_par['hi'] + n11_par['hii'])),label='N11',color='darkgray',linewidth=0.5,zorder=3)
    plt.errorbar(1.275,-2.65,yerr=[0.3],color='black',capsize=2,linewidth=0.5)

    plt.scatter(np.log10(66),np.log10(n79_sfe),marker='*',s=50,color='blue')
    plt.scatter(np.log10(102),np.log10(dor30_sfe),marker='*',s=50,color='green')
    plt.scatter(np.log10(150),np.log10(n11_sfe),marker='*',s=50,facecolor='darkgray',edgecolor='darkgray')

    plt.errorbar(np.log10(66),np.log10(n79_sfe),yerr=[np.sqrt(0.24**2+0.3**2)],color='blue',capsize=2,linewidth=0.5)
    plt.errorbar(np.log10(102),np.log10(dor30_sfe),yerr=[np.sqrt(0.13**2+0.3**2)],color='green',capsize=2,linewidth=0.5)
    plt.errorbar(np.log10(150),np.log10(n11_sfe),yerr=[np.sqrt(0.165**2+0.3**2)],color='darkgray',capsize=2,linewidth=0.5)

    axes = plt.gca()
    axes.set_ylim([-3.0,-0.5])
    axes.set_xlim([1.2,2.4])
    plt.legend(loc='upper right',prop={'size':8},frameon=False)
    plt.xlabel(r'log[$R$ (pc)]')
    plt.ylabel(r'log[SFR$_\mathrm{mYSO}$/$M_\mathrm{H_2 + HI + HII}$ (Myr$^{-1}$)]')
    plt.title(r'SFE = SFR/$M_\mathrm{gas}$')

    plt.savefig(output+'/figure4.pdf', bbox_inches='tight',dpi=75)


def plot_s1(tryout,n79_vir,dor30_vir,n11_vir,output):

    ####################################
    # VIRIAL PARAMETER PLOT (FIGURE S1)
    ####################################

    plt.close('all')

    fig = plt.figure(figsize=(3, 3))

    plt.step(np.log10(tryout*5),np.log10(n79_vir),label='N79',color='blue') # 1 pix = 10" = 2.5 pc
    plt.step(np.log10(tryout*5),np.log10(dor30_vir),label='30Dor',color='green')
    plt.step(np.log10(tryout*5),np.log10(n11_vir),label='N11',color='darkgray')

    plt.axhline(y=0.3, color='k', linestyle='--')

    axes = plt.gca()
    axes.set_ylim([-0.25,2.5])
    axes.set_xlim([1,2.4])
    plt.legend(loc='upper right',prop={'size':8},frameon=False)
    plt.xlabel(r'log[$R$ (pc)]')
    plt.ylabel(r'log[$\alpha_\mathrm{vir}$]')
    plt.title('Virial parameter')

    plt.savefig(output+'/figure_s1.pdf', bbox_inches='tight',dpi=75)
    
    
def plot_s2(tryout,n79_par,dor30_par,n11_par,output):
    
    ####################################
    # ISM PLOT (FIGURE S2)
    ####################################

    fig = plt.figure(figsize=(8, 8))
    font = {'family' : 'normal','weight' : 'normal', 'size'   : 10}
    plt.rc('font', **font)

    plt.subplots_adjust(left=0.2, bottom=0.2, right=0.875, top=0.85, wspace=0.35, hspace=0.4)

    plt.subplot(221)
    plt.step(np.log10(tryout*5),np.log10(n79_par['co']),label='N79',color='blue') # 1 pix = 10" = 2.5 pc
    plt.step(np.log10(tryout*5),np.log10(dor30_par['co']),label='30Dor',color='green')
    plt.step(np.log10(tryout*5),np.log10(n11_par['co']),label='N11',color='darkgray')
    plt.errorbar(1.1,-2.65,yerr=[0.3],color='black',capsize=2,linewidth=0.5)
    plt.errorbar(2.3,3.4,yerr=[0.3],color='black',capsize=2,linewidth=0.5)

    axes = plt.gca()
    axes.set_ylim([3,7])
    axes.set_xlim([1,2.4])
    plt.legend(loc='upper left',prop={'size':8},frameon=False)
    plt.xlabel(r'log[$R$ (pc)]')
    plt.ylabel(r'log[$M_\mathrm{H2,CO}$ (M$_\odot$)]')
    plt.title('Molecular gas (CO)')

    plt.subplot(222)
    plt.step(np.log10(tryout*5),np.log10(n79_par['j16']),label='N79',color='blue') # 1 pix = 56" = 14 pc
    plt.step(np.log10(tryout*5),np.log10(dor30_par['j16']),label='30Dor',color='green')
    plt.step(np.log10(tryout*5),np.log10(n11_par['j16']),label='N11',color='darkgray')
    plt.errorbar(2.3,3.4,yerr=[0.3],color='black',capsize=2,linewidth=0.5)

    axes = plt.gca()
    axes.set_ylim([3,7])
    axes.set_xlim([1,2.4])
    plt.legend(loc='upper left',prop={'size':8},frameon=False)
    plt.xlabel(r'log[$R$ (pc)]')
    plt.ylabel(r'log[$M_\mathrm{H2,dust}$ (M$_\odot$)]')
    plt.title('Molecular gas (dust-based)')

    plt.subplot(223)
    plt.step(np.log10(tryout*5),np.log10(n79_par['hi']),label='N79',color='blue')
    plt.step(np.log10(tryout*5),np.log10(dor30_par['hi']),label='30Dor',color='green')
    plt.step(np.log10(tryout*5),np.log10(n11_par['hi']),label='N11',color='darkgray')
    plt.errorbar(2.3,3.4,yerr=[0.1],color='black',capsize=2,linewidth=0.5)

    axes = plt.gca()
    axes.set_ylim([3,7])
    axes.set_xlim([1,2.4])
    plt.legend(loc='upper left',prop={'size':8},frameon=False)
    plt.xlabel(r'log[$R$ (pc)]')
    plt.ylabel(r'log[$M_\mathrm{HI}$ (M$_\odot$)]')
    plt.title('Neutral atomic gas')

    plt.subplot(224)
    plt.step(np.log10(tryout*5),np.log10(n79_par['hii']),label='N79',color='blue')
    plt.step(np.log10(tryout*5),np.log10(dor30_par['hii']),label='30Dor',color='green')
    plt.step(np.log10(tryout*5),np.log10(n11_par['hii']),label='N11',color='darkgray')
    plt.errorbar(2.3,3.4,yerr=[0.175],color='black',capsize=2,linewidth=0.5)
    
    axes = plt.gca()
    axes.set_ylim([3,7])
    axes.set_xlim([1,2.4])
    plt.legend(loc='upper left',prop={'size':8},frameon=False)
    plt.xlabel(r'log[$R$ (pc)]')
    plt.ylabel(r'log[$M_\mathrm{HII}$ (M$_\odot$)]')
    plt.title('Ionized gas')

    plt.savefig(output+'/figure_s2.pdf', bbox_inches='tight',dpi=75)
    
