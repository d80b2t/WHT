"""
This is a script to make several plots of the WHT IR 'CLQ' data from
2019June28 and 2019June29.
"""

import random
import datetime
import numpy   as np
import seaborn as sns
import matplotlib.pyplot   as plt
import matplotlib.patches  as patches
import matplotlib.lines    as mlines
import matplotlib.gridspec as gridspec

from matplotlib                      import colors as mcolors
from matplotlib.ticker               import ScalarFormatter
from matplotlib.patches              import Rectangle
from matplotlib.gridspec             import GridSpec
from matplotlib.backends.backend_pdf import PdfPages

from astropy.io          import ascii
from astropy.io          import fits
from astropy.convolution import convolve, Box1DKernel

#import pylustrator
#pylustrator.start()

##  R E A D I N G   I N   T H E   D A T A

##  L I G H T   C U R V E S
print(' --- Reading in AllWISE, NEOWISE-R data.... ' )

## Directory with the data 
path    = '/cos_pc19a_npr/programs/quasars/WHT/data/NEOWISER_LCs/'

## AllWISE
infile    =  'WHT_IR_AllWISE.tbl'
AllWISE   = ascii.read(path+infile)
upto_this = AllWISE['cntr_01'].max()
## NEOWISE-R
infile    = 'WHT_CLQs_NEOWISER_LCs.tbl'
NEOWISER  = ascii.read(path+infile)
## Averages of the e.g. NEOWISE-R L1b light curve points
infile  = 'WHT_CLQs_NEOWISER_LCs_averaged.tbl'
NEOWISER_aver = ascii.read(path+infile)


##  S P E C T R A 
path    = '/cos_pc19a_npr/programs/quasars/WHT/data/spectra/'
#infile  = 'spectra_20190628.lis'
infile = 'spectra.lis'
objects = ascii.read(path+infile)

object_name = objects['ObjectName']
redshift    = objects['Redshift']

## Emission Line list
linelist_file = 'emission_lines.dat'
linelist = ascii.read(path+linelist_file)


## Create the PdfPages object to which we will save the pages:
## The with statement makes sure that the PdfPages object is closed properly at
## the end of the block, even if an Exception occurs.
##pdf_pages = PdfPages('my-fancy-document.pdf')
#pdfs       = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
#out_pdf     = 'SDSS_WHT_IR_CLQs_temp.pdf'
#pdf_pages   = PdfPages(out_pdf)


#with PdfPages('SDSS_WHT_IR_CLQs_temp.pdf') as pdf
for ii in range(len(objects)):
    print(ii, 'SDSS_spectra/'+str(object_name[ii])+'.fits')
    #print(ii, 'WHTSpec_20190628/'+str(object_name[ii])+'.fits')
    print(ii, 'WHT_spectra/'+str(object_name[ii])+'.fits')   

            
    sdssname        = 'SDSS_spectra/'+str(object_name[ii])+'.fits'
    sdss_data       = fits.open(path+sdssname)
    sdss_spectrum   = sdss_data[1].data
    sdss_flux       = sdss_spectrum.flux
    sdss_loglam     = sdss_spectrum.loglam
    sdss_wavelength = 10**(sdss_loglam)
    sdss_mjd        = sdss_data[2].data['MJD'][0]

    print(sdss_flux.mean())
    
    #whtname        = 'WHTSpec_20190628/'+str(object_name[ii])+'.fits'
    whtname        = 'WHT_spectra/'+str(object_name[ii])+'.fits'
    wht_data       = fits.open(path+whtname)
    wht_flux       = wht_data[0].data / (1e-17)
    wht_len        = len(wht_data[0].data)
    # CRVAL1  =       3100.041015625 / Reference value on 1st axis in primary WCS     
    # CRPIX1  =                   1. / Reference pixel on 1st axis in primary WCS
    # CD1_1   =    0.861641567346351 / Transformation matrix for primary WCS
    #
   	# iraf.imgets(image, "CRVAL1")
	# start = iraf.imgets.value
	# iraf.imgets(image, "CD1_1")
	# delta = iraf.imgets.value
	# iraf.imgets(image, "CRPIX1")
	# offset = iraf.imgets.value
	# start, delta, offset = float(start), float(delta), float(offset)
    # wavelength += [ start + delta * ( i + 1 - offset ) ]

    start  = wht_data[0].header['CRVAL1']
    delta  = wht_data[0].header['CD1_1']
    offset = wht_data[0].header['CRPIX1']
    wavelength = []
	#for lam in range(wht_len):
    #    wavelength += [ start + delta * ( lam + 1 - offset ) ]

    ## CDELT1  =    0.861700942293611                                                  
    wht_wavelength = (np.arange(wht_len) + 3600.)*0.861700942293611


    ##  S E T T I N G     U P    T H E    P L O T 
    #fig, ax = plt.subplots(figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')
    fig = plt.figure(figsize=(12, 4.3))
    gs = GridSpec(2, 5, wspace=0.6, hspace=0.3)

    ax1 = fig.add_subplot(gs[0:, 0:2]) 
    ax2 = fig.add_subplot(gs[0, 2:5])
    ax3 = fig.add_subplot(gs[-1, 2])
    ax4 = fig.add_subplot(gs[-1, 3])
    ax5 = fig.add_subplot(gs[-1, -1])

    #fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
    #                               gridspec_kw={'height_ratios': [2, 1]}, 
    #                               figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')
    
    ## Adjusting the Whitespace for the plots
    left   = 0.06   # the left side of the subplots of the figure
    right  = 0.96   # the right side of the subplots of the figure
    bottom = 0.16   # the bottom of the subplots of the figure
    top    = 0.96   # the top of the subplots of the figure
    wspace = 0.26   # the amount of width reserved for blank space between subplots
    hspace = 0.26   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    ## Some NPR defaults
    lw              = 1.0
    ls              = 'solid'
    ms              = 5.
    ms_big          = ms*6.
    ls              = 'solid'
    alpha           = 1.0
    fontsize        = 16
    labelsize       = fontsize
    tickwidth       = 2.0
    linewidth       = 2.4
    tickwidth       = 2.0
    ticklength      = 6.0
    ticklabelsize   = labelsize
    majorticklength = 12
    minorticklength = 6

    ##  T H E    L I G H T     C U R V E S 
    xmin = 51400; xmax = 59000   ## SDSS EDR MJD.min() = 51608 ## was previously xmin = 53500. 
    ymin = 18.29; ymax = 14.71   
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin, ymax])

    ## Plotting the SPECTRA as vertical lines
    lw = 2.0
    ax1.axvline(x=58662, linewidth=lw, linestyle='dashed', color='r')
    
    jj = ii+1
    AllWISE_plot        = AllWISE[np.where(AllWISE['cntr_01']   == jj)]
    AllWISE_plot_W1_AB  = AllWISE_plot['w1mpro'] + 2.673
    AllWISE_plot_W2_AB  = AllWISE_plot['w2mpro'] + 3.313
    
    NEOWISER_plot       = NEOWISER[np.where(NEOWISER['cntr_01'] == jj)]
    NEOWISER_W1_AB      = NEOWISER_plot['w1mpro'] + 2.673
    NEOWISER_W2_AB      = NEOWISER_plot['w2mpro'] + 3.313
    
    NEOWISER_aver_plot  = NEOWISER_aver[np.where(NEOWISER_aver['cntr_01'] == jj)]
    NEOWISER_aver_W1_AB = NEOWISER_aver_plot['w1mpro_wgt'] + 2.673
    NEOWISER_aver_W2_AB = NEOWISER_aver_plot['w2mpro_wgt'] + 3.313
    
    ## NEOWISER W1/2 (AB)
    ax1.scatter(AllWISE_plot['w1mjdmean'],  AllWISE_plot_W1_AB, color='k', alpha=alpha, s=ms_big*2.2)
    ax1.scatter(AllWISE_plot['w1mjdmean'],  AllWISE_plot_W1_AB, color='r', alpha=alpha, s=ms_big)
    ax1.scatter(AllWISE_plot['w2mjdmean'],  AllWISE_plot_W2_AB, color='k', alpha=alpha, s=ms_big*2.2)
    ax1.scatter(AllWISE_plot['w2mjdmean'],  AllWISE_plot_W2_AB, color='c', alpha=alpha, s=ms_big)

    ax1.scatter(NEOWISER_plot['mjd'], NEOWISER_W1_AB, color='k', alpha=alpha, s=ms*1.8)
    ax1.scatter(NEOWISER_plot['mjd'], NEOWISER_W1_AB, color='r', alpha=alpha, s=ms, label='NEOWISE W1')
    ax1.scatter(NEOWISER_plot['mjd'], NEOWISER_W2_AB, color='k', alpha=alpha, s=ms*1.8)
    ax1.scatter(NEOWISER_plot['mjd'], NEOWISER_W2_AB, color='c', alpha=alpha, s=ms, label='NEOWISE W2')

    ax1.scatter(NEOWISER_aver_plot['mean_mjd'], NEOWISER_aver_W1_AB, color='k', alpha=alpha, s=ms_big*2.2)
    ax1.scatter(NEOWISER_aver_plot['mean_mjd'], NEOWISER_aver_W1_AB, color='r', alpha=alpha, s=ms_big)
    ax1.scatter(NEOWISER_aver_plot['mean_mjd'], NEOWISER_aver_W2_AB, color='k', alpha=alpha, s=ms_big*2.2)
    ax1.scatter(NEOWISER_aver_plot['mean_mjd'], NEOWISER_aver_W2_AB, color='c', alpha=alpha, s=ms_big)

    
    ##  P L O T T I N G    T H E      S P E C T R A
    ## full, main spectral plot
    xmin_ax2 =   2250; xmax_ax2 = 7600
    ymin_ax2 = -4.9
    ymax_ax2 = max(sdss_flux.max(), wht_flux.max()) * 1.35
    ax2.set_xlim([xmin_ax2, xmax_ax2])
    ax2.set_ylim([ymin_ax2, ymax_ax2])
    
    ## MgII
    mgii_offset = 175
    xmin_ax3 = 2800-mgii_offset; xmax_ax3 = 2800+mgii_offset
    ymin_ax3 = -2.9; ymax_ax3 = 34.
    ax3.set_xlim([xmin_ax3, xmax_ax3])
    ax3.set_ylim([ymin_ax3, ymax_ax3])
    
    ## H-beta/OIII
    Hbeta_offset = 175
    xmin_ax4 = 4862.-Hbeta_offset; xmax_ax4 = 4862.+Hbeta_offset
    ymin_ax4 = -2.9; ymax_ax4 = 39.9   
    ax4.set_xlim([xmin_ax4, xmax_ax4])
    ax4.set_ylim([ymin_ax4, ymax_ax4])
    
    ## H-alpha
    Halpha_offset = 175
    xmin_ax5 =  6564-Halpha_offset; xmax_ax5 = 6564.61+Halpha_offset
    ymin_ax5 = -2.9; ymax_ax5 = ymax_ax2
    ax5.set_xlim([xmin_ax5, xmax_ax5])
    ax5.set_ylim([ymin_ax5, ymax_ax5])


    ##  Plotting  the  spectra!!
    ax2.plot(sdss_wavelength /(1. + redshift[ii]), sdss_flux,  '-k', lw=linewidth/4.0, label='SDSS MJD '+str(sdss_mjd))
    ax2.plot( wht_wavelength /(1. + redshift[ii]),  wht_flux,  '-r', lw=linewidth/3.4, label='WHT  MJD 58662/3')
    for ll in range(len(linelist)):
        ax2.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
        label =       linelist['LineName'][ll]
        xylabel = (  (linelist['Wavelength'][ll]), ymax_ax2/1.15)
        ax2.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )

    ## ``Zooming-in'' in various lines and line complexes...
    ax3.plot(sdss_wavelength / (1+redshift[ii]), sdss_flux, '-k', lw=linewidth/4.0)
    ax3.plot( wht_wavelength / (1+redshift[ii]),  wht_flux, '-r', lw=linewidth/4.0)
    for ll in range(len(linelist)):
        ax3.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
        label   =     linelist['LineName'][ll]
        xylabel = ((  linelist['Wavelength'][ll]), ymax_ax3/1.1)
        #ax3.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )
        
    ax4.plot(sdss_wavelength / (1+redshift[ii]), sdss_flux, '-k', lw=linewidth/4.0)
    ax4.plot( wht_wavelength / (1+redshift[ii]),  wht_flux, '-r', lw=linewidth/4.0)
    for ll in range(len(linelist)):
        ax4.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
        label   =     linelist['LineName'][ll]
        xylabel =   ((linelist['Wavelength'][ll]),  ymax_ax4/1.15)
        #ax4.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )

    ax5.plot(sdss_wavelength / (1+redshift[ii]), sdss_flux, '-k', lw=linewidth/4.0)
    ax5.plot( wht_wavelength / (1+redshift[ii]),  wht_flux, '-r', lw=linewidth/4.0)
    #ax5.plot(LRIS_r_wavelength/(1+redshift), LRIS_r_flux_smoothed, '-k', lw=linewidth/4.0)
    for ll in range(len(linelist)):
        ax5.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
        label   =     linelist['LineName'][ll]
        xylabel =   ((linelist['Wavelength'][ll]), ymax_ax5/1.15)
        #ax5.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )
        
        
    ## Making all the nice plots boxes the right dimensions and placings
    #% start: automatic generated code from pylustrator
    plt.figure(1).ax_dict = {ax.get_label(): ax for ax in plt.figure(1).axes}
    import matplotlib as mpl
    plt.figure(1).axes[0].get_xaxis().get_label().set_fontsize(fontsize)
    plt.figure(1).axes[0].get_yaxis().get_label().set_fontsize(fontsize)
    plt.figure(1).axes[0].set_position([0.060000, 0.136250, 0.345, 0.823750])
    plt.figure(1).axes[0].xaxis.labelpad = 10.0
    
    plt.figure(1).axes[1].get_xaxis().get_label().set_fontsize(fontsize/1.2)
    plt.figure(1).axes[1].get_yaxis().get_label().set_fontsize(fontsize)
    plt.figure(1).axes[1].get_yaxis().get_label().set_rotation(90.0)
    x_nudge = 0.012
    plt.figure(1).axes[1].set_position([0.464541+x_nudge, 0.560924, 0.499459, 0.399076])
    plt.figure(1).axes[1].xaxis.labelpad = 323.200000
    plt.figure(1).axes[1].yaxis.labelpad = 4.720000
    plt.figure(1).axes[1].yaxis.set_label_coords(-0.045, -0.1)
    
    plt.figure(1).axes[2].get_xaxis().get_label().set_text("")
    plt.figure(1).axes[2].set_position([0.464541+x_nudge, 0.135000, 0.146836, 0.372826])
    plt.figure(1).axes[3].set_position([0.643228+x_nudge, 0.135000, 0.146836, 0.372826])
    plt.figure(1).axes[3].xaxis.labelpad = 10.0
    plt.figure(1).axes[4].set_position([0.817164+x_nudge, 0.135000, 0.146836, 0.372826])
    #% end: automatic generated code from pylustrator
    #plt.show()
    
        
    ##  A X E S   L A B E L S
    plt.rcParams['text.usetex'] = True
    ax1.set_xlabel(r'MJD',                                                          fontsize=fontsize/1.2)
    ax1.set_ylabel(r'AB Magnitude',                                                 fontsize=fontsize/1.2)
    ax2.set_xlabel(r'Rest Wavelength / ${ \rm \AA }$',                              fontsize=fontsize/1.2)
    ax2.set_ylabel(r'Flux,  $f_{\lambda}$ / 10$^{-17}$ erg/s/cm$^2$/${ \rm \AA }$', fontsize=fontsize/1.2)
    ax4.set_xlabel(r'Restframe wavelength / ${ \rm \AA }$',                         fontsize=fontsize/1.2)
    ## Title label
    #ax2.set_title(str(object_name[ii]), fontsize=fontsize*1.4)

    ## Axes style
    #ax2.tick_params(axis='both', which='major', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength,   width=tickwidth)
    #ax2.tick_params(axis='both', which='minor', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength/2, width=tickwidth)
    ax2.tick_params(labelsize=labelsize/1.25, direction='in', length=ticklength,   width=tickwidth)
    ax3.tick_params(labelsize=labelsize/1.5,  direction='in', length=ticklength,   width=tickwidth)
    ax4.tick_params(labelsize=labelsize/1.5,  direction='in', length=ticklength,   width=tickwidth)
    ax5.tick_params(labelsize=labelsize/1.5,  direction='in', length=ticklength,   width=tickwidth)

    print(ii, object_name[ii], type(object_name[ii]),  ymax_ax3)
    ## `Polishing; individual objects
    ##  First five...
    if str(object_name[ii]) == 'J1348+2456':
        ax1.axvline(x=53535, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1348+2456',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 53535', color='k')
    if str(object_name[ii]) == 'J1555+2119':
        ax1.axvline(x=53557, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1555+2119',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 53557', color='k')
    if str(object_name[ii]) == 'J1600+3345':
        ax1.axvline(x=53142, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1600+3345',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 53142', color='k')
        ymax_ax3 = 79.9
    if str(object_name[ii]) == 'J1601+4745':
        ax1.axvline(x=52354, linewidth=lw, linestyle='dashed', color='k')
        ax1.axvline(x=57895, linewidth=lw, linestyle='dashed', color='c')
        obj_label = mlines.Line2D([], [], label='J1601+4745',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 52354', color='k')
        mmt_line  = mlines.Line2D([], [], label=' MMT MJD 57895', color='c', )
        ## Plotting the MMT data from Chelsea::
        mmt_infile = 'MMT_spectra/J160111.coaddmean'
        MMT_J1601p4745_data       = ascii.read(path+mmt_infile)
        MMT_J1601p4745_flux       = MMT_J1601p4745_data['flux'] * 1e17
        MMT_J1601p4745_wavelength = MMT_J1601p4745_data['wavelength'] 
        ax2.plot((MMT_J1601p4745_wavelength / (1.+redshift[ii])), MMT_J1601p4745_flux, '-c', lw=linewidth/4.0)
        ax3.plot((MMT_J1601p4745_wavelength / (1.+redshift[ii])), MMT_J1601p4745_flux, '-c', lw=linewidth/4.0)
        ax4.plot((MMT_J1601p4745_wavelength / (1.+redshift[ii])), MMT_J1601p4745_flux, '-c', lw=linewidth/4.0)
        ax5.plot((MMT_J1601p4745_wavelength / (1.+redshift[ii])), MMT_J1601p4745_flux, '-c', lw=linewidth/4.0)
        #leg = ax2.legend(loc='lower right', fontsize=fontsize*1.2)
        #for line in leg.get_lines(): line.set_linewidth(3.6)
        ymax_ax2 = 28.
        ymin_ax2 = -2.9
        ymax_ax4 = 22.
    if str(object_name[ii]) == 'J1603+1531':
        ax1.axvline(x=53555, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1603+1531',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 53555', color='k')
        ymax_ax4 = 18.

        
    ## Second five objects...
    if str(object_name[ii]) == 'J1605+4834':
        ax1.axvline(x=52054, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1605+4834',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 52054', color='k')
        ymin_ax2 =   -24.9
    if str(object_name[ii]) == 'J1605+2309':
        ax1.axvline(x=53524, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1605+2309',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 53524', color='k')
        ymax_ax2 = 20.
    if str(object_name[ii]) == 'J1610+1525':
        ax1.axvline(x=53918, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1610+1525',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 53918', color='k')
        #leg = ax2.legend(loc='lower left', fontsize=fontsize*1.2)
        #for line in leg.get_lines(): line.set_linewidth(3.6)
    if str(object_name[ii]) == 'J1615+2507':
        ax1.axvline(x=53493, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1615+2507',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 53493', color='k')
    if str(object_name[ii]) == 'J1632+2403':
        ax1.axvline(x=53177, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1632+2403',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 53177', color='k')

    ## 'Middle four' objects
    if str(object_name[ii]) == 'J1704+3331':
        ax1.axvline(x=52427, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1704+3331',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 52427', color='k')
        ymin_ax2 = -8.
    if str(object_name[ii]) == 'J1709+3421':
        ax1.axvline(x=54586, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1709+3421',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 54586', color='k')
        ymin_ax2 = -8.
    if str(object_name[ii]) == 'J1713+2736':
        ax1.axvline(x=52410, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J1713+2736',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 52410', color='k')
        ymin_ax2 = -14.9
        
    ##  RA > 300 degs; first five
    if str(object_name[ii]) == 'J2102-0645':
        ax1.axvline(x=52174, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J2102-0645',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 52174', color='k')
    if str(object_name[ii]) == 'J2222+1733':
        ax1.axvline(x=56944, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J2222+1733',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 56944', color='k')
    if str(object_name[ii]) == 'J2223+2101':
        ax1.axvline(x=56960, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J2223+2101',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 56960', color='k')
        ymin_ax2 = -19.9
    if str(object_name[ii]) == 'J2232-0806':
        ax1.axvline(x=56960, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label=r'The BIG DIPPER!',     color='k', linestyle='None')
        #sdss_line = mlines.Line2D([], [], label=r'SDSS MJD 56960', color='k')
    if str(object_name[ii]) == 'J2251+2419':
        ax1.axvline(x=56238, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J2251+2419',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 56238', color='k')
        ymax_ax2 = 15.
        ymin_ax2 = -7.49
        
    ##  RA > 300 degs; last four
    if str(object_name[ii]) == 'J2256+1450':
        ax1.axvline(x=52263, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J2256+1450',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 52263', color='k')
        ymax_ax2 = 60.
    if str(object_name[ii]) == 'J2307+1901':
        ax1.axvline(x=56946, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J2307+1901',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 56946', color='k')
        ymax_ax2 = 120.
        ymin_ax2 = -34.9
    if str(object_name[ii]) == 'J2320+2305':
        ax1.axvline(x=56535, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J2320+2305',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 56535', color='k')
        ymin_ax2 = -14.9
    if str(object_name[ii]) == 'J2322+2235':
        ax1.axvline(x=56977, linewidth=lw, linestyle='dashed', color='k')
        obj_label = mlines.Line2D([], [], label='J2322+2235',     color='k', linestyle='None')
        sdss_line = mlines.Line2D([], [], label='SDSS MJD 56977', color='k')
        ax2.plot(sdss_wavelength/(1.+ redshift[ii]), sdss_flux,  '-k', lw=linewidth/4.0)
        #leg = ax2.legend(loc='lower right', fontsize=fontsize*1.2)
        #for line in leg.get_lines(): line.set_linewidth(3.6)
        ymin_ax2 = -34.9
        
    ## x- and y-limits for the spectra    
    ax2.set_xlim([xmin_ax2, xmax_ax2])
    ax2.set_ylim([ymin_ax2, ymax_ax2])
    ax3.set_xlim([xmin_ax3, xmax_ax3])
    ax3.set_ylim([ymin_ax3, ymax_ax3])
    ax4.set_xlim([xmin_ax4, xmax_ax4])
    ax4.set_ylim([ymin_ax4, ymax_ax4])
    ax5.set_xlim([xmin_ax5, xmax_ax5])
    ax5.set_ylim([ymin_ax5, ymax_ax5])

    ##  LIGHT CURVE Legend
    neo_w1 = mlines.Line2D([], [], label='NEOWISE W1', color='red',
                               marker="o", markeredgecolor='k', markeredgewidth=1.6, markersize=8, linestyle='None')
    neo_w2 = mlines.Line2D([], [], label='NEOWISE W2', color='c',
                               marker="o", markeredgecolor='k', markeredgewidth=1.6, markersize=8, linestyle='None')
    wht_line  = mlines.Line2D([], [], label=' WHT MJD 58662 ', color='r', )
    
    handles=[obj_label, neo_w1, neo_w2, sdss_line, wht_line]
    
    if str(object_name[ii]) == 'J1601+4745':
        handles=[obj_label, neo_w1, neo_w2, sdss_line, mmt_line, wht_line]
        
    leg = ax1.legend(loc='upper left', fontsize=fontsize/1.35, handles=handles,
                    frameon=True, framealpha=1.0, fancybox=True)

    print()    
    plt.savefig(str(object_name[ii])+'_landscape_temp.png', format='png')
    plt.close(fig)


