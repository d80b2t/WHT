"""
This is a script to make several plots of the WHT IR 'CLQ' data from
2019June28 and 2019June29.
"""

import numpy  as np
import random
import datetime
import matplotlib.pyplot  as plt
import matplotlib.patches as patches

from astropy.io import ascii
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker               import ScalarFormatter


path    = '/cos_pc19a_npr/programs/quasars/WHT/data/spectra/'
#infile  = 'spectra_20190628.lis'
infile = 'spectra.lis'
objects = ascii.read(path+infile)

object_name     = objects['ObjectName']
object_redshift = objects['Redshift']

## Emission Line list
linelist_file = 'emission_lines.dat'
linelist = ascii.read(path+linelist_file)


## Create the PdfPages object to which we will save the pages:
## The with statement makes sure that the PdfPages object is closed properly at
## the end of the block, even if an Exception occurs.
##pdf_pages = PdfPages('my-fancy-document.pdf')
#pdfs       = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
out_pdf     = 'SDSS_WHT_IR_CLQs_temp.pdf'
pdf_pages   = PdfPages(out_pdf)


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
        
    ## Setting up the plot
    fig, ax = plt.subplots(figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')
    #fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
    #                               gridspec_kw={'height_ratios': [2, 1]}, 
    #                               figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')
    
    ## Adjusting the Whitespace for the plots
    left   = 0.12   # the left side of the subplots of the figure
    right  = 0.96   # the right side of the subplots of the figure
    bottom = 0.16   # the bottom of the subplots of the figure
    top    = 0.92   # the top of the subplots of the figure
    wspace = 0.16   # the amount of width reserved for blank space between subplots
    hspace = 0.06   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    ## Some NPR defaults
    alpha           = 1.0
    fontsize        = 22
    labelsize       = fontsize
    tickwidth       = 2.0
    linewidth       = 2.4
    tickwidth       = 2.0
    ticklength      = 6.0
    ticklabelsize   = labelsize
    majorticklength = 12
    minorticklength = 6
    
    ## Axes limits
    #xmin =  3000.  ## Angstroms  OBSERVED
    #xmax = 10000.  ## Angstroms  OBSERVED
    xmin =  2250.  ## Angstroms  REST
    xmax =  7600.  ## Angstroms  REST
    ymin =   -9.99  ## f_lambda  /  10^-17  erg/s/cm^2/Ang  
    ymax = max(sdss_flux.max(), wht_flux.max()) * 1.05
      
    #  Plotting  the  spectra!!
    ax.plot(sdss_wavelength/(1.+object_redshift[ii]), sdss_flux,  '-k', lw=linewidth/4.0, label='SDSS MJD '+str(sdss_mjd))
    ax.plot( wht_wavelength/(1.+object_redshift[ii]),  wht_flux,  '-r', lw=linewidth/3.4, label='WHT  MJD 58662/3')
    #ax.plot( wavelength,      wht_flux,  '-r', lw=linewidth/3.4, )

    
    ## Axes labels
    plt.rcParams['text.usetex'] = True
    ax.set_xlabel(r'Rest Wavelength / ${ \rm \AA }$',                                   fontsize=fontsize)
    ax.set_ylabel(r'Flux,  $f_{\lambda}$ / 10$^{-17}$ erg/s/cm$^2$/${ \rm \AA }$', fontsize=fontsize)
    ## Title label
    ax.set_title(str(object_name[ii]), fontsize=fontsize*1.4)
    
    ## Axes style
    ax.tick_params(axis='both', which='major', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength,   width=tickwidth)
    ax.tick_params(axis='both', which='minor', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength/2, width=tickwidth)

    ## Setting the legend linewidths
    leg = ax.legend(loc='lower left', fontsize=fontsize*1.2)
    for line in leg.get_lines(): line.set_linewidth(3.6)

    ## Polishing individual objects...
    if str(object_name[ii]) == 'J1555+2119':
        ymin =  -19.9
    if str(object_name[ii]) == 'J1601+4745':

        ## plot the MMT data from Chelsea::
        mmt_infile = 'MMT_spectra/J160111.coaddmean'
        MMT_J1601p4745_data       = ascii.read(path+mmt_infile)
        MMT_J1601p4745_flux       = MMT_J1601p4745_data['flux'] * 1e17
        MMT_J1601p4745_wavelength = MMT_J1601p4745_data['wavelength'] 
        ax.plot((MMT_J1601p4745_wavelength/(1.+object_redshift[ii])), MMT_J1601p4745_flux,
                '-c', lw=linewidth/4.0, label='MMT MJD 57895')
        
        leg = ax.legend(loc='lower right', fontsize=fontsize*1.2)
        for line in leg.get_lines(): line.set_linewidth(3.6)
        ymax = 28.
        ymin = -9.9
            
    if str(object_name[ii]) == 'J1605+2309':
        ymax = 20.
    if str(object_name[ii]) == 'J1605+4834':
        ymin =   -24.9
    if str(object_name[ii]) == 'J1610+1525':
        leg = ax.legend(loc='lower left', fontsize=fontsize*1.2)
        for line in leg.get_lines(): line.set_linewidth(3.6)
    if str(object_name[ii]) == 'J1709+3421':
        ymin = -8.
    if str(object_name[ii]) == 'J1713+2736':
        ymin = -14.9
    if str(object_name[ii]) == 'J2223+2101':
        ymin = -19.9
    if str(object_name[ii]) == 'J2251+2419':
        ymax = 15.
        ymin = -7.49
    if str(object_name[ii]) == 'J2256+1450':
        ymax = 60.
    if str(object_name[ii]) == 'J2307+1901':
        ymax = 120.
        ymin = -34.9
    if str(object_name[ii]) == 'J2320+2305':
        ymin = -14.9
    if str(object_name[ii]) == 'J2322+2235':
        ax.plot(sdss_wavelength/(1.+object_redshift[ii]), sdss_flux,  '-k', lw=linewidth/4.0)
        leg = ax.legend(loc='lower right', fontsize=fontsize*1.2)
        for line in leg.get_lines(): line.set_linewidth(3.6)
        ymin = -34.9

    #ax.plot( wht_wavelength,  wht_flux,  '-r', lw=linewidth/3.4, label='WHT  MJD 58662')

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin, ymax])

    ## Putting on the Emission Line names and vertical lines
    ## e.g. https://scipython.com/blog/rotating-text-onto-a-line-in-matplotlib/
    for ll in range(len(linelist)):
        plt.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
        label = linelist['LineName'][ll]
        #print(label)
        label_tex = '$'+label+'$'
        xylabel = ((linelist['Wavelength'][ll]), (ymax/1.10))
        ax.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )

                
    ## Done with the page
    pdf_pages.savefig(fig)

## Write the PDF document to the disk
pdf_pages.close()
  
plt.close()
        
    

'''    
def getwave(image, length):
	iraf.imgets(image, "CRVAL1")
	start = iraf.imgets.value
	iraf.imgets(image, "CD1_1")
	delta = iraf.imgets.value
	iraf.imgets(image, "CRPIX1")
	offset = iraf.imgets.value
	start, delta, offset = float(start), float(delta), float(offset)
	wavelength = []
	for i in range(length):
		wavelength += [ start + delta * ( i + 1 - offset ) ]
	return(wavelength)

'''
