'''
http://blog.marmakoide.org/?p=94
'''

from astropy.io import ascii
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches

import numpy as np
import random

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

plot_the_avergs = 'y' 
## averages of the e.g. NEOWISE-R L1b light curve points
if plot_the_avergs == 'y':
    file  = 'WHT_CLQs_NEOWISER_LCs_averaged.tbl'
    NEOWISER_aver = ascii.read(path+file)


##pdf_pages = PdfPages('my-fancy-document.pdf')
#pdfs       = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
out_pdf     = 'WHT_CLQs_NEOWISER_LCs_temp.pdf'
pdf_pages   = PdfPages(out_pdf)

## http://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec1_2.html#phases
##
AllWISE_MJD_min = 55203.     #  2010-January-07 
AllWISE_MJD_max = 55593.     #  2011-February-01
mjd_range_AllWISE=[AllWISE_MJD_min, AllWISE_MJD_max]

## MJDs in NEOWISE-R data:
##    NEOWISE 2015    December 13, 2013 and December 13, 2014 UTC
##    NEOWISE 2016    December 13, 2014 and December 13, 2015 UTC
##    NEOWISE 2017    December 13, 2015 and December 13, 2016 UTC
##    NEOWISE 2018    December 13, 2016 and December 13, 2017 UTC
##    NEOWISE 2019    December 13, 2017 and December 13, 2018 UTC
## from NEOWISER-R_SingleExposure_L1bs.tbl:
##    data['mjd'].min()   56639.790   which is 2013-Dec-13
##    data['mjd'].max()   58465.286   which is 2018-Dec-13
## 58654 is 2019-Jun-20

##  Reliability/Completeness from Assef et al. (2018), ApJS, 234, 23
## 90% Reliability  for:: (alpah_R90, beta_R90, gamma_R90) = (0.650, 0.153, 13.86)
## 75% Reliability  for:: (alpah_R75, beta_R75, gamma_R75) = (0.486, 0.092, 13.07)
## 90% Completeness for:: delta_C90 = 0.50
## 75% completeness for:: delta_C75 = 0.71

alpha_R   = 0.650
beta_R    = 0.153
gamma_R   = 13.86
delta_C90 = 0.50
w2_stepsize = 0.02
w2_array_hi = (np.arange(100) * w2_stepsize ) + gamma_R
w2_array_lo = (np.arange(70) * (w2_stepsize*2) ) + 10.9

## eqn. (4) from  Assef et al. (2018)
w1w2limit_hi = alpha_R * (np.exp(beta_R*((w2_array_hi - gamma_R)**2)))  
w1w2limit_lo = w2_array_lo


## define the colormap
cmap     = plt.cm.inferno_r
ls       = 'solid'
lw       = 1.0
ms       = 30.
ms_big   = ms*6.
ms_large = ms*8.
fontsize = 28
alpha    = 1.00
totes    = 0
lw       = 2.0

for x in range(upto_this):
    x=x+1
    AllWISE_plot       = AllWISE[np.where(AllWISE['cntr_01'] == x)]
    NEOWISER_plot      = NEOWISER[np.where(NEOWISER['cntr_01'] == x)]
    if plot_the_avergs == 'y': NEOWISER_aver_plot = NEOWISER_aver[np.where(NEOWISER_aver['cntr_01'] == x)]

    if len(AllWISE_plot) < 1:
          print('0')
          
    #if len(data_one) > 0:  
        #print(x, data_one['ra'][0], data_one['dec'][0], len(data_one))
          
    #if len(AllWISE_plot) > 0:
    if len(AllWISE_plot) > 0 and (len(NEOWISER_plot)>50):  
        print(x, AllWISE_plot['ra'][0], AllWISE_plot['dec'][0], 'len(NEOWISER_plot):', len(NEOWISER_plot))
#        print(len(data_one))
        totes = totes + len(AllWISE_plot)
        ra = AllWISE_plot['ra'][0]
        dec = AllWISE_plot['dec'][0]

        ## Setting up a 3 panel, "triptych" style plot
        #fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(21.5, 8.5)) # inches
        gridspec = dict(wspace=0.0, width_ratios=[3, 0.65, 1.2])
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(21.5, 8.5), gridspec_kw = gridspec ) # inches
        
        left   = 0.10   # the left side of the subplots of the figure
        right  = 0.98   # the right side of the subplots of the figure
        bottom = 0.16   # the bottom of the subplots of the figure
        top    = 0.92   # the top of the subplots of the figure
        wspace = 0.38   # the amount of width reserved for blank space between subplots
        hspace = 0.06   # the amount of height reserved for white space between subplots
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        AllWISE_MJD_mid = ((AllWISE_MJD_max + AllWISE_MJD_min))/2.
        
        ## W1/2 Light-curves, y-axis reverse
        ax1.scatter(AllWISE_plot['w1mjdmean'],  AllWISE_plot['w1mpro'], color='k', alpha=alpha, s=ms_large*1.8)
        ax1.scatter(AllWISE_plot['w1mjdmean'],  AllWISE_plot['w1mpro'], color='r', alpha=alpha, s=ms_large)
        ax1.scatter(AllWISE_plot['w2mjdmean'],  AllWISE_plot['w2mpro'], color='k', alpha=alpha, s=ms_large*1.8)
        ax1.scatter(AllWISE_plot['w2mjdmean'],  AllWISE_plot['w2mpro'], color='c', alpha=alpha, s=ms_large)
        
        ax1.scatter(NEOWISER_plot['mjd'],  NEOWISER_plot['w1mpro'], color='k', alpha=alpha, s=ms*1.8)
        ax1.scatter(NEOWISER_plot['mjd'],  NEOWISER_plot['w1mpro'], color='r', alpha=alpha, s=ms)
        ax1.scatter(NEOWISER_plot['mjd'],  NEOWISER_plot['w2mpro'], color='k', alpha=alpha, s=ms*1.8)
        ax1.scatter(NEOWISER_plot['mjd'],  NEOWISER_plot['w2mpro'], color='c', alpha=alpha, s=ms)
        
        if plot_the_avergs == 'y':
            ax1.scatter(NEOWISER_aver_plot['mean_mjd'],  NEOWISER_aver_plot['w1mpro_wgt'], color='k', alpha=alpha, s=ms_big*1.8)
            ax1.scatter(NEOWISER_aver_plot['mean_mjd'],  NEOWISER_aver_plot['w1mpro_wgt'], color='r', alpha=alpha, s=ms_big)
            ax1.scatter(NEOWISER_aver_plot['mean_mjd'],  NEOWISER_aver_plot['w2mpro_wgt'], color='k', alpha=alpha, s=ms_big*1.8)
            ax1.scatter(NEOWISER_aver_plot['mean_mjd'],  NEOWISER_aver_plot['w2mpro_wgt'], color='c', alpha=alpha, s=ms_big)

        ax1.set_xlim((55000,58600))
        ax1.set_ylim(ax1.get_ylim()[::-1])
        ax1.set_xlabel('MJD',                   fontsize=fontsize)
        ax1.set_ylabel('WISE  W1/2  magnitude', fontsize=fontsize)
        ax1.tick_params('x', direction='in', which='both', bottom='True', top='True', left='True', right='True', labelsize=fontsize)
        ax1.tick_params('y', direction='in', which='both', bottom='True', top='True', left='True', right='True', labelsize=fontsize)
        ax1.set_title('(R.A., Decl.) = {} {} {}'.format(ra,',',dec), fontsize=fontsize)

        ## Make a gap between the LC and color plots
        ax2.set_visible(False)

        ## W2 vs. W1-W2 color; to look for color-changes and cf. Assef (2013) Reliability/completeness criteria
        ax3.scatter(AllWISE_plot['w2mpro'],  (AllWISE_plot['w1mpro']-AllWISE_plot['w2mpro']),   color='k', alpha=alpha, s=ms_large*1.8)
        ax3.scatter(AllWISE_plot['w2mpro'],  (AllWISE_plot['w1mpro']-AllWISE_plot['w2mpro']),   color='c', alpha=alpha, s=ms_large*1.8)
        ax3.scatter(NEOWISER_plot['w2mpro'], (NEOWISER_plot['w1mpro']-NEOWISER_plot['w2mpro']), color='k', alpha=alpha, s=ms*1.8)
        ax3.scatter(NEOWISER_plot['w2mpro'], (NEOWISER_plot['w1mpro']-NEOWISER_plot['w2mpro']), color='c', alpha=alpha, s=ms)
        ## Assef et al. (2018) 
        ax3.plot(w2_array_hi, w1w2limit_hi, color='gray', linestyle='--', linewidth=lw/3.4)
        ax3.plot(w2_array_lo, w1w2limit_lo, color='gray', linestyle='--', linewidth=lw/3.4)
        ax3.hlines(alpha_R, 10.8, 13.86,    color='gray', linestyle='--', linewidth=lw/3.4)
        #ax3.hlines(delta_C90, 11.5, 14.2, linestyle='--', linewidth=lw/3.4)

        ax3.set_xlabel('W2 magnitude', fontsize=fontsize)
        ax3.set_ylabel(r'W1 $-$ W2',      fontsize=fontsize)
        ax3.tick_params('x', direction='in', which='both', bottom='True', top='True', left='True', right='True',  labelsize=fontsize)
        ax3.tick_params('y', direction='in', which='both', bottom='True', top='True', left='True', right='False', labelsize=fontsize)
        ax3.set_xlim(10.8, 14.4)
        ax3.set_ylim(-0.3, 1.55)

        ## Done with the page
        pdf_pages.savefig(fig)
    
## Write the PDF document to the disk
pdf_pages.close()

plt.close()







