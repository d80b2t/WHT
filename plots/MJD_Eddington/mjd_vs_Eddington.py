'''

'''
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from matplotlib import colors as mcolors
from astropy.io import ascii
from astropy.io import fits

## The Shen et al. (2011) DR7Q "VAC"
path = '/cos_pc19a_npr/data/SDSS/DR7Q/'
filename ='Shen_dr7_bh_May_2010.fits'
infile = path+filename

## Okay, what we really want to do... ;-)
#data_table = pyfits.getdata(infile)
dr7q       = fits.open(infile)
data_table = dr7q[1].data

## Quick check on the format/dimensions of the FITS table file...
print(type(data_table), '\n')
print('The number of rows of is.... ', data_table.shape, '\n')
print('The number of columns is...  ', len(data_table.names), '\n\n')

## Now getting into it some more...

mjd_full = data_table.field('MJD')
logtau_full  = data_table.field('LOGEDD_RATIO')

data = data_table[np.where(data_table['LOGEDD_RATIO'] > -9.0)]
mjd    = data['MJD']
logtau = data['LOGEDD_RATIO']
tau    = (10**data['LOGEDD_RATIO'])*100


## The data for our interesting NEOWISE-R IR LC quasars
path = '/Users/npr1/Dropbox/WHT-CLQ/2019A-proposal/MJD_Eddington/'
risers_data = ascii.read(path+'name_mjd_eddingtons_risers.dat')
faders_data = ascii.read(path+'name_mjd_eddingtons_faders.dat')


##
## Making the plot
##
fig, ax = plt.subplots(figsize=(14.0, 8.0))
left   = 0.08   # the left side of the subplots of the figure
right  = 0.96   # the right side of the subplots of the figure
bottom = 0.12   # the bottom of the subplots of the figure
top    = 0.96   # the top of the subplots of the figure
wspace = 0.22   # the amount of width reserved for blank space between subplots
hspace = 0.06   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)


lw       = 5.0
ls       = 'solid'
ms       = 80.
ms_large = ms*3.
fontsize = 24
alpha    = 1.00
totes    = 0
gridsize  = 30  #use log bins
mincnt    = .1
#color_map = plt.cm.get_cmap('viridis')
color_map = plt.cm.get_cmap('Greys')

ax.hexbin(mjd, logtau,         bins='log',  gridsize=gridsize, cmap=color_map, mincnt=mincnt)
ax.scatter(faders_data['mjd'], np.log10(faders_data['EDD_RATIO']/100),  marker='o', s=ms_large*1.6, color='black')
ax.scatter(faders_data['mjd'], np.log10(faders_data['EDD_RATIO']/100),  marker='o', s=ms_large,     color='dodgerblue')
ax.scatter(risers_data['mjd'], np.log10(risers_data['EDD_RATIO']/100),  marker='o', s=ms_large*1.6, color='black')
ax.scatter(risers_data['mjd'], np.log10(risers_data['EDD_RATIO']/100),  marker='o', s=ms_large, color='red')

## Tidy up the figure
xmin     = 51400      ## 51100
xmax     = 58664.+360
ymin     =   -4.1     ## -3.4 
ymax     =   1.4      ##  1.2
ymin_log =   0.1      ##  0.1
ymax_log = 100.0    

ax.set_xlim((xmin, xmax))
#ax.set_xlim((54000, xmax))
ax.set_ylim((ymin, ymax))
#ax.set_yscale('log')
#ax.set_ylim((ymin_log, ymax_log))

ax.tick_params('x', direction='in')
ax.tick_params('y', direction='in')
#ax.minorticks_on('x', direction='in')
#ay.minorticks_on()
ax.tick_params('x', direction='in', which='major', bottom='True', top='True', left='True', right='True', labelsize=fontsize/1.2)
ax.tick_params('x', direction='in', which='minor', bottom='True', top='True', left='True', right='True', labelsize=fontsize/1.2)
ax.tick_params('y', direction='in', which='both', bottom='True', top='True', left='True', right='True', labelsize=fontsize)
ax.minorticks_on()

## ALLWISE   timespan
ALLWISE_MJD_min = 55210.     #  2010-January-14 
ALLWISE_MJD_max = 55593.      
ax.axvspan(ALLWISE_MJD_min, ALLWISE_MJD_max, alpha=0.9, color='red')

## NEOWISE-R timespan
NEOWISER_MJD_min = 56639.
NEOWISER_MJD_max = 58100.   ## 2017-Dec-13
NEOWISER_MJD_max = 58465.   ## 2018-Dec-13
NEOWISER_MJD_max = 58484.   ## start of 2019A
ax.axvspan(NEOWISER_MJD_min, NEOWISER_MJD_max, alpha=0.5, color='red')

## 2019A
Semester_2019A_MJD_min = 58484.0
Semester_2019A_MJD_max = 58664.0
ax.axvspan(Semester_2019A_MJD_min, Semester_2019A_MJD_max, alpha=0.5, color='blue')

ax.text( ALLWISE_MJD_min,          0.8, 'ALLWISE',   style='italic',    fontsize=fontsize/1.2, rotation=270)
ax.text(NEOWISER_MJD_min,          0.8, 'NEOWISE-R', style='italic',    fontsize=fontsize/1.2, rotation=270)
#ax.text(NEOWISER_MJD_min,          0.8, 'NEOWISE-R', style='italic',    fontsize=fontsize/1.4)
ax.text(Semester_2019A_MJD_min-20, 0.8, '2019A',     fontweight='bold', fontsize=fontsize/1.4, rotation=270)  

ax.set_xlabel('MJD',                       fontsize=fontsize)
ax.set_ylabel(r'log$_{10}$ Eddington Ratio', fontsize=fontsize)

##plt.show()
plt.savefig('MJD_vs_Eddington_temp.pdf', format='pdf')
#plt.savefig('bias_with_redshift_temp.png',format='png')

#plt.show()
plt.close(fig)

