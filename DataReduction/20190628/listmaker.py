import numpy as np

# produces lists from the edited runlog log.csv which is a stripped-down version with no spaces and no unnecessary images, at present different slit widths
# require different log.csv files so should be treated as a spearate reduction process. Remember to remove spaces from the runlog spreadsheet before saving as a csv
# prefixes: r = raw images, b = bias corrected, f = bias corrected & flat-fielded

def saver(filename, files, pre, post):
	new = list(files)
	for i in range(len(files)):
		new[i] = pre + new[i] + post
	np.savetxt(filename, new, fmt="%s")

data = np.loadtxt('log.csv', dtype='string')

allfiles = data[:,0]

rbiasfiles = data[np.where(data[:,1] == 'isisredbias')][:,0]
bbiasfiles = data[np.where(data[:,1] == 'isisbluebias')][:,0]

rnonbiasfiles = data[np.where( (data[:,7] == 'ISISR') & (data[:,1] != 'isisredbias') )][:,0]
bnonbiasfiles = data[np.where( (data[:,7] == 'ISISB') & (data[:,1] != 'isisbluebias') )][:,0]

rflatfiles = data[np.where(data[:,1] == 'isisredflat')][:,0]
bflatfiles = data[np.where(data[:,1] == 'isisblueflat')][:,0]

rdatafiles = data[np.where( (data[:,7] == 'ISISR') & (data[:,1] != 'isisredbias') & (data[:,1] != 'isisredflat') & (data[:,1] != 'isisredarc') )][:,0]
bdatafiles = data[np.where( (data[:,7] == 'ISISB') & (data[:,1] != 'isisbluebias') & (data[:,1] != 'isisblueflat') & (data[:,1] != 'isisbluearc') )][:,0]

saver('All.csv', allfiles, 'r', '.fit[1]')
saver('AllFits.csv', allfiles, 'r', '.fits')
saver('RBias.csv', rbiasfiles, 'r', '.fits')
saver('BBias.csv', bbiasfiles, 'r', '.fits')

saver('RNonBias.csv', rnonbiasfiles, 'r', '.fits')
saver('bRNonBias.csv', rnonbiasfiles, 'b', '.fits')
saver('BNonBias.csv', bnonbiasfiles, 'r', '.fits')
saver('bBNonBias.csv', bnonbiasfiles, 'b', '.fits')

saver('RFlat.csv', rflatfiles, 'b', '.fits')
saver('BFlat.csv', bflatfiles, 'b', '.fits')

saver('RData.csv', rdatafiles, 'b', '.fits')
saver('fRData.csv', rdatafiles, 'f', '.fits')
saver('BData.csv', bdatafiles, 'b', '.fits')
saver('fBData.csv', bdatafiles, 'f', '.fits')





