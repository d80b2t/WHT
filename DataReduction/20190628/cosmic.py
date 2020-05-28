from pyraf import iraf
import numpy as np

## attempts to remove cosmic rays from each image and outputs a bad
## pixel map, usually works well but images must be visually inspected
## for signs of fitting issues. In case of any problems, try reducing
## the x/y orders until the desired result is achieved.  Also,
## double-check exposure times will catch all images in the appropriate
## lists.  Using _alternative_ for the Standard stars (e.g. shorter,
## brighter objects)


def cosmic(ffile):
	iraf.lacos_spec(ffile, 'c'+ffile[1:], 'p'+ffile[1:], gain="1.16", readn="5", xorder="5", yorder="50", sigclip="10", sigfrac="0.1", objlim="1", niter="4", verbose="yes", mode="ql")

def alternative(ffile):
	iraf.lacos_spec(ffile, 'c'+ffile[1:], 'p'+ffile[1:], gain="1.16", readn="5", xorder="5", yorder="50", sigclip="100", sigfrac="0.01", objlim="1", niter="4", verbose="yes", mode="ql")

    
data = np.loadtxt('log.csv', dtype='string')

## This logic basically aims to use attributes in the log.csv file to
## idenfity which files are Biases, Flat Fields and Standard Star
## files.
## For 20190628 and 20190629 all Science Frames are 1800 secs
std = data[ np.where( (data[:,8] != '900') & (data[:,8] != '1800') & (data[:,2] != '0') & (data[:,8] != '300') & ((data[:,7] == 'ISISR') | (data[:,7] == 'ISISB')) ) ][:,0]
sci = data[ np.where( ((data[:,8] == '900') | (data[:,8] == '1800') | (data[:,8] == '300')) & ((data[:,7] == 'ISISR') | (data[:,7] == 'ISISB')) ) ][:,0]


for i in std:
	alternative('f'+i)
for i in sci:
	cosmic('f'+i)

## Can use this script in "single exposure" mode::
# cosmic('f2364346')
