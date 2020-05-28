from pyraf import iraf
import numpy as np

## starts the aperture extraction process for every non-calibration
## image. Set nsum to 20/100/200 and use 'z' to find the FWHM at lines
## 1000/2022/3000 take the maximum value and multiply by 3/2 to get
## the upper/lower aperture limits, set using ':upp/low'.
## 'c' will center the aperture,
## 'd' to delete and
## 'm' to mark new.
## Define ## appropriate background regions,
## 'b' to enter background mode,
## 't' to delete regions and
## 's' to define new ones,
## 'f' to fit.
## Use 'q' to leave background/aperture mode and the next screen allows you to
## fit the trace.
## Use 'w'+'e'+'e' to zoom,
## 'w'+'a' to return to full-screen.
##
## Use 's'/'t' to define trace region, avoid 'kinks'
## towards the wings, ':o <n>' sets the fit order and 'f' to fit.
## Generally a 4th order fit is all that is needed.  also extracts the
## identical aperture from the arc image to be used in the wavelength
## calibration process

## 
## for the data from 20190628, then e.g. lower and upper might be more like -3.5 and (+)3.5
##



def aperture(cfile):
	#Massey changes saturation=32400, readnoise="12", gain="1", t_order="2", nsum="10", tsum="10", t_step="10", b_naverage="-100", b_niterate="0", line="INDEF"
	iraf.apall('c'+cfile[1:], output='a'+cfile[1:], apertures="", format="multispec", references="", profiles="", interactive="yes", find="yes", recenter="yes", resize="no",
			   edit="yes", trace="yes", fittrace="yes", extract="yes", extras="yes", review="yes", line="2022", nsum="20", lower="-5.5", upper="5.5", apidtable="",
			   b_function="chebyshev", b_order="2", b_sample="-38:-8,8:38", b_naverage="-3", b_niterate="1", b_low_reject="3", b_high_reject="3", b_grow="0", width="8",
			   radius="8", threshold="0", nfind="1", minsep="5", maxsep="1000", order="increasing", aprecenter="", npeaks="INDEF", shift="yes", llimit="INDEF", ulimit="INDEF",
			   ylevel="0.5", peak="yes", bkg="yes", r_grow="0", avglimits="no", t_nsum="20", t_step="20", t_nlost="3", t_function="legendre", t_order="3", t_sample="*",
			   t_naverage="1", t_niterate="1", t_low_reject="3", t_high_reject="3", t_grow="0", background="fit", skybox="1", weights="variance", pfit="fit1d", clean="yes",
			   saturation="60000", readnoise="READNOIS", gain="GAIN", lsigma="4", usigma="4", nsubaps="1", mode="ql")

def arcaperture(arc, cfile):
	iraf.apall(arc, output='arc'+cfile[1:], apertures="", format="multispec", references='c'+cfile[1:], profiles="", interactive="no", find="yes", recenter="no", resize="no",
			   edit="yes", trace="no", fittrace="yes", extract="yes", extras="yes", review="yes", line="2022", nsum="20", lower="-5.5", upper="5.5", apidtable="",
			   b_function="chebyshev", b_order="1", b_sample="-38:-8,8:38", b_naverage="-3", b_niterate="1", b_low_reject="3", b_high_reject="3", b_grow="0", width="8", radius="8",
			   threshold="0", nfind="1", minsep="5", maxsep="1000", order="increasing", aprecenter="", npeaks="INDEF", shift="yes", llimit="INDEF", ulimit="INDEF", ylevel="0.05",
			   peak="yes", bkg="yes", r_grow="0", avglimits="no", t_nsum="20", t_step="20", t_nlost="3", t_function="legendre", t_order="3", t_sample="*", t_naverage="1",
			   t_niterate="1", t_low_reject="3", t_high_reject="3", t_grow="0", background="none", skybox="1", weights="variance", pfit="fit1d", clean="no", saturation="60000",
			   readnoise="READNOIS", gain="GAIN", lsigma="4", usigma="4", nsubaps="1", mode="ql")

data = np.loadtxt('log.csv', dtype='string')
rdatafiles = data[np.where( (data[:,7] == 'ISISR') & (data[:,1] != 'isisredbias') & (data[:,1] != 'isisredflat') & (data[:,1] != 'isisredarc') )][:,0]
bdatafiles = data[np.where( (data[:,7] == 'ISISB') & (data[:,1] != 'isisbluebias') & (data[:,1] != 'isisblueflat') & (data[:,1] != 'isisbluearc') )][:,0]
Rarcs = data[np.where(data[:,1] == 'isisredarc')][:,0]
Barcs = data[np.where(data[:,1] == 'isisbluearc')][:,0]

for i in rdatafiles[15:]:
	aperture('c'+i)
	arcaperture('b'+Rarcs[0], 'c'+i)
for i in bdatafiles:
	aperture('c'+i)
	arcaperture('b'+Barcs[0], 'c'+i)

#aperture('c2036153')
#arcaperture('b2036087', 'c2036153')
