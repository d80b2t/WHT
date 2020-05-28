from pyraf import iraf

# performs the initial reduction steps: bias, trim and flat-fielding.  First creates .fits images of all files and then asks you to check overscan regions.
# a 5th order fit seems to work well, check a few images with 'q'+<enter> and then use 'q'+'NO'+<enter> to apply the same to the rest. Next is the flat-field
# normalisation where an ':o 200' spline3 works well, fit with 'f' and 'q' out if happy. Perform the same steps for the B images.
# if using narrow window, comment out badpixmasks and fixfile=""; fixpix="no"

def wordy(output):
	print
	print
	print
	print('***   '+output+'   ***')
	print

def getfits():
	wordy('Converting all to .fits')
	iraf.imcopy('@All.csv', '@AllFits.csv')

def getbias(c):
	wordy('Getting Master Bias in '+c)
	#Massey had ccdtype="zero", nkeep=NA, rdnoise=0, gain=1
	iraf.zerocombine('@'+c+'Bias.csv', output=c+'Bias', combine="average", reject="minmax", ccdtype="", process="no", delete="no", clobber="no", scale ="none", statsec="", nlow="0",
					 nhigh="1", nkeep="1", mclip="yes", lsigma="3", hsigma="3", rdnoise="READNOIS", gain="GAIN", snoise="0", pclip="-0.5", blank="0", mode="ql")
	iraf.imstat('@'+c+'Bias.csv')
	iraf.imstat(c+'Bias')
	
def getflat(c):
	wordy('Getting Master Flat in '+c)
	#Massey has ccdtype="flat", process="yes", subsets="yes", nkeep=NA, rdnoise="rdnoise", gain="gain"
	iraf.flatcombine('@'+c+'Flat.csv', output=c+'Flat', combine="average", reject="crreject", ccdtype="", process="no", subsets="no", delete="no", clobber="no", scale ="mode",
					 statsec="", nlow="1", nhigh="1", nkeep="1", mclip="yes", lsigma="3", hsigma="3", rdnoise="READNOIS", gain="GAIN", snoise="0", pclip="-0.5", blank="1", mode="ql")
	iraf.imstat('@'+c+'Flat.csv')
	iraf.imstat(c+'Flat')
	
def biasandtrim(c):
	wordy('Performing trim and bias correction in '+c)
	if c == 'R':
		bsec, tsec = '[20:483,4110:4190]', '[20:483,5:4096]'
		fixer = 'Rbadpix_stdwin_21'
	if c == 'B':
		bsec, tsec = '[20:483,4110:4190]', '[20:483,1:4050]'
		fixer = 'Bbadpix_stdwin_21'
	#Massey has readaxis="line", biassec="image", trimsec="[1:3064, 10:101]", zero="Zeron3", minreplace="1", order="1"
	iraf.ccdproc('@'+c+'NonBias.csv', output='@b'+c+'NonBias.csv', ccdtype="", max_cache="0", noproc="no", fixpix="yes", overscan="yes", trim="yes", zerocor="yes", darkcor="no",
				 flatcor="no", illumcor="no", fringecor="no", readcor="no", scancor="no", readaxis="column", fixfile=fixer, biassec=bsec, trimsec=tsec, zero=c+'Bias', dark="", flat="",
				 illum="", fringe="", minreplace="-1000", scantype="shortscan", nscan="1", interactive="yes", function="chebyshev", order="5", sample="*", naverage="1", niterate="1",
				 low_reject="3", high_reject="3", grow="1", mode="ql")
	iraf.imstat('@'+c+'NonBias.csv')
	iraf.imstat('@b'+c+'NonBias.csv')

def normflat(c):
	wordy('Normalising flat in '+c)
	iraf.response(c+'Flat', c+'Flat', 'n'+c+'Flat', interactive="yes", threshold="INDEF", sample="*", naverage="1", function="spline3", order="1", low_reject="3", high_reject="3",
				  niterate="1", grow="0", graphics="stdgraph", cursor="", mode="ql")
	iraf.imstat('n'+c+'Flat')

def flatfield(c):
	wordy('Performing flat correction in '+c)
	#Massey has overscan="yes", trim="yes", zerocor="yes", readaxis="line", biassec="image", trimsec="[1:3064, 10:101]", zero="Zeron3", flat="nFlat0", order="1"
	iraf.ccdproc('@'+c+'Data.csv', output='@f'+c+'Data.csv', ccdtype="", max_cache="0", noproc="no", fixpix="no", overscan="no", trim="no", zerocor="no", darkcor="no",
				 flatcor="yes", illumcor="no", fringecor="no", readcor="no", scancor="no", readaxis="column", fixfile="", biassec="", trimsec="", zero="", dark="", flat='n'+c+'Flat',
				 illum="", fringe="", minreplace="-1000", scantype="shortscan", nscan="1", interactive="yes", function="chebyshev", order="5", sample="*", naverage="1", niterate="1",
				 low_reject="3", high_reject="3", grow="1", mode="ql")
	iraf.imstat('@'+c+'Data.csv')
	iraf.imstat('@f'+c+'Data.csv')

getfits()
process = [ 'R' , 'B' ]
for i in process:
	getbias(i)
	biasandtrim(i)
	getflat(i)
	normflat(i)
	flatfield(i)
wordy('End of line...')
