'''
'''

from astropy.io import fits
from astropy.io import ascii
from astropy    import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, MaskedColumn

import numpy as np
import pandas as pd               

data = ascii.read("WHT_IR_CLQ_positions.tbl") 
data

len(data)

name = np.array(data['ObjectName'])
ra   = np.array(data['RA'])
dec  = np.array(data['Decl'])
print(type(ra))

d = Table([name, ra, dec], names=['name', 'ra', 'dec'])

ascii.write(d, 'WHT_IRCLQs_4IPAC.dat', format='ipac') 

