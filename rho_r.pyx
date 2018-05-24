import numpy as np
cimport numpy as np
import fitsio
from scipy.integrate import quad
from template_fit import NX
from template_fit import NY
from template_fit import NZ
from template_fit import XMIN
from template_fit import XMAX
from template_fit import YMIN
from template_fit import YMAX
from template_fit import ZMIN
from template_fit import ZMAX
from template_fit import dx
from template_fit import dy
from template_fit import dz

def integrand(double z):
	cdef double Om_m = 0.31
	cdef double Om_v = 0.69
	return 2997.92458/np.sqrt((1.+z)*(1.+z)*(1.+z)*Om_m+Om_v)

def rho_r(char* randomfilename, char* outputfilename):
	cdef np.ndarray rhor=np.zeros(NX*NY*NZ)
	print "reading Random data"
	cdef np.ndarray Randomdata = fitsio.read(randomfilename)
	print "Random data read, folding onto grid"
	cdef np.ndarray gal
	cdef double dp, x, y, z
	cdef int ix, iy, iz, ipos
	for gal in Randomdata:
		if gal['Z']>0.43 and gal['Z']<0.70:
			dp=quad(integrand,0.,gal['Z'])[0]
			x=dp*np.cos(np.radians(gal['DEC']))*np.cos(np.radians(gal['RA']))
			y=dp*np.cos(np.radians(gal['DEC']))*np.sin(np.radians(gal['RA']))
			z=dp*np.sin(np.radians(gal['DEC']))
			ix=int((x-XMIN)/dx)
			iy=int((y-YMIN)/dy)
			iz=int((z-ZMIN)/dz)
			ipos = iz+NZ*(iy+NY*ix)
			rhor[ipos]+=gal['WEIGHT_FKP']
	np.savetxt(outputfilename,rhor)
