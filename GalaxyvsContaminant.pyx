'''
    Programme to compute the ratio of observed galaxies to expected
    galaxies in contaminant bins
    Copyright (C) 2018-2019 Benedict Kalus
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import numpy as np
cimport numpy as np
import healpy as hp
import fitsio
from cpython cimport bool
from cython.parallel import prange

DTYPE = np.int
ctypedef np.int_t DTYPE_t

# Find HealPix index corresponding to Equatorial coordinates
cdef int DeclRaToIndex(double decl,double RA) nogil:
	with gil:
		return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(RA))

cdef int NSIDE = 64

def DeclRaToIndex_py(double decl,double RA):
	return DeclRaToIndex(decl,RA)

# Contaminant bin corresponding to contamination level f
cpdef int binid(double f, np.ndarray binning):
	return (np.abs(f-binning)).argmin()

#main routine, fsys is a healpix map of the contaminant, mask is an auxiliary map that allows masking out pixels that are not part of the survey or shouldn't be considered for other reasons, minfsys and maxfsys correspond to the contaminant values of the minimum and maximum bin, Galaxyfiles is a list containing the filenames of galaxy catalogue FITS files (can be more than one to allow e.g. for data on different hemispheres), randoms provides the filenames for the corresponding random catalogues, redshiftmin and redshiftmax are the minimal and maximal redshift of the survey
cpdef np.ndarray main(np.ndarray fsys, np.ndarray mask, int numbins, np.ndarray Galaxyfiles, np.ndarray randoms, double redshiftmin=0.43, double redshiftmax=0.7, bint flag_noz=True, bint flag_cp=True, bint flag_fkp=True, bint flag_seeing=True):
	cdef int cell
	cdef np.ndarray nsys=np.zeros((numbins,len(Galaxyfiles)))
	cdef np.ndarray nbarbynsys=np.zeros((numbins,len(Galaxyfiles)))
	cdef np.ndarray sigmanbarbynsys=np.zeros((numbins,len(Galaxyfiles)))
	cdef np.ndarray ngal=np.zeros(len(Galaxyfiles))

	cdef int i, nonzerocells, bin, igal, len_Galaxydata, pos
	cdef double ng, nr, nbarsum, nbarrsum, alpha
	cdef np.ndarray Galaxydata, Randomdata, nbar, nbarr, nbarb, nbarrb, numfsys
	cdef double[:] nbar_view, nbarr_view
	cdef np.ndarray fsysnums
	cdef np.ndarray binning
	cdef double galweight
	cdef double ra,dec,red,weight_noz,weight_cp,weight_fkp,weight_seeing
	cdef bint in_mask

	for i in np.arange(len(Galaxyfiles)):
		ng=0
		Galaxydata = fitsio.read(Galaxyfiles[i],columns=['RA','DEC','Z','WEIGHT_NOZ','WEIGHT_CP','WEIGHT_FKP','WEIGHT_SEEING'])
		print "opened "+Galaxyfiles[i]
		nbar=np.zeros(hp.nside2npix(NSIDE))
		nbarr=np.zeros(hp.nside2npix(NSIDE))
		nbarb=np.zeros(numbins)
		nbarrb=np.zeros(numbins)
		nonzerocells=0
		nr=0
		len_Galaxydata=len(Galaxydata)
		nbar_view=nbar
		for igal in prange(len_Galaxydata,nogil=True):
			with gil:
				ra=Galaxydata[igal][0]
				dec=Galaxydata[igal][1]
				red=Galaxydata[igal][2]
				weight_noz=Galaxydata[igal][3]
				weight_cp=Galaxydata[igal][4]
				weight_fkp=Galaxydata[igal][5]
				weight_seeing=Galaxydata[igal][6]
			pos=int(DeclRaToIndex(dec,ra))
			with gil:
				in_mask=(mask[pos])
			if red>redshiftmin and red<redshiftmax and in_mask:
				galweight=((flag_noz*weight_noz+1.-flag_noz)+(flag_cp*weight_cp+1.-flag_cp)-1.)*(flag_fkp*weight_fkp+1.-flag_fkp)*(flag_seeing*weight_seeing+1.-flag_seeing)
				nbar_view[pos]+=galweight
				ng+=galweight
		nbarsum=0
		for cell in np.arange(0,hp.nside2npix(NSIDE)):
			nbarsum+=nbar[cell]
		print "ng:",ng,"sum of galaxies in cells:",nbarsum
		print "read "+Galaxyfiles[i]
		ngal[i]=Galaxydata.shape[0]
		del Galaxydata
		Randomdata = fitsio.read(randoms[i],columns=['RA','DEC','Z','WEIGHT_FKP'])
		print "opened "+randoms[i]
		print(Randomdata)
		len_Galaxydata=len(Randomdata)
		nbarr_view=nbarr
		for igal in prange(len_Galaxydata,nogil=True):
			with gil:
				ra=Randomdata[igal][0]
				dec=Randomdata[igal][1]
				red=Randomdata[igal][2]
				weight_fkp=Randomdata[igal][3]
			pos=int(DeclRaToIndex(dec,ra))
			with gil:
				in_mask=(mask[pos])
			if red>redshiftmin and red<redshiftmax and in_mask:
				nbarr_view[int(DeclRaToIndex(dec,ra))]+=(flag_fkp*weight_fkp+1.-flag_fkp)
				nr+=weight_fkp
		nbarrsum=0
		for cell in np.arange(0,hp.nside2npix(NSIDE)):
			nbarrsum+=nbarr[cell]
		print "nr:",nr,"sum of randoms in cells:",nbarrsum
		print "read "+randoms[i]
		numfsys=np.zeros(numbins)
		alpha=nr/ng
		del Randomdata
		print fsys[(fsys!=0) & (fsys>-1e30) & (nbarr!=0)]
		binning=np.percentile(fsys[(fsys!=0) & (fsys>-1e30) & (nbarr!=0)],np.linspace(0,100,numbins))
		print "Binning:",binning
		print "ratio of randoms per galaxy:",alpha
		for cell in np.arange(0,hp.nside2npix(NSIDE)):
			if fsys[cell]>-1e30 and not np.isnan(fsys[cell]) and binid(fsys[cell],binning)<numbins:
				nbarb[binid(fsys[cell],binning)]+=nbar[cell]
				nbarrb[binid(fsys[cell],binning)]+=nbarr[cell]
				nsys[binid(fsys[cell],binning),i]+=fsys[cell]
				numfsys[binid(fsys[cell],binning)]+=1
		for bin in range(0,numbins):
			if numfsys[bin]>0 and nbarb[bin]>0:
				nbarbynsys[bin,i]=alpha*nbarb[bin]/nbarrb[bin]
				sigmanbarbynsys[bin,i]=np.sqrt(nbarb[bin])/nbarrb[bin]*alpha
				nsys[bin,i]/=numfsys[bin]
	cdef np.ndarray nsystot=np.zeros(numbins)
	cdef np.ndarray nbarbynsystot=np.zeros(numbins)
	cdef np.ndarray sigmanbarbynsystot=np.zeros(numbins)
	cdef double N=0
	for i in np.arange(len(Galaxyfiles)):
		nsystot+=nsys[:,i]*ngal[i]
		nbarbynsystot+=nbarbynsys[:,i]*ngal[i]
		sigmanbarbynsystot+=sigmanbarbynsys[:,i]*ngal[i]
		N+=ngal[i]
		print ngal[i],nsys[1,i],nsystot[1]/N,nbarbynsys[1,i],nbarbynsystot[1]/N,sigmanbarbynsys[1,i],sigmanbarbynsystot[1]/N
	nsystot/=N
	nbarbynsystot/=N
	sigmanbarbynsystot/=N
	return np.transpose([nsystot,nbarbynsystot,sigmanbarbynsystot])
