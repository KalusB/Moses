import numpy as np
import fitsio
import sys

if not len(sys.argv)==5:
	print("Error: Wrong number of arguments. Call fitstolist.py infile weight zmin zmax")
	sys.exit()
if sys.argv[2] not in ['WEIGHT_SEEING','WEIGHT_SYSTOT','WEIGHT_STAR','None']:
	print("Error: Weight has to be WEIGHT_SEEING, WEIGHT_SYSTOT, WEIGHT_STAR or None")
	sys.exit()
Galaxydata = fitsio.read(sys.argv[1])
fieldid=[]
run=[]
if sys.argv[2]=='None':
	for gal in Galaxydata:
		if gal['Z']>float(sys.argv[3]) and gal['Z']<float(sys.argv[4]):
			print(gal['RA'],gal['DEC'],gal['Z'],gal['WEIGHT_FKP'],gal['WEIGHT_CP'],gal['WEIGHT_NOZ'])
			#print gal['RA'],gal['DEC'],gal['Z'],gal['WEIGHT_FKP']
			#print(gal['RA'],gal['DEC'],gal['Z'],gal['WEIGHT_FKP'],gal['WEIGHT_CP'],gal['WEIGHT_NOZ'],gal[sys.argv[2]])
