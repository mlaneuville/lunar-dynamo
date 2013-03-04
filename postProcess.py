import numpy as np
import os

dir = '../out/'
for file in os.listdir(dir):
	str = []
	if file.find('dat') != -1:
		str.append(file)
		spl = file.split('_')
		spl[-1] = str[-1][:-4]
		str.append(spl[0])
		str.append(spl[1])
		str.append(spl[2])

		[t,Q,ri,diss,B] = np.loadtxt(dir+file, unpack=True)

		rimax = '%5.2f' % ri[-1]

		imin = np.nonzero(B)[0][0]
		imax = imin + np.nonzero(B[imin:]==0)[0][0]
		B_avg = '%5.2f' % (sum(B[imin:imax])/(imax-imin))
		tmin =	'%5.2f' % (4.5-t[imin])
		tmax =	'%5.2f' % (4.5-t[imax])
		status = False
		if B[-1] > 0:
			status = True

		str.append(rimax)
		str.append(tmin)
		str.append(tmax)
		str.append(B_avg)
		str.append(status)

		print str
