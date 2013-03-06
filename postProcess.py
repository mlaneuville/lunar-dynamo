import numpy as np
import makePlots as p
import ConfigParser as c
import os, glob

dir = '../out/'
stats = dict()
stats['TWP0LW'] = dict()
stats['TWP0LD'] = dict()
stats['TWP0LB'] = dict()
stats['AWP0SB'] = dict()
params = ['Rc','k','alpha','rho','drho','LH','x0','delta','fudge']
for item in params:
	stats['TWP0LW'][item] = []
	stats['TWP0LD'][item] = []
	stats['TWP0LB'][item] = []
	stats['AWP0SB'][item] = []

for file in os.listdir(dir):
	str = []
	if file.find('dat') != -1:
		str.append(file)
		spl = file.split('_')
		spl[-1] = str[-1][:-4]
		str.append(spl[0])
		str.append(spl[1])
		str.append(spl[2])

		[t,Q,ri,diss,B,compo,Qs,Qg,Ql] = np.loadtxt(dir+file, unpack=True)

		rimax = '%5.2f' % ri[-1]

		imin = np.nonzero(B)[0][0]
		dynamo = np.nonzero(B[imin:]==0)
		if len(dynamo[0]) != 0:
			imax = imin + np.nonzero(B[imin:]==0)[0][0]
		else:
			imax = len(B)-1
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

		if float(B_avg) > 1.0 and status == False:
			print str
			config = c.RawConfigParser()
			config.read('../out/'+str[3]+'.cfg')
			for item in params:
				stats[str[1]][item].append(float(config.get('normalized',item)))
				
			if len(glob.glob1('../fig/',str[3]+"_*.eps")) == 5:
				continue
			p.makeRunPlots(str[3])

for i in stats.keys():
	p.makeStatPlots(i,stats[i])
