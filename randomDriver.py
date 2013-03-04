import thermo as t
import makeConfig as c
import numpy as np
import math as m
import random as r
import sys, glob

if len(sys.argv) != 2:
	print "one argument needed"
	sys.exit()



N = int(sys.argv[1])

param_space = dict()

param_space['runs'] = ['TWP0LB','TWP0LD','AWP0SB']
param_space['previous'] = dict()

param_space['Rc'] = np.linspace(300e3,400e3,11)
param_space['k'] = np.linspace(30,40,11)
param_space['alpha'] = np.linspace(5e-5,15e-5,11)
param_space['rho'] = np.linspace(6000,8000,21)
param_space['drho'] = np.linspace(0.01,0.10,10)
param_space['LH'] = np.linspace(3e5,7.5e5,10)
param_space['x0'] = np.linspace(0.01,0.10,10)
param_space['delta'] = np.linspace(1.01,1.30,30)
param_space['fudge'] = np.linspace(5e-2,30e-2,11)

for j in param_space['runs']:
	param_space['previous'][j] = len(glob.glob1('../out/',j+"_random_*.dat"))
	print j,param_space['previous'][j]


for j in range(N):

	kw = dict()
	for i,v in enumerate(param_space):
		if v == 'runs':
			run = r.sample(param_space[v],1)[0]
			affix = 'random_%05d' % int(param_space['previous'][run]+1)
			param_space['previous'][run] += 1
		elif v == 'previous':
			continue
		else:
			kw[v] = "%5.3e" % float(r.sample(param_space[v],1)[0])

	cfile = '../out/'+run+'_'+affix+'.cfg'
	dfile = '../out/'+run+'_'+affix+'.dat'

	n = int(m.ceil(m.log10(N)))
	print str(j+1).zfill(n)+"/"+str(N),cfile

	c.makeConfig(run,affix,**kw)
	t.readConfig(cfile)
	t.timeEvolution()
