import thermo as t
import makeConfig as c
#import makePlots as p
import numpy as np
import random as r
import sys

if len(sys.argv) != 2:
	print "one argument needed"
	sys.exit()

N = int(sys.argv[1])

param_space = dict()

param_space['runs'] = ['TWP0LB']
param_space['Rc'] = np.linspace(300e3,400e3,11)
param_space['k'] = np.linspace(30,40,11)
param_space['alpha'] = np.linspace(5e-5,15e-5,11)
param_space['rho'] = np.linspace(6000,8000,21)
param_space['drho'] = np.linspace(0.01,0.10,10)
param_space['LH'] = np.linspace(3e5,7.5e5,10)
param_space['x0'] = np.linspace(0.01,0.10,10)
param_space['delta'] = np.linspace(1.01,1.30,30)
param_space['fudge'] = np.linspace(5e-2,30e-2,11)



for j in range(N):

	affix = 'random_%03d' % j
	kw = dict()
	for i,v in enumerate(param_space):
		if v == 'runs':
			run = r.sample(param_space[v],1)[0]
		else:
			kw[v] = "%5.3e" % float(r.sample(param_space[v],1)[0])

	cfile = '../out/'+run+'_'+affix+'.cfg'
	dfile = '../out/'+run+'_'+affix+'.dat'

	c.makeConfig(run,affix,**kw)
	t.readConfig(cfile)
	t.timeEvolution()
