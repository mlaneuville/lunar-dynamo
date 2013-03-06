import thermo as t
import makeConfig as c
import ConfigParser as p
import numpy as np
import math as m
import random as r
import sys, glob

if len(sys.argv) != 2:
	print "one argument needed"
	sys.exit()



N = int(sys.argv[1])
config = p.RawConfigParser()
config.read('main.cfg')

param_space = dict()

param_space['runs'] = ['TWP0LB','TWP0LD','AWP0SB','TWP0LW']
param_space['previous'] = dict()

param_space['Rc'] = np.linspace(config.getfloat('var_min','Rc'),config.getfloat('var_max','Rc'),100)
param_space['k'] = np.linspace(config.getfloat('var_min','k'),config.getfloat('var_max','k'),100)
param_space['alpha'] = np.linspace(config.getfloat('var_min','alpha'),config.getfloat('var_max','alpha'),100)
param_space['rho'] = np.linspace(config.getfloat('var_min','rho'),config.getfloat('var_max','rho'),100)
param_space['drho'] = np.linspace(config.getfloat('var_min','drho'),config.getfloat('var_max','drho'),100)
param_space['LH'] = np.linspace(config.getfloat('var_min','LH'),config.getfloat('var_max','LH'),100)
param_space['x0'] = np.linspace(config.getfloat('var_min','x0'),config.getfloat('var_max','x0'),100)
param_space['delta'] = np.linspace(config.getfloat('var_min','delta'),config.getfloat('var_max','delta'),100)
param_space['fudge'] = np.linspace(config.getfloat('var_min','fudge'),config.getfloat('var_max','fudge'),100)

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
