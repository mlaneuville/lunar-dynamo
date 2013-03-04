import thermo as t
import makeConfig as c
import makePlots as p
import numpy as np

runs = []

# [run,variable,value_min,value_max,#pts,value_scaling]
runs.append(['TWP0LB','k',			30,		40,		11,		1])
runs.append(['AWP0SB','k',			30,		40,		11,		1])
runs.append(['DWP0SB','Rc',			300,	400,	11,		1e3])
runs.append(['TWP0LB','rho',		6000,	8000,	21,		1])
runs.append(['TWP0LB','drho',		1,		10,		10,		1e-3])
runs.append(['TWP0LB','alpha',	5,		15,		11,		1e-5])
runs.append(['TWP0LB','LH',			300,	750,	10,		1e3])
runs.append(['TWP0LB','x0',			1,		10,		10,		1e-2])
runs.append(['TWP0LB','delta',	101,	130,	30,		1e-2])
runs.append(['TWP0LB','fudge',	5,		30,		26,		1e-2])

for r in runs:
	print r

	vrange = np.linspace(r[2],r[3],r[4])
	for i in range(r[4]):
		run = r[0]
		affix = r[1]+'_%04d' % vrange[i]
		cfile = '../out/'+run+'_'+affix+'.cfg'
		dfile = '../out/'+run+'_'+affix+'.dat'

		kw = {r[1]:vrange[i]*r[5]}
		c.makeConfig(run,affix,**kw)
		t.readConfig(cfile)
		t.timeEvolution()

		p.makeRunPlots(run+'_'+affix)
