def makeRunPlots(run_id):
	import matplotlib.pyplot as p
	import numpy as np
	from matplotlib import rcParams
	import os

	rcParams['font.sans-serif']='Arial'
	rcParams['font.size']='7'
	fig_folder = '../fig/'
	dat_folder = '../out/'
	filename = dat_folder+run_id+'.dat'

	prefix = run_id+'_'

	[t,Qcmb,inner,diss,B,comp] = np.loadtxt(filename,unpack=True)

	p.figure(figsize=(3.34,2.56), dpi=300)
	p.plot(t,B,'k',linewidth=2)
	p.xlabel('Time [Ga]')
	p.ylabel('Surface magnetic field [$\mu$T]')
	p.savefig(fig_folder+prefix+'magnetic.eps', format='eps', bbox_inches='tight')
	
	p.figure(figsize=(3.34,2.56), dpi=300)
	p.xlabel('Time [Ga]')
	p.ylabel('Relative inner core size [$r_c$]')
	p.plot(t,inner,'k',linewidth=2)
	p.savefig(fig_folder+prefix+'innerCore.eps', format='eps', bbox_inches='tight')
	
	p.figure(figsize=(3.34,2.56), dpi=300)
	p.plot(t,Qcmb,'k',linewidth=2)
	p.xlabel('Time [Ga]')
	p.ylabel('$Q_{CMB}$ [W]')
	p.savefig(fig_folder+prefix+'Qcmb.eps', format='eps', bbox_inches='tight')

	p.figure(figsize=(3.34,2.56), dpi=300)
	p.plot(t,diss,'k',linewidth=2)
	p.plot(t,diss*0,'k--',linewidth=2)
	p.xlabel('Time [Ga]')
	p.ylabel('Dissipation []')
	p.savefig(fig_folder+prefix+'diss.eps', format='eps', bbox_inches='tight')

def makeCompPlots(*data):
	import matplotlib.pyplot as p
	from matplotlib import rcParams
	import numpy as np

	fig_folder = '../fig/'
	files = data[0]
	leg = data[1]

	p.figure(figsize=(3.34,2.56), dpi=300)

	for i in range(len(files)):
		[t,Qcmb,inner,diss,B] = np.loadtxt(files[i],unpack=True)
		p.plot(t,inner,linewidth=2)
		p.xlabel('Time [Ga]')
		p.ylabel('Relative inner core size [$r_c$]')

	p.legend(leg,loc='best')
	p.savefig(fig_folder+'comparison.eps', format='eps', bbox_inches='tight')

	return

import sys

if len(sys.argv) != 1:
	run_id = sys.argv[1]
	makeRunPlots(run_id)
