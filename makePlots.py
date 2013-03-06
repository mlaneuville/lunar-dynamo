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

	[t,Qcmb,inner,diss,B,comp,Qs,Qg,Ql] = np.loadtxt(filename,unpack=True)

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
	p.plot(t,Qcmb*0,'k--',linewidth=2)
	p.xlabel('Time [Ga]')
	p.ylabel('$Q_{CMB}$ [W]')
	p.savefig(fig_folder+prefix+'Qcmb.eps', format='eps', bbox_inches='tight')

	p.figure(figsize=(3.34,2.56), dpi=300)
	p.plot(t,diss,'k',linewidth=2)
	p.plot(t,diss*0,'k--',linewidth=2)
	p.xlabel('Time [Ga]')
	p.ylabel('Dissipation []')
	p.savefig(fig_folder+prefix+'diss.eps', format='eps', bbox_inches='tight')

	p.figure(figsize=(3.34,2.56), dpi=300)
	p.plot(t,Qs,linewidth=2)
	p.plot(t,Qg,linewidth=2)
	p.plot(t,Ql,linewidth=2)
	p.xlabel('Time [Ga]')
	p.ylabel('$Q_{i}$ [W]')
	p.legend(('Qs','Qg','Ql'),loc='best')
	p.yscale('log')
	p.savefig(fig_folder+prefix+'Qi.eps', format='eps', bbox_inches='tight')

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

def makeStatPlots(run,data):
	import matplotlib.pyplot as p
	import numpy as np
	print run

	if len(data['k']) == 0:
		return

	var = []
	variables = []
	for item in data.keys():
		var.append(data[item])
		variables.append(item)

	fig = p.figure(figsize=(10,6))
	ax1 = fig.add_subplot(111)
	p.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
	bp = p.boxplot(var, notch=0, sym='+', vert=1, whis=1.5)
	for i in range(len(data.keys())):
		med = bp['medians'][i]
		p.plot([np.average(med.get_xdata())], [np.average(var[i])], color='w', marker='*', markeredgecolor='k')
	p.setp(bp['boxes'], color='black')
	p.setp(bp['whiskers'], color='black')
	p.setp(bp['fliers'], color='red', marker='+')
	ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
	ax1.set_axisbelow(True)
	ax1.set_xlabel('Variable parameters',fontsize=14)
	ax1.set_ylabel('Relative value in range',fontsize=14)
	xtickNames = p.setp(ax1, xticklabels=variables)
	p.setp(xtickNames, fontsize=11)

	p.savefig(run)



import sys

if len(sys.argv) != 1:
	run_id = sys.argv[1]
	makeRunPlots(run_id)
