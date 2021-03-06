def readConfig(filename):
	import ConfigParser as p
	config = p.RawConfigParser()
	config.read(filename)

	global run, Rc, Rp, G, mu, k, alpha, rho, cp, drho, LH
	global Tm0, Tm1, Tm2, alphac, P0, X0, Xmax, delta, fudge
	global data_file, rho0, rhos0, alphad

	run = config.get('general','run')
	data_file = run+'.dat'

	Rc = float(config.get('general','Rc'))
	Rp = float(config.get('general','Rp'))
	G = float(config.get('general','G'))
	mu = float(config.get('general','mu'))
	k = float(config.get('general','k'))
	alpha = float(config.get('general','alpha'))
	rho = float(config.get('general','rho'))
	cp = float(config.get('general','cp'))
	drho = float(config.get('general','drho'))
	LH = float(config.get('general','LH'))
	Tm0 = float(config.get('general','Tm0'))
	Tm1 = float(config.get('general','Tm1'))
	Tm2 = float(config.get('general','Tm2'))
	alphac = float(config.get('general','alphac'))
	alphad = float(config.get('general','alphad'))
	P0 = float(config.get('general','P0'))
	X0 = float(config.get('general','X0'))
	Xmax = float(config.get('general','Xmax'))
	delta = float(config.get('general','delta'))
	fudge = float(config.get('general','fudge'))

	rho0 = rho
	rhos0 = 7.5e3


def timeEvolution():
	import numpy as np
	import math as m
	from scipy.optimize import fsolve

	global run, Rc, Rp, G, mu, k, alpha, rho, cp, drho, LH
	global Tm0, Tm1, Tm2, alphac, P0, X0, Xmax, delta, fudge
	global D, Mc, T0, dt, Qc, t, B, inner, Qcmb
	global data_file
	
	
	indat_folder = '../dat/'
	outdat_folder = '../out/'
	debug = 1

	# second order variables (depend on others)
	D = m.sqrt(3*cp/(2*m.pi*alpha*rho*G))	# scale height (depend on the other parameters actually)
	Mc = 4*m.pi*rho*Rc**3/3								# mass of the whole core
	
	filename = indat_folder+run+'.dat'
	f = open(filename,'r')
	f.readline()
	f.readline()
	dR = float(f.readline())
	f.close
	
	### STEP 1 ###
	#
	# load Qcmb(t)
	#
	
	data = np.loadtxt(filename,skiprows=3)
	t = np.linspace(0,0,len(data))
	Qcmb = np.linspace(0,0,len(data))
	k_mantle = 4
	
	for i,row in enumerate(data):
		t[i] = row[0]
		Tcm = row[1]
		Tm = row[2]
		Q = k_mantle*(Tcm-Tm)/(dR/2)*4*m.pi*Rc**2/2
		Qcmb[i] = Q
		if i==0:
			T0 = Tcm
	
	Qcmb[0] = Qcmb[1]
	#dt = (t[-1]-t[-2])*1e9*365*24*3600.0
	
	
	
	### STEP 2 ###
	#
	# function definitions
	#
	
	
	# pressure as a function of r
	def pres(r):
		global P0,G,rho
		return P0-2*m.pi*r**2*rho**2*G/3
	
	# outer core light element content as a function of ri
	def compo(ri):
		global X0,Xmax,Rc
		return min(X0/(1-(ri/Rc)**3),Xmax)
	
	# solidus temperature as a function of r and ri (due to light elements)
	def Tsol(r,ri):
		global Tm0,alphac,Tm1,Tm2
		return (1-alphac*compo(ri))*(Tm0+Tm1*pres(r))
		#return Tm0*(1-alphac*compo(ri))*(1+Tm1*pres(r)+Tm2*pres(r)**2)
	
	# adiabatic temperature profile
	def Tad(T0,r):
		global D,Rc
		return T0*m.exp(-(r**2-Rc**2)/D**2)
	
	# gravity factor from Nimmo's calculation
	def factor(f):
		a = (0.2+2*f**5/15-f**2/3)
		b = f/(1-f**3)
		return a*b
	
	# secular coefficient (see Nimmo09's table 1)
	def Q_secular():
		global Mc,cp,Rc,D
		return Mc*cp*(1+2*Rc**2/(5*D**2))
	
	# latent heat coefficient (see Nimmo09's table 1)
	def Q_latent(f,Tc):
		global Mc,LH,D,Rc,delta
		return 1.5*Mc*f*LH*D**2/(Tc*Rc**2*(delta-1))
	
	# compositional coefficient (see Nimmo09's table 1)
	def Q_compo(f,Tc):
		global G,rho,Mc,drho,D,delta
		return 3*m.pi*G*rho*Mc*factor(f)*drho*D**2/(Tc*(delta-1))
	
	# dissipation term (see Nimmo09's table 1)
	def E_phi(f,dTc):
		global G,rho,Mc,drho,D,delta,T0,dt,k,LH,cp,Rc
		g = 3*m.pi*G*rho*Mc*factor(f)*drho*D**2/(delta-1)/T0**2*dTc/dt
		l = 1.5*Mc*f*(1-f**2)*LH*dTc/dt/((delta-1)*T0**2)
		s = Mc*cp*2*Rc**2/(5*D**2*T0)*dTc/dt
		adiabat = Mc*12*k*Rc**2/(5*rho*D**4)
		return s+g+l-adiabat
	
	# function to solve to obtain ri at a given timestep
	def calcInnerCore(ri):
		global T0,Rc,dt,Qc
		T = T0-dt*Qc/(Q_secular()+Q_latent(ri/Rc,T0)+Q_compo(ri/Rc,T0))
		return Tsol(ri,ri)-Tad(T,ri)
	
	def update_parameters(ri,Tc):
		global rho0, rhos0, alphac, cp, alpha, G, Tm0 # regular params
		global rho, drho, D, delta # to be updated
		rho = rho0 - 4e3*compo(ri)
		rhos = rhos0 - alphad*1e3*compo(ri)
		drho = (rhos-rho)/rhos
		D = m.sqrt(3*cp/(2*m.pi*alpha*rho*G))
		delta = Tm1*(1-alphac*compo(ri))*rho*cp/(alpha*Tc)
		print drho

	### STEP 3 ###
	#
	# evolution calculations
	#
	
	# start of the evolution loop
	Tcore = np.linspace(T0,T0,len(t))
	inner = np.linspace(0,0,len(t))
	c = np.linspace(X0,X0,len(t))
	diss = np.linspace(0,0,len(t))
	B = np.linspace(0,0,len(t))
	Qsec = np.linspace(0,0,len(t))
	Qg = np.linspace(0,0,len(t))
	QL = np.linspace(0,0,len(t))

	i = 0
	while i < len(t)-1:

		update_parameters(inner[i],Tcore[i])

		Tcore[i] = T0
		if debug == 1:
			print i,t[i],T0,Tcore[i],Qcmb[i],inner[i]/Rc,c[i]
		Qc = Qcmb[i]

		if Qc < 0.0:
			ri = inner[i]
			inner[i+1] = ri
			diss[i+1] = E_phi(ri/Rc,0)*Tad(T0,ri)/(4*m.pi*(Rc**3-ri**3)/3)
			B[i+1] = 0.0
			c[i+1] = c[i]
			i += 1
			continue
	
		dt = (t[i+1]-t[i])*1e9*365*24*3600.0
		dT = dt*Qc/(Q_secular()+Q_latent(inner[i]/Rc,T0)+Q_compo(inner[i]/Rc,T0))
		dist_m = Tad(T0,0)-Tsol(0,0)

		if dist_m < dT and dist_m > 0:
			dQ = dist_m/dT
			dT = dist_m
			#Qcmb[i] = (1-dQ)*Qcmb[i] # Qcmb is already 'per second'
			t[i] = t[i] + dt*dQ/(1e9*365*24*3600.0)
			T0 = Tsol(0,0.0)/m.exp(Rc**2/D**2)-1e-4
			#print 'done secular only'
			continue

		import sys
		import numpy as np


		if dist_m <= 0:
			#print 'crystallizing'

			tmp = np.linspace(0,Rc,100)
			tmp2 = np.linspace(0,Rc,100)
			for z in range(len(tmp)):
				tmp2[z] = calcInnerCore(tmp[z])

			#print tmp
			#print tmp2
	
			#sys.exit()

			ri = fsolve(calcInnerCore,Rc/2)
			dT = dt*Qc/(Q_secular()+Q_latent(ri/Rc,T0)+Q_compo(ri/Rc,T0))
			if dT < 0:
				print 'halt'
				print ri,T0,dt,Qc
			diss[i+1] = E_phi(ri/Rc,dT)*Tad(T0-dT,ri)/(4*m.pi*(Rc**3-ri**3)/3)
			print E_phi(ri/Rc,dT)
			if diss[i+1] > 0:
				# scaling law from Aubert & Christensen 09
				B[i+1] = fudge*m.pow(rho*mu**3*diss[i+1]**2*(Rc-ri)**2,1./6.)*1e6*(Rc/Rp)**3
			else:
				B[i+1] = 0.0

			inner[i+1] = ri[0]
			c[i+1] = compo(ri[0])

		Qsec[i+1] = Q_secular()*dT/dt
		Qg[i+1] = Q_compo(inner[i]/Rc,T0)*dT/dt
		QL[i+1] = Q_latent(inner[i]/Rc,T0)*dT/dt

		T0 = T0 - dT
		i +=1

	fname = outdat_folder+data_file
	head = 't, qcmb, Tcore, inner/Rc, diss, B, c, Qsec, Qg, QL'
	data = (t,Qcmb,Tcore,inner/Rc,diss,B,c,Qsec,Qg,QL)
	np.savetxt(fname, np.column_stack(data), fmt='%5.3e', comments='# ',header= head, delimiter='\t')
