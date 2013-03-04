import ConfigParser as p
import shutil as sh
import math as m

def makeConfig(run,affix,Rc=390e3,k=30,alpha=9.2e-5,rho=7019,cp=835,drho=0.05,LH=750e3, \
								Tm0=1809,Tm1=1.54e-11,Tm2=-1.17e-22,alphac=2.4,x0=0.04,Xmax=0.4,delta=1.2,fudge=0.1):


	config = p.RawConfigParser()
	config.add_section('general')

	config.set('general','run',			run)			# run to study
	config.set('general','Rp',			1740e3)				# lunar radius
	config.set('general','G',				6.67e-11)			# gravity constant
	config.set('general','mu',			4*m.pi*1e-7)  # permeability of space
	config.set('general','P0',			8e9)					# assumed central pressure

	config.set('general','Rc',			Rc)						# core radius
	config.set('general','k',				k)						# core thermal conductivity
	config.set('general','alpha',		alpha)				#	core thermal expansivity
	config.set('general','rho',			rho)					# core density (assumed constant)
	config.set('general','cp',			cp)						# core specific heat
	config.set('general','drho',		drho)					# density change upon crystallisation (not clear)
	config.set('general','LH',			LH)						# latent heat of cristallisation
	config.set('general','Tm0',			Tm0)					# coef for melting curve
	config.set('general','Tm1',			Tm1)					# coef for melting curve
	config.set('general','Tm2',			Tm2)					# coef for melting curve
	config.set('general','alphac',	alphac)				# coef for melting curve
	config.set('general','X0',			x0)						# light element initial content
	config.set('general','Xmax',		Xmax)					# max light element in outer core (never reached)
	config.set('general','delta',		delta)				# dTdP/dTmdP (needs to be self consistent)
	config.set('general','fudge',		fudge)

	with open('../out/'+run+'_'+affix+'.cfg','w') as configfile:
		config.write(configfile)
