#export HEADAS=/home/arpan/heasoft-6.19/x86_64-unknown-linux-gnu-libc2.23-0
#. $HEADAS/headas-init.sh

import glob
import numpy as np
import math                                                                                                        
import scipy
from scipy import interpolate
from xspec import *
import matplotlib.pyplot as plt
import matplotlib
from astropy import *


AllData.dummyrsp(0.3,10,990,"lin")
m1=Model("mekal")
par1=m1(1)
par1.values=[1,.01,0.008617328149741,0.008617328149741,79.9,79.9]
Plot.xAxis = "keV"
Plot("model")
x=Plot.x()
kev=np.array(x) 

PI=3.14
hplank=6.62606896e-27
NU_over_EV=1.60217646e-12 / hplank
NUIONIZATION=13.60*NU_over_EV
HeI_NUIONIZATION=24.59*NU_over_EV
HeII_NUIONIZATION=(NUIONIZATION*4.0) 
h=4.135667662e-18
c=299792458
hertz=2.417990504024e+17*kev

#####################################################
#Define the Hydrogen, Helium and metal cross sections
######################################################
def HI(frequency,result=[]):
	for nu in frequency:
		if (nu < NUIONIZATION):
			result.append(0)	    
		else:
			epsilon = math.sqrt( nu/NUIONIZATION - 1)
			#epsilon =math.sqrt( e/0.0136 - 1)
			result.append((6.3e-18) * math.pow(NUIONIZATION/nu, 4) * math.exp(4-(4*math.atan(epsilon)/epsilon)) / (1-math.exp(-2*PI/epsilon)))
	return result

def HeI(frequency,result=[]):
	for nu in frequency:
		if (nu < HeI_NUIONIZATION):
			result.append(0)
		else:
			x = nu/NU_over_EV/13.61 - 0.4434
			y = math.sqrt(x*x + math.pow(2.136, 2))
			result.append(9.492e-16*((x-1)*(x-1) + 2.039*2.039) * math.pow(y, (0.5 * 3.188 - 5.5))* math.pow(1.0 + math.sqrt(y/1.469), -3.188))

	return result

def metals(energy,result=[]):
	for e in energy:
		if(0.030<= e < 0.100):
			result.append((17.3 + 608.1 * e - 2150.0 * e * e)*math.pow(e,-3)*1.0e-24)
		elif(0.100<= e <.284):
			result.append((34.6 + 267.9 * e- 476.1 * e * e)*math.pow(e,-3)*1.0e-24)
		elif(0.284<= e <0.400):
			result.append((78.1+18.8*e+4.3*e*e)*math.pow(e,-3)*1.0e-24)
		elif(0.400<= e <0.532):
			result.append((71.4+66.8*e-51.4*e*e)*math.pow(e,-3)*1.0e-24)
		elif(0.533<= e <0.707):
			result.append((95.5+145.8*e-61.1*e*e)*math.pow(e,-3)*1.0e-24)
		elif(0.707<= e <0.867):
			result.append((308.9-380.6*e+294.0*e*e)*math.pow(e,-3)*1.0e-24)
		elif(0.867<= e <1.303):
			result.append((120.6+169.3*e-47.7*e*e)*math.pow(e,-3)*1.0e-24)
		elif(1.303<= e <1.840):
			result.append((141.3+146.8*e-31.5*e*e)*math.pow(e,-3)*1.0e-24)
		elif(1.840<= e <2.471):
			result.append((202.7+104.7*e-17.0*e*e)*math.pow(e,-3)*1.0e-24)
		elif(2.471<= e <3.210):
			result.append((342.7+18.7*e)*math.pow(e,-3)*1.0e-24)
		elif(3.210<= e <4.038):
			result.append((352.2+18.7*e)*math.pow(e,-3)*1.0e-24)
		elif(4.038<= e <7.111):
			result.append((433.9 - 2.4*e + 0.75*e*e)*math.pow(e,-3)*1.0e-24)
		elif(7.111<= e <8.331):
			result.append((629.0 + 30.9*e)*math.pow(e,-3)*1.0e-24)
		else:
			result.append((701.2 + 25.2*e)*math.pow(e,-3)*1.0e-24)

	return result

cross_section_HI=np.array(HI(hertz))
cross_section_HeI=np.array(HeI(hertz))
cross_section_metals=np.array(metals(kev))-cross_section_HI-cross_section_HeI/10


indx=np.array([12])
rvir=np.array([1.19089])
sfr=1775/659.7e+6
for v in range(len(indx)):
	emissivity=np.array([0]*len(kev))
	list_of_files = glob.glob('/home/arpan/Desktop/hotgas/hotgas/redshift8/halo{0}/density/*.dat'.format(indx[v]))
	for x in range(1,len(list_of_files)):
		print x
		data1=np.loadtxt('/home/arpan/Desktop/hotgas/hotgas/redshift8/halo{0}/density/los{1}.dat'.format(indx[v],x), skiprows=5)
		data2=np.loadtxt('/home/arpan/Desktop/hotgas/hotgas/redshift8/halo{0}/hiifraction/los{1}.dat'.format(indx[v],x), skiprows=5)
		data3=np.loadtxt('/home/arpan/Desktop/hotgas/hotgas/redshift8/halo{0}/metallicity/los{1}.dat'.format(indx[v],x), skiprows=5)
		data4=np.loadtxt('/home/arpan/Desktop/hotgas/hotgas/redshift8/halo{0}/temperature/los{1}.dat'.format(indx[v],x), skiprows=5)
		density=np.array(data1[:,5])
                hiifraction=np.array(data2[:,5])
                metallicity=np.array(data3[:,5])
		temperature=np.array(data4[:,5])
                length=np.array(data1[:,7])*1000/(1+14.78)
                cut=[]
                for k in range(len(length)):
                        if length[k]<=rvir[v]:
                                cut.append(k)
                cut=np.array(cut)
                length_new=length[:np.amax(cut)]
                density_new=density[:np.amax(cut)]
                hiifraction_new=hiifraction[:np.amax(cut)]
                metallicity_new=metallicity[:np.amax(cut)]
		temperature_new=temperature[:np.amax(cut)]
                for m in range(len(metallicity_new)):
			if metallicity_new[m]<1.0e-5:
                                metallicity_new[m]=0.0
                        else:
                                metallicity_new[m]=metallicity_new[m]

		stepsize=[]
		for p in range(len(length_new)-1):
			stepsize.append(length_new[p+1]-length_new[p])
		stepsize=np.array(stepsize)
		opa=[]
		for q in range(len(stepsize)):
			column_density=density_new[q]*(1-hiifraction_new[q])*stepsize[q]*3.086e+21*0.92
			column_hydrogen=density_new[q]*stepsize[q]*3.086e+21*0.92
			column_metal=density_new[q]*metallicity_new[q]*length_new[q]*3.086e+21*0.92
			tau=column_density*cross_section_HI+column_hydrogen*cross_section_HeI/10+column_metal*cross_section_metals
			opacity=np.exp(-tau)
			opa.append(opacity)
		opa=np.array(opa)


		blank=[0]*len(kev)
		blankarray=np.array(blank)
		luminosity=[]	
		for index in range(len(opa)):
			temp_kev=temperature_new[index]*8.617328149741e-8
			if temp_kev>0.0087:
				m1.setPars(temp_kev,density_new[index]*0.92,metallicity_new[index],0,0,1)
				Plot("model")
				y=Plot.model()
				norm=density_new[index]*0.92*density_new[index]*hiifraction_new[index]*0.92*math.pow((stepsize[index]*3.086e+21),3.0)
				luminosity.append(norm*1.e-14*1.60218e-9*kev*np.array(y)/sfr)
			else:
				luminosity.append(blankarray)
		luminosity=np.array(luminosity)

		total=blankarray
		for ind in range(len(luminosity)):
			total=total*opa[ind]+luminosity[ind]
		emissivity=emissivity+total

	final=emissivity
	f = open('halo{0}.dat'.format(indx[v]),'w' )
	for g in range(len(kev)):
		f.write(str(kev[g]) + " " + str(final[g]) + "\n")	
	f.close()
