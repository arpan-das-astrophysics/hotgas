import glob
import numpy as np
import math
#from scipy import interpolate 
import scipy
import pylab
from scipy import interpolate
import h5py
from scipy import stats

PI=3.14
hplank=6.62606896e-27
NU_over_EV=1.60217646e-12 / hplank
NUIONIZATION=13.60*NU_over_EV
HeI_NUIONIZATION=24.59*NU_over_EV
HeII_NUIONIZATION=(NUIONIZATION*4.0) 

kev = np.arange(0.1,10,.01)
pi= 3.14159
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

with h5py.File('redshift14.h5', 'w') as hf:
	index=np.array([0,4])
	rvir=np.array([0.324008,0.324008])
	for v in range(len(index)):
		list_of_files = glob.glob('/astro/home/arpan.das/enzo2/redshift14/halo{0}/density/*.dat'.format(index[v]))
		print len(list_of_files)
		opacity_final=[]
		cd=[]
		for x in range(1,len(list_of_files)):

			data1=np.loadtxt('/astro/home/arpan.das/enzo2/redshift14/halo{0}/density/los{1}.dat'.format(index[v],x), skiprows=5)
			data2=np.loadtxt('/astro/home/arpan.das/enzo2/redshift14/halo{0}/hiifraction/los{1}.dat'.format(index[v],x), skiprows=5)
			data3=np.loadtxt('/astro/home/arpan.das/enzo2/redshift14/halo{0}/metallicity/los{1}.dat'.format(index[v],x), skiprows=5)
			density=np.array(data1[:,5])
			hiifraction=np.array(data2[:,5])
			metallicity=np.array(data3[:,5])
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
			for m in range(len(metallicity_new)):
				if metallicity_new[m]<1.0e-4:
					metallicity_new[m]=0.0
				else:
					metallicity_new[m]=metallicity_new[m]


			column_density=np.trapz(density_new*(1-hiifraction_new),length_new)*3.086e+21*0.92
			cd.append(column_density)
			column_hydrogen=np.trapz(density_new,length_new)*3.086e+21*0.92
			column_metal=np.trapz(density_new*metallicity_new,length_new)*3.086e+21*0.92
			tau=column_density*cross_section_HI+column_hydrogen*cross_section_HeI/10+column_metal*cross_section_metals
			opacity=np.exp(-tau)
			opacity_final.append(opacity)

		cd=np.array(cd)
		opacity_final=np.array(opacity_final)
		cd_avg=np.sum(cd)/len(list_of_files)
#		opacity_avg=np.sum(opacity_final,axis=0)/len(list_of_files)
		g1 = hf.create_group('halo{0}'.format(index[v]))
		g1.create_dataset('opacity_los', data = opacity_final)
		g1.create_dataset('NHI', data = cd_avg)
