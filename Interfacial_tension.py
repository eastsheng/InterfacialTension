'''
# calculating the interfacial tension by the pxx, pyy, pzz obtained from MD simulation
# interfacial tension: gamma = 0.5*Lz*(pzz-0.5*(pxx+pyy))
# unit for real 
'''
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def ifxyz(line,xyz):
	if xyz in line:
		# print(line)
		line = line.strip().split()
		lo = float(line[0])
		hi = float(line[1])
		return lo,hi
def read_size(lammpsdata):
	with open(lammpsdata,'r') as d:
		for index, line in enumerate(d,1):
			if "xlo" in line:
				xlo,xhi = ifxyz(line,"xlo")
			elif "ylo" in line:
				ylo,yhi = ifxyz(line,"ylo")
			elif "zlo" in line:
				zlo,zhi = ifxyz(line,"zlo")
		Lx = xhi-xlo
		Ly = yhi-ylo
		Lz = zhi-zlo
	return Lx,Ly,Lz

# def interfacial_tension(p_data,lammpsdata,zbin=1.0):
# 	atm2pa = 101325
# 	A2m = 1e-10
# 	m,n = p_data.shape
# 	# print(p_data.shape)
# 	Lx,Ly,Lz = read_size(lammpsdata)
# 	# Lz = (p_data[m-1,1]-p_data[0,1])
# 	# print(Lx,Ly,Lz)
# 	dv = Lx*Ly*zbin
# 	# print(dv)
# 	n_bin = p_data[:,2]
# 	pxx = -p_data[:,4]*n_bin/dv
# 	pyy = -p_data[:,5]*n_bin/dv
# 	pzz = -p_data[:,6]*n_bin/dv
# 	pT = (pxx + pyy)*0.5
# 	r = (pzz - pT)*0.5
# 	gamma = np.sum(r)*A2m*atm2pa #pa*m = pa.m*2/m --> N/m*1000=mN/m
# 	print("Interface tension =",round(gamma*1000,2),"mN/m")
# 	x = p_data[:,1]
# 	return pzz,pT,r,x,gamma

def interfacial_tension_2(p_data,lammpsdata,zbin=1.0):
	atm2Mpa = 0.101325
	A2m = 1e-10
	m,n = p_data.shape
	Lx,Ly,Lz = read_size(lammpsdata)
	# Lz = (p_data[m-1,1]-p_data[0,1])
	# print(Lx,Ly,Lz)
	dv = Lx*Ly*zbin
	# print(dv)
	n_bin = p_data[:,2]
	pxx = -p_data[:,4]*atm2Mpa*n_bin/dv
	pyy = -p_data[:,5]*atm2Mpa*n_bin/dv
	pzz = -p_data[:,6]*atm2Mpa*n_bin/dv
	pT = (pxx + pyy)*0.5
	r = (pzz - pT)*0.5
	gamma_z = np.cumsum(r)*A2m*1e9 # 1 Mpa*m = 1e6 N/m = 1e9 mN/m
	gamma = np.sum(r)*A2m*1e9
	print("Interface tension =",round(gamma,2),"mN/m")
	x = p_data[:,1]
	# print(gamma_z)
	return pzz,pT,r,x,gamma_z,gamma

def plot_pressure():
	color = ["k:","r:","g:","b:","c:","o:"]
	plt.rc('font', family='Times New Roman', size=18)
	fig = plt.figure(figsize=(12,8))
	fig.subplots_adjust(bottom=0.15,left=0.15)
	# ax=fig.add_subplot(111)
	# Path("./oil_water/").mkdir(parents=True,exist_ok=True)
	lammpsdata = "./oil_water/1_nvt_300.data"
	
	p_data = np.loadtxt("./oil_water/1_press_all.profile",skiprows=4)
	pN,pT,r,x,gamma_z,gamma = interfacial_tension_2(p_data,lammpsdata,1)

	ax=fig.add_subplot(111)
	ax.plot(x,pN,color[1],label="pN",linewidth=1)
	ax.plot(x,pT,color[2],label="pT",linewidth=1)
	ax.plot(x,pN-pT,color[3],label="pN-pT",linewidth=1)
	ax.legend(loc="best")
	ax.set_xlabel("z (Å)",fontweight="bold",fontsize=22)
	ax.set_ylabel("Pressure (MPa)")

	ay = ax.twinx()
	ay.plot(x,gamma_z,"r-",label="gamma_z",linewidth=2)
	# ay.set_xlabel("z (Å)")
	ay.set_ylabel("IFT (mN/m)",color="r",fontweight="bold",fontsize=22)
	ay.tick_params(axis='y',colors='red')
	ay.spines['right'].set_color('red')
	ay.legend(loc="upper center")
	plt.savefig("./oil_water/interfacial_tension.png",dpi=300)
	plt.show()
	
	return

if __name__ == '__main__':
	plot_pressure()

	