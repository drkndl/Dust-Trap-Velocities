import numpy as np 
import matplotlib.pyplot as plt
import os
from matplotlib.font_manager import FontProperties
from find_ring import load_dust_outputs, load_gas_outputs, get_dust_trap



G = 6.67e-11		# SI Gravitational Constant
M = 1.989e30		# mass of the Sun in kg (the default MSTAR in FARGO3D)
R = 5.2*1.4959e11  	# 5.2 AU
GAMMA = 1.6667		# Adiabatic index
CS = 1000			# Speed of sound taken as 1 km/s



def omega_kepler(rmin, rmax, Ny):
	"""
	Function to calculate Kepler velocities
	"""

	omega_k = []
	v_k = []
	rarr = np.linspace(rmin, rmax, Ny)

	for r in rarr:
		omega_k.append(np.sqrt(G*M/(r*R)**3))
		v_k.append(np.sqrt(G*M/(r*R)))

	return omega_k, v_k




def v_rad_phi(path, r, phi, n_species, n_out):
	"""
	Getting radial and azimuthal velocities separately
	"""

	vphic = []
	vradc = []

	for i in range(n_species):

		# Loading radial and azimuthal velocities for each dust species
		vphi = np.fromfile(os.path.join(path, "dust"+str(i+1)+"vx"+str(n_out)+".dat")).reshape(len(r)-1, len(phi)-1)
		vrad = np.fromfile(os.path.join(path, "dust"+str(i+1)+"vy"+str(n_out)+".dat")).reshape(len(r)-1, len(phi)-1)
		
		# Centering since the grid is staggered
		vphic.append(0.5*(vphi[:-1,1:]+vphi[:-1,:-1]))
		vradc.append(0.5*(vrad[1:,:-1]+vrad[:-1,:-1]))

	return vphic, vradc 




def plots(vphic, vradc):

	stokes = np.logspace(-4, 2, 127)

	for i in range(len(vphic)):

		plt.plot(stokes, vradc[i].mean(axis=1))
		plt.xscale('log')

	plt.show()




def plot_vel(vk, vdust, vgas, Ny, path):
	"""
	Plotting the actual output velocities and the calculated Keplerian velocities of the gas and the dust species
	"""

	plt.plot(range(Ny), vk, alpha=0.5, label="Keplerian velocity")

	for i in range(len(vdust)):
		plt.plot(range(Ny-1), vdust[i].mean(axis=1), alpha=0.5, label="Species %s"% str(i+1))

	plt.plot(range(Ny-1), vgas.mean(axis=1), alpha=0.5, label="Gas velocity")

	plt.xlabel("Radius (Ny)")
	plt.ylabel("Velocity (m/s)")
	plt.title("Azimuthally Averaged Velocities: Keplerian, Gas & Dust")
	plt.legend(prop={'size':7})
	plt.savefig(path + "kep_vs_actual.png")
	plt.show()




def stopping_time(stokes, omega_k, dt):

	t_stop = []
	ok_trap = omega_k[dt]

	for i in stokes:
		t_stop.append(i/ok_trap)

	return t_stop



def pressure(sigma, Ny, energy, rmin, rmax):

	# Calculating gas pressure
	rarr = np.linspace(rmin, rmax, Ny)
	rho = []
	pressure = []

	for i in range(len(sigma)):
		rho.append(sigma[i].mean(axis=1)/(np.sqrt(2)*np.pi*0.05*rarr))		# 0.05 is the Aspect Ratio I used in FARGO3D
		pressure.append(rho[i]*CS**2)

	gas_pressure = (GAMMA-1)*np.asarray(energy)

	plt.plot(range(Ny), gas_pressure.mean(axis=1), linewidth=10, alpha=0.5, label="Gas")

	# for i in range(len(pressure)):
	# 	plt.plot(range(Ny), pressure[i], label=f"Dust {i+1}")

	plt.title("Azimuthally Averaged Gas Pressure")
	plt.xlabel("Radius [Ny]")
	plt.ylabel("Pressure (Pa)")
	plt.legend()
	plt.show()



if __name__ == "__main__":

	rmin = 0.4
	rmax = 2.5
	Ny = 128
	path = "./Set 7/fargo_multifluid/"
	pic_path = "./Set 7/pics/"
	output_number = 200
	species_number = 15
	p_bumps = 2

	r, phi, sigma, vel, energy = load_dust_outputs(path, species_number, output_number)
	gas_sig, gas_vel, gas_energy = load_gas_outputs(r, phi, path, output_number)

	dt = get_dust_trap(sigma, pic_path, p_bumps=p_bumps)

	stoke = np.logspace(-5, 3, 15)

	ok, vk = omega_kepler(rmin, rmax, Ny)
	plot_vel(vk, vel, gas_vel, Ny, pic_path)

	pressure(sigma, Ny, gas_energy, rmin, rmax)

	vphic, vradc = v_rad_phi(path, r, phi, species_number, output_number)
	# plots(vphic, vradc)

	for i in dt:

		stoptime = stopping_time(stoke, ok, i)
		print(i, stoptime)


