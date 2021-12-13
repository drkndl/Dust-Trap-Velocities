import numpy as np 
import matplotlib.pyplot as plt 
import os 
from scipy.signal import argrelextrema



def load_dust_outputs(path, n_species, n_out):
	"""
	This functions loads all the required fields such as surface density, velocities, energy of the dust species and the domains from the output files of the FARGO3D simulations. 

	Parameters:

	path:		path to the output folder of the simulations
	n_species:	number of dust species (integer)
	n_out:		output number (integer)

	Returns r, phi, sigma, velocity v, energy
	"""

	sigma = []
	energy = []
	v = []

	# load cell interface locations
	r = np.loadtxt(os.path.join(path, "domain_y.dat"))[3:-3] # ignore ghost cells
	phi = np.loadtxt(os.path.join(path, "domain_x.dat"))

	for i in range(n_species):

	    # Loading surface density for each dust species
		sigma.append(np.fromfile(os.path.join(path, "dust" + str(i+1) + "dens" + str(n_out) + ".dat")).reshape(len(r) - 1, len(phi) - 1))

		# Loading radial and azimuthal velocities for each dust species
		vphi = np.fromfile(os.path.join(path, "dust"+str(i+1)+"vx"+str(n_out)+".dat")).reshape(len(r)-1, len(phi)-1)
		vrad = np.fromfile(os.path.join(path, "dust"+str(i+1)+"vy"+str(n_out)+".dat")).reshape(len(r)-1, len(phi)-1)
		
		# Centering since the grid is staggered
		vphic = 0.5*(vphi[:-1,1:]+vphi[:-1,:-1])
		vradc = 0.5*(vrad[1:,:-1]+vrad[:-1,:-1])

		# Calculating total velocity v
		v.append(np.sqrt(vphic**2+vradc**2))

		# Loading energy for each dust species
		energy.append(np.fromfile(os.path.join(path, "dust" + str(i+1) + "energy" + str(n_out) + ".dat")).reshape(len(r) - 1, len(phi) - 1))

	# The shape of the 3 lists are:
	# sigma = (n_species, Ny, Nx)		e.g. (10, 128, 384)
	# v = (n_species, Ny-1, Nx-1)
	# energy = (n_species, Ny, Nx)

	return r, phi, sigma, v, energy



def load_gas_outputs(r, phi, path, n_out):
	"""
	This function loads all the required fields such as surface densities, energies and velocities of the gas from the output files of FARGO3D simulations. 

	Parameters:

	r:		y domain (a list of length Ny + 1)
	phi:	x domain (a list of length Nx + 1)
	path:	path to the output folder of the simulations
	n_out:	output number (integer)

	Returns gas surface density, gas velocity and gas energy.
	"""

	# load azimuthal and radial gas velocities
	vphi = np.fromfile(os.path.join(path, "gasvx" + str(n_out) + ".dat")).reshape(len(r)-1, len(phi)-1)
	vrad = np.fromfile(os.path.join(path, "gasvy" + str(n_out) + ".dat")).reshape(len(r)-1, len(phi)-1)

	# Centering since the grid is staggered
	vphic = 0.5*(vphi[:-1,1:]+vphi[:-1,:-1])
	vradc = 0.5*(vrad[1:,:-1]+vrad[:-1,:-1])

	# Calculating total velocity 
	vgas = (np.sqrt(vphic**2+vradc**2))

	# load gas surface density 
	sigma_gas = np.fromfile(os.path.join(path, "gasdens" + str(n_out) + ".dat")).reshape(len(r) - 1, len(phi) - 1)

	# load gas energy
	energy_gas = np.fromfile(os.path.join(path, "gasenergy" + str(n_out) + ".dat")).reshape(len(r) - 1, len(phi) - 1)

	# Shape of vgas = (Ny-1, Nx-1)
	# Shape of sigma_gas = (Ny, Nx)
	# Shape of energy_gas = (Ny, Nx)

	return sigma_gas, vgas, energy_gas




def get_dust_trap(sigm, pic_path, n_species, p_bumps=1):
	"""
	This function plots the average surface density radially for each dust species. The points of maxima are considered as the dust trapped regions. These points is returned. The number of dust traps are determined by the manually input number p_bumps.
	
	Parameters:

	sigm:		3-D list of the 2 dimensional surface densities of the various species
	pic_path:	output directory of the radial plot of average surface density
	n_species:	the total number of dust species called for no reason but to create the satisfying colour coordinated plots 
	p_bumps:	number of dust trapped regions to be considered (integer). Default is 1.

	Returns the radial indices of the dust trapped regions.
	"""

	dust_traps = []				# A list containing the maxima of surface density for each species
	dt = []

	# Creating a fancy colour map
	colormap = plt.cm.jet
	plt.gca().set_prop_cycle(plt.cycler('color', colormap(np.linspace(0, 1, n_species))))

	for i in range(len(sigm)):

		# Plots the mean surface density for each species
		plt.semilogy(sigm[i].mean(axis=1), label="S %s"% str(i+1))

		# Appending the radial index of the maximum point of surface density for each species  
		dust_traps.append(np.argmax(sigm[i].mean(axis=1)))

	if p_bumps==1:

		# Takes the dust trap outside the planet's orbit. It is usually the dust trap with the maximum density because the curve shoots up at St = 1
		dt.append(np.max(dust_traps))

	else:

		# Counts the unique locations of dust traps	
		uni, counts = np.unique(dust_traps, return_counts=True)

		# Sorts all the unique counts and returns p_bumps number of dust traps with the highest count
		indices = counts.argsort()[-p_bumps:][::-1]
		dt = uni[indices]

	# Indicating the dust trapped regions and their surrounding region by an offset of 5 radial indices
	for i in range(len(dt)):
		plt.axvline(dt[i], color='black', alpha=0.5)
		plt.axvline(dt[i]-5, linestyle='--', color='black', alpha=0.5)
		plt.axvline(dt[i]+5, linestyle='--', color='black', alpha=0.5)

	plt.title("Radial distribution of surface density of the disk")
	plt.xlabel("Radius (Ny)")
	plt.ylabel("Surface density ($kg/m^2$)")
	plt.legend(bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
	plt.savefig(pic_path+"surf_dens.png", bbox_inches='tight')
	plt.show()

	return dt



def cell_volume(r, phi, dt):
	"""
	This function calculates the cell volume in the dust trapped ring. 
	Cell volume (expectedly) changes across the rings but remains equal azimuthally. 

	Parameters:

	r:		y-domain of the output files (a list of 129 elements)
	phi:	x-domain of the output files (a list of 385 elements)
	dt:		radial index of the dust trapped ring (integer)

	Returns the cell volume in a ring (assuming the cell volume is azimuthally equal)
	"""

	r_m, phi_m = np.meshgrid(r, phi)		# Their shape is (385, 129)

	vol = np.pi*(r_m[:, dt+1]**2 - r_m[:, dt]**2)

	# vol is an array of areas around the ring. Each element in the array corresponds to the area in each of the 385 azimuthal cells Since the volume here is azimuthally symmetric, I consider the volume of every cell to be the same and therefore take the first cell.

	return vol[0]



def get_ring_param(sigma, v, sigma_g, v_g, dt, thick=True):
	"""
	This function finds the velocities and surface densities of the dust species and the gas only for the dust trapped ring. 

	Parameters:

	sigma:		3-D list of 2 dimensional surface densities for dust species
	v: 			3-D list of 2 dimensional velocities for dust species
	sigma_g:	2-D list of 2 dimensional gas surface density
	v_g:		2-D list of 2 dimensional gas velocities
	dt:			radial index of the dust trapped region (integer)
	thick:		toggles thick or thin ring (Boolean). A thick ring considers 5 radial indices on either side in the region 					surrounding the pressure maxima. A thin ring (thick=False) considers only the radial index of pressure maxima. 				By default, thick=True.

	Returns ring velocities and ring surface densities for the gas and the dust species
	"""

	ring_vs = []
	ring_sigma = []
	ring_gsigma = []
	ring_gvs = []
	

	if thick:

		# Getting ring parameters for gas
		ring_gsigma.append(sigma_g[dt - 5: dt + 5][:])	# Shape: (1, 10, Nx)
		ring_gvs.append(v_g[dt - 5: dt + 5][:])			# Shape: (1, 10, Nx-1)

		# Getting ring parameters for dust
		for i in range(len(sigma)):

			ring_vs.append(v[i][dt - 5: dt + 5][:])
			ring_sigma.append(sigma[i][dt - 5: dt + 5][:])
			# The shape of the ring velocity list is (n_species, 10, Nx-1)
			# The shape of the ring sigma list is (n_species, 10, Nx)

	else:

		# Getting ring parameters for gas
		ring_gsigma.append(sigma_g[dt][:])	# Shape: (Nx)
		ring_gvs.append(v_g[dt][:])			# Shape: (Nx-1)

		# Getting ring parameters for dust
		for i in range(len(sigma)):

			ring_vs.append(v[i][dt])
			ring_sigma.append(sigma[i][dt])
			# The shape of the ring velocity list is (n_species, Nx-1)
			# The shape of the ring sigma list is (n_species, Nx)

	return ring_vs, ring_sigma, ring_gvs, ring_gsigma




def vdiff_wavg(ring_vs, ring_sigma, cell_vol, thick=True):
	"""
	This function calculates the velocity differences between each set of species as a weighted average of both the species' surface densities.

	Parameters:

	ring_vs:	velocities of the dust species in the dust trapped region (3-D list for thick ring, 2-D list for thin ring)
	ring_sigma:	surface densities of the dust species in the dust ring (3-D list for thick ring, 2-D list for thin ring)
	cell_vol:	cell volume at the dust trapped region (float number)
	thick:		toggles thick or thin ring (Boolean). A thick ring considers 5 radial indices on either side in the region 					surrounding the pressure maxima. A thin ring (thick=False) considers only the radial index of pressure maxima. 				By default, thick=True. NEEDS TO BE THE SAME AS VALUE FOR get_ring_param.

	Returns v_wavg, a 1-D list of velocity differences for the various sets of dust species e.g (1,2), (1,3), (8, 10) etc
	"""		

	v_wavg = []

	for i in range(len(ring_vs)):

		for j in range(i+1, len(ring_vs)):

			if thick:

				int_upper = (abs(ring_vs[i] - ring_vs[j]) * ring_sigma[i][:, :-1] * ring_sigma[j][:, :-1] * cell_vol).sum()
				int_lower = (ring_sigma[i][:, :-1] * ring_sigma[j][:, :-1] * cell_vol).sum()

				wavg = int_upper/int_lower
			
				print(f"V{i+1}_{j+1} = {round(wavg, 3)} m/s")

				v_wavg.append(wavg)

			else:

				int_upper = (abs(ring_vs[i] - ring_vs[j]) * ring_sigma[i][:-1] * ring_sigma[j][:-1] * cell_vol).sum()
				int_lower = (ring_sigma[i][:-1] * ring_sigma[j][:-1] * cell_vol).sum()

				wavg = int_upper/int_lower
				
				print(f"V{i+1}_{j+1} = {round(wavg, 3)} m/s")

				v_wavg.append(wavg)

	return v_wavg


def vdiff_gas(ring_vs, ring_sigma, ring_gasv, ring_gassig, cell_vol, n_species, dt, path, thick=True):

	gasdiff = []
	print("\n \n")

	for i in range(len(ring_vs)):

		if thick:

			int_upper = (abs(ring_vs[i] - ring_gasv) * ring_sigma[i][:, :-1] * ring_gassig[0][:, :-1] * cell_vol).sum()
			int_lower = (ring_sigma[i][:, :-1] * ring_gassig[0][:, :-1] * cell_vol).sum()

			wavg = int_upper/int_lower

			print(f"V{i+1}_gas = {round(wavg, 3)} m/s")

			gasdiff.append(wavg)

		else:

			int_upper = (abs(ring_vs[i] - ring_gasv) * ring_sigma[i][:-1] * ring_gassig[:][:-1] * cell_vol).sum()
			int_lower = (ring_sigma[i][:-1] * ring_gassig[:][:-1] * cell_vol).sum()

			wavg = int_upper/int_lower

			print(f"V{i+1}_gas = {round(wavg, 3)} m/s")

			gasdiff.append(wavg)

	# Plotting the velocity differences
	labels = []

	for i in range(n_species):

		labels.append([i+1, "gas"])

	# plt.figure(figsize=(15, 6))

	x = range(len(gasdiff))
	plt.xticks(x, labels, rotation=60)
	plt.plot(x, gasdiff, marker='o', label=f"Dust trap at {dt}")
	plt.title(f"Velocity differences between gas and the dust species")
	plt.xlabel("Species")
	plt.ylabel("Velocity (m/s)")
	# plt.grid()
	plt.legend()
	plt.tight_layout()
	plt.savefig(path+f"veldiff_gas_{dt}.png")
	# plt.show()

	return gasdiff 




def plot_vdiff(vdiff, n_species, dts, path):
	""" 
	This function plots the velocity differences for the various sets of dust species for all the dust trapped regions

	Parameters:

	vdiff:		2-D list of velocity differences for the various sets of dust species for all the dust traps
	n_species:	number of dust species (integer)
	dts:		1-D list of the radial indices of all the dust trapped regions
	path:		path of the output directory of the plot
	"""

	labels = []

	for i in range(n_species):

		for j in range(i+1, n_species):

			labels.append(f"{i+1}-{j+1}")

	plt.figure(figsize=(25, 10))

	x = range(len(vdiff[0]))
	plt.xticks(x, labels, rotation=60)
	for i in range(len(vdiff)):
		plt.plot(x, vdiff[i], marker='o', label=f"Dust trap at {dts[i]}")
		plt.title(f"Velocity differences for {n_species} species")
		plt.xlabel("Species")
		plt.ylabel("Velocity (m/s)")
	plt.grid()
	plt.legend()
	plt.xticks(fontsize=6)
	plt.tight_layout()
	plt.savefig(path+"vel_diff.png")
	plt.show()



#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------		

if __name__ == "__main__":

	path = "../fargo3d/outputs/fargo_multifluid/"
	pic_path = "temp_pics/"
	output_number = 7
	species_number = 15
	p_bumps = 1
	vdiffs = []

	r, phi, sigma, vel, energy = load_dust_outputs(path, species_number, output_number)
	gas_sig, gas_vel, gas_energy = load_gas_outputs(r, phi, path, output_number)

	dt = get_dust_trap(sigma, pic_path, species_number, p_bumps=p_bumps)
	plt.clf()

	# Obtaining velocity differences at each dust trapped region

	for i in range(len(dt)):

		volume = cell_volume(r, phi, dt[i])

		v_ring, sigma_ring, vgas_ring, sigmagas_ring = get_ring_param(sigma, vel, gas_sig, gas_vel, dt[i])

		vdiffs.append(vdiff_wavg(v_ring, sigma_ring, volume))
		gasdiffs = vdiff_gas(v_ring, sigma_ring, vgas_ring, sigmagas_ring, volume, species_number, dt[i], pic_path)

		print(f"Average velocity difference for dust trapped region at {dt[i]} : {np.average(vdiffs)}")
		print(f"Average velocity difference between gas and dust at {dt[i]} : {np.average(gasdiffs)}")

	plot_vdiff(vdiffs, species_number, dt, pic_path)

	

