import numpy as np 
import matplotlib.pyplot as plt 
from find_ring import load_dust_outputs, load_gas_outputs, vdiff_wavg, get_dust_trap, get_ring_param, cell_volume
from scipy.spatial.distance import squareform


def vel_hist(vdiffs, dt, path):
	"""
	This function plots a histogram of the relative velocities between the various dust species

	Parameters:

	vdiffs: 	the velocity differences calculated by vdiff_wavg (1D list)
	dt:			dust trap value (int)
	path: 		the path to where the histogram plot is to be saved

	Does not return anything, but plots the histogram
	"""

	# Calculating the maximum relative velocity value to determine the bins for the histogram
	maximum = 5 + int(np.max(vdiffs))
	bins = np.arange(0, maximum, 5)

	# Plotting the histogram
	arr = plt.hist(vdiffs, bins, edgecolor='black', linewidth=1.2)

	# Above the bars in the plot, the count value for the particular bin is printed if it is not zero
	for i in range(len(bins)-1):

		if arr[0][i] != 0:

			# Improving the formatting 
			if arr[0][i] < 10:
				plt.text(arr[1][i] + 1, arr[0][i] + 1, str(int(arr[0][i])))
			else:
				plt.text(arr[1][i], arr[0][i] + 1, str(int(arr[0][i])))

	plt.title(f"Distribution of velocity differences for dust trap at {dt}")
	plt.xlabel("Velocity (m/s)")
	plt.ylabel("Frequency")
	plt.savefig(path + f"hist_{dt}.png")
	plt.show()


def vel_contour(vdiffs, dt, n_species, path):
	"""
	This function creates a contour plot of relative velocities through the varying Stokes numbers of the dust species

	Parameters:

	vdiffs: 	the velocity differences calculated by vdiff_wavg (1D list)
	dt:			dust trap value (int)
	n_species:	the total number of dust species (int)
	path: 		the path to where the histogram plot is to be saved

	Does not return anything, but plots the contour
	"""

	# Converting the list of relative velocities into a symmetric square matrix 
	vel_matrix = squareform(vdiffs)

	# Determining the Stokes numbers and creating a mesh of them
	xarr = np.logspace(-5, 3, n_species)
	X, Y = np.meshgrid(xarr, xarr)

	# Determining the levels of contours
	maximum = 5 + int(np.max(vdiffs))
	ticks1 = np.arange(0, 1, 0.1)
	ticks2 = np.arange(1, maximum, 2)
	ticks = np.concatenate((ticks1, ticks2))

	# Defining the colour map
	colorinterpolation = 50
	colourMap = plt.cm.jet

	fig = plt.figure()

	# Plotting the contour plot on a log-log scale
	plt.xscale('log')
	plt.yscale('log')
	plt.contourf(X, Y, vel_matrix, colorinterpolation, levels=ticks, cmap=colourMap)
	plt.title("Contour Plot of Velocity Differences (m/s)")
	plt.xlabel("Grain Stokes numbers")
	plt.ylabel("Grain stokes numbers")
	plt.colorbar()
	plt.savefig(path+f"velcontour_{dt}.png")
	plt.show()


def main():

	path = "../fargo3d/outputs/fargo_multifluid/" 		# Path of simulation output files
	pic_path = "temp_pics/" 							# Path where the plots are to be saved
	output_number = 8									# Number to be considered for the simulation output files
	species_number = 15 								# Total number of dust species
	p_bumps = 1 										# Number of dust traps to be considered
	r_in = 0.4
	r_out = 2.5
	Ny = 128

	r, phi, sigma, vel, energy = load_dust_outputs(path, species_number, output_number)
	gas_sig, gas_vel, gas_energy = load_gas_outputs(r, phi, path, output_number)

	# Determining location of dust traps
	dt = get_dust_trap(sigma, pic_path, species_number, r_in, r_out, Ny, p_bumps=p_bumps)

	# Clear the radial surface density plot
	plt.clf()

	# Obtaining velocity differences at each dust ring
	for i in range(len(dt)):

		volume = cell_volume(r, phi, dt[i])
		v_ring, sigma_ring, vgas_ring, sigmagas_ring = get_ring_param(sigma, vel, gas_sig, gas_vel, dt[i])
		vdiffs = vdiff_wavg(v_ring, sigma_ring, volume)
		vel_hist(vdiffs, dt[i], pic_path)
		vel_contour(vdiffs, dt[i], species_number, pic_path)


if __name__ == "__main__":

	main()
