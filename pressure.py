import numpy as np 
import matplotlib.pyplot as plt 
import os

def pressure(path, n_species, n_out):

	# load cell interface locations
	r = np.loadtxt(os.path.join(path, "domain_y.dat"))[3:-3] # ignore ghost cells
	phi = np.loadtxt(os.path.join(path, "domain_x.dat"))

	r_m, phi_m = np.meshgrid(r, phi)		# Their shape is (385, 129)

	for i in range(len(r)):
		vol

	# vol = np.pi*(r_m**2 - r_m**2)
	# print(vol)
	# print(np.shape(vol))

path = "./Set 4/Output 7/fargo_multifluid/"
output_number = 175
species_number = 10

pressure(path, species_number, output_number)