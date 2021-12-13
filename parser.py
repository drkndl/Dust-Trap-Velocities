import argparse 

def parse():

	parser = argparse.ArgumentParser(description = 'Calculates the relative velocities of gas and dust in the dust trap of a protodisk')

	parser.add_argument("-p", "--path", help = "Path of the simulation output files", default = "../fargo3d/outputs/fargo_multifluid/")
	parser.add_argument("-pp", "--pic_path", help = "Path where the plots are to be saved", default = "temp_pics/")
	parser.add_argument("-o", "--output", help = "Number to be considered for the simulation output files", default = 7, type = int)
	parser.add_argument("-s", "--species", help = "Total number of dust species", default = 15, type = int)
	parser.add_argument("-b", "--p_bumps", help = "Number of dust traps to be considered", default = 1, type = int)
	parser.add_argument("--rin", help = "Inner simulation radius", default = 0.4, type = float)
	parser.add_argument("--rout", help = "Outer simulation radius", default = 2.5, type = float)
	parser.add_argument("--ny", help = "Number of radial zones", default = 128, type = int)

	args = parser.parse_args()

	return args.path, args.pic_path, args.output, args.species, args.p_bumps, args.rin, args.rout, args.ny
