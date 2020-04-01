import MDAnalysis as md
import sys
import numpy as np

if __name__=='__main__':
	psf = sys.argv[1]
	traj = sys.argv[2]
	ofname = sys.argv[3]
	U = md.Universe(psf,traj)
	sel = U.select_atoms('protein')

	U.trajectory.rewind()
	rog = [sel.radius_of_gyration() for ts in U.trajectory]

	rog = np.array(rog)
	np.savetxt(ofname, rog)
