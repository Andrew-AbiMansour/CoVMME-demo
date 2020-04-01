import MDAnalysis as md
import sys
import numpy as np

if __name__=='__main__':
	psf = sys.argv[1]
	traj = sys.argv[2]
	ofname = sys.argv[3]
	trajfname = '1ohf.MD.dcd'
	U = md.Universe(psf,traj)
	sel = U.selectAtoms('all')
	writer = md.Writer(trajfname,U.atoms.numberOfAtoms())
	rgyr = []

	U.trajectory.rewind()
	for ts in U.trajectory: 
		print('Calculating ROG from frame {}'.format(U.trajectory.frame))
		rgyr.append(sel.radiusOfGyration())	
		writer.write(sel)

	ROG = np.array(rgyr)
	np.savetxt(ofname,ROG)
		
	
