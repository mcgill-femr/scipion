from prody import *
import matplotlib.pyplot as plt

#ANM Calculations
GluA2_sim_ca = parsePDB('/home/javiermota/ProDy/initial_ionized.pdb', subset='ca')
GluA2_em_ca = parseCIF('4uqj', subset='ca')
anm = ANM('GluA2 AMPAR sim model')
anm.buildHessian(GluA2_sim_ca)
anm.calcModes()
slowest_mode = anm[0]
print( slowest_mode.getEigval().round(3) )
abs(anm[0] * anm[1]).round(10)
writeNMD('ini_ampar_anm.nmd', anm, GluA2_sim_ca)
saveModel(anm, 'ampar')
anm = loadModel('ampar.anm.npz')

#PCA Calculations
GluA2_sim = parsePDB('/home/javiermota/ProDy/initial_ionized.pdb')
combined_traj = Trajectory('/home/javiermota/ProDy/initr.dcd')
combined_traj.setCoords(GluA2_sim)
combined_traj.setAtoms(GluA2_sim.ca)
combined_traj.addFile('/home/javiermota/ProDy/fintr.dcd')
pca = PCA('AMPAR trajectories')
pca.buildCovariance(combined_traj)
pca.calcModes()

for mode in pca[:4]:
    print(calcFractVariance(mode).round(2))

ini_traj = Trajectory('/home/javiermota/ProDy/initr.dcd')
ini_traj.setCoords(GluA2_sim)   # Set the initial structure as the reference
ini_traj.setAtoms(GluA2_sim.ca) # A shortcut for .select('ca')

fin_traj = Trajectory('/home/javiermota/ProDy/fintr.dcd')
fin_traj.setCoords(GluA2_sim)   # Set the initial structure as the reference
fin_traj.setAtoms(GluA2_sim.ca) # A shortcut for .select('ca')

showProjection(ini_traj, pca[:2], color='g', new_fig=True)
showProjection(fin_traj, pca[:2], color='r', new_fig=False)

showProjection(ini_traj, anm[:3], color='g', new_fig=True)
showProjection(fin_traj, anm[:3], color='r', new_fig=False)

writeNMD('combined_pca.nmd', pca, GluA2_sim_ca)
saveModel(pca, 'combined', matrices=True) # By default we don't save the covariance matrix
pca = loadModel('combined.pca.npz')
showOverlapTable(pca, anm)
printOverlapTable(pca[:7], anm[:7])
overlaps = calcOverlap(pca, anm)