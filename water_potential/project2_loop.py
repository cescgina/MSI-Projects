# PROJECT MSI 2: Water chemical potential ####

#   Joan Francesc Gilabert  ##############


import htmd as ht
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def scatter3D(cord):
    zs = np.array(cord[:, 2])
    ys = np.array(cord[:, 1])
    xs = np.array(cord[:, 0])
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(xs, ys, zs, 'ok')
    plt.show()

if __name__ == '__main__':
    mol = ht.Molecule('data/1x14/structure.pdb')
    mol.read('data/1x14/traj.xtc')
    mol.wrap('protein')
    mol.align('protein')

    (natom, ndim, nframe) = mol.coords.shape
    mask = mol.get('serial', 'water and name OH2')
    traj_O = mol.coords[mask-1, :, :]
    nwat = traj_O.shape[0]
    # Reshape the coordinates matrix in order to have a 2 dimension array which
    # can be used to calculate the histogram
    traj_flat = traj_O.transpose(2, 0, 1).reshape(-1, traj_O.shape[1])
    maxval = traj_flat.max(axis=0)
    maxval = maxval[maxval < 100].max()
    minval = traj_flat.min(axis=0)
    minval = minval[minval > -100].max()
    maxvec = np.array([maxval+15]*3)
    minvec = np.array([minval-15]*3)
    res = np.array([1, 1, 1])
    bins = np.ceil((maxvec-minvec)/res)
    # Tuple of ints needed to pass to histogram
    bins = tuple(bins.astype(int))
    edges = [np.arange(minvec[0], maxvec[0], res[0])]*3
    H_loop = np.zeros(bins)
    traj_flat += abs(minvec[0])
    traj_flat = np.mod(traj_flat, maxvec[0]+abs(minvec[0]))
    traj_flat = np.rint(traj_flat).astype(int)
    for row in traj_flat:
        H_loop[row[0], row[1], row[2]] += 1

    Kb = 0.001987191  # kcal/mol/K
    T = 298  # K, value not confirmed
    H_loop /= (nwat*nframe)  # Normalize so that it corresponds to a probability
    H_loop += 1e-20  # Add a very small quantity to avoid taking the log of 0
    G_loop = - Kb*T*np.log(H_loop)
    ht.molecule.util.writeVoxels(G_loop, 'voxel_loop.cube', minvec, maxvec, res)
    sys.exit()
