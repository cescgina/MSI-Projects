# PROJECT MSI 2: Water chemical potential ####

#   Joan Francesc Gilabert  ##############


import htmd as ht
import numpy as np
import glob as gb
import scipy.ndimage as nd
import sys
import copy


class Simulation():
    """The Simulation class loads a folder with trajectory data and creates
    an HTMD molecule that contains the pdb and the xtc trajectory.
    Additionally, it also has attributes needed to calculate the histogram of
    water molecules
    """
    def __init__(self, folder, ref_sim, ref_mol):
        """The init method loads the pdb and trajectory files into an HTMD
        molecule, and prepares the molecule for the posterior work (align and
        wrap, see HTMD documentation for a description of the methods)"""
        self.mol = ht.Molecule('{:s}/structure.pdb'.format(folder))
        self.mol.read('{:s}/traj.xtc'.format(folder))
        self.mol.wrap('protein')
        self.mol.align('protein', refmol=ref_mol)
        self.maxvec = ref_sim.maxvec
        self.minvec = ref_sim.minvec
        self.bins = ref_sim.bins
        self.edges = ref_sim.edges
        self.res = ref_sim.res

    def prepare_coords(self):
        """Filters the molecule object and extracts its water coordinates to
        calculate the histogram"""
        self.mol.filter('water and name OH2')
        traj_O = self.mol.coords
        # Reshape the coordinates matrix in order to have a 2 dimension array
        # which can be used to calculate the histogram
        self.traj_flat = traj_O.transpose(2, 0, 1).reshape(-1, traj_O.shape[1])

    def water_hist(self):
        """Calculates the histogram of water molecules to latter use it to
        estimate the free energy"""
        self.prepare_coords()
        self.H, foo = np.histogramdd(self.traj_flat, bins=self.edges)


class SimulationRef(Simulation):
    """The SimulationRef class is a special case of the Simulation class which
    processes the first folder, considered as reference"""
    def __init__(self, folder, res):
        self.mol = ht.Molecule('{:s}/structure.pdb'.format(folder))
        self.mol.read('{:s}/traj.xtc'.format(folder))
        self.mol.wrap('protein')
        self.mol.align('protein')
        self.res = res

    def water_hist(self):
        """Calculates the histogram of water molecules to latter use it to
        estimate the free energy"""
        self.prepare_coords()
        maxval = self.traj_flat.max(axis=0)
        maxval = maxval[maxval < 100].max()
        minval = self.traj_flat.min(axis=0)
        minval = minval[minval > -100].max()
        self.maxvec = np.array([maxval+15]*3)
        self.minvec = np.array([minval-15]*3)
        self.bins = np.ceil((self.maxvec-self.minvec)/self.res)
        # Tuple of ints needed to pass to histogram
        self.bins = tuple(self.bins.astype(int))
        self.edges = [np.arange(self.minvec[0]-self.res[0]/2,
                                self.maxvec[0]+self.res[0]/2,
                                self.res[0])]*3
        self.H, foo = np.histogramdd(self.traj_flat, bins=self.edges)


def add_gaussian_noise(H):
    """Adds gaussian noise to the nonzero elements of the array resulting
        from the histogram. This is old code currently not used
    """
    ind = np.array(np.nonzero(H))
    npoints = ind.shape[1]
    size_gauss = 3**3
    for i in range(npoints):
        occup_value = int(H[tuple(ind[:, i])])
        temp = abs(np.random.normal(0, 1.5, (occup_value, size_gauss)))
        temp = temp.mean(axis=0)
        temp = temp.reshape((-1, 3, 3))
        assert temp.shape == (3, 3, 3), "Reshape wrongly broadcasted"
        temp[1, 1, 1] = 0
        loc_ind = ind[:, i]
        try:
            H[loc_ind[0]-1:loc_ind[0]+2, loc_ind[1]-1:loc_ind[1]+2,
                loc_ind[2]-1:loc_ind[2]+2] += temp
        except IndexError:
            # TODO: Add catching errors for borders
            # For now just catches the exception by ignoring the problematic
            # point
            continue
    return H

if __name__ == '__main__':

    datalist = gb.glob('data/*')
    nfolders = len(datalist)
    # Pop the first folder that will be considered the reference
    ref_folder = datalist.pop(0)
    print('Procesing folder {:s}'.format(ref_folder))
    res = np.array([1, 1, 1])  # Resolution of the grid (1 amstrong)
    sim_ref = SimulationRef(ref_folder, res)
    ref_mol = copy.deepcopy(sim_ref.mol)
    sim_ref.water_hist()
    H = sim_ref.H
    for folder in datalist:
        # Iterate over all the folder containing trajectory data, counting now
        # with a reference structure and histogram data (bins edges)
        print('Procesing folder {:s}'.format(folder))
        sim = Simulation(folder, sim_ref, ref_mol)
        sim.water_hist()
        H += sim.H
    Kb = 0.001987191  # kcal/mol/K
    T = 298  # K, value not confirmed
    # Apply a gaussian filter to the histogram data to have the effect of
    # interpolating the water occupancy
    H = nd.gaussian_filter(H, 1.5)
    H /= H.sum()  # Normalize so that it corresponds to a probability
    H += 1e-40  # Add a very small quantity to avoid taking the log of 0
    G = -Kb*T*np.log(H)
    ht.molecule.util.writeVoxels(G, 'voxel.cube', sim_ref.minvec,
                                 sim_ref.maxvec, res)
    sys.exit()
