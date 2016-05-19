#### PROJECT MSI 1: MEMBRANE ADDITION IN PDB INTERMEMBRANE PROTEIN ####

#############   Joan Francesc Gilabert & Andreu Bofill   ##############


import htmd as ht
import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D


def dummy_leaflet(pmin, pmax, z, spacing=3.0, name='O'):
    """ This function allows to create the membrane with slab of dummy atoms """
    # Suposing that pmax and pmin are numpy arrays of length 2
    sizes = (pmax - pmin)/spacing + 1
    (Ni, Nj) = sizes.astype(int)
    mem = ht.Molecule()
    mem.empty(Ni*Nj)
    mem.set('resname', 'DUM')
    mem.set('record', 'HETATM')
    mem.resid = np.arange(Ni*Nj)
    mem.set('name', name, 'name UNK')
    x = pmin[0]
    y = pmin[1]
    for i in range(Ni*Nj):
        if y == Nj*spacing + pmin[1]:
            y = pmin[1]
            x += spacing
        mem.set('coords', np.array([x, y, z]), 'resid ' + str(i))
        y += spacing

    return mem


def scatter3D(PDB, selection):
    """Creates a scatter plot in 3D, not currently used"""
    mol = ht.Molecule(PDB)
    cord = mol.get('coords', selection)
    zs = np.array(cord[:, 2])
    ys = np.array(cord[:, 1])
    xs = np.array(cord[:, 0])
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(xs, ys, zs, 'ok')
    plt.show()


def create_planes(mol, zlim=10000.0):
    """Calculates the mean of the normal vectors to the aromatic rings"""
    sel = mol.atomselect("aromatic and z < {:f}".format(zlim))
    res = np.unique(mol.get('resid', sel))
    n = len(res)
    if n == 0:
        raise ValueError('No aromatic residues found')
    vec = np.array([0.0, 0.0, 0.0])
    for r in res:
        names = mol.get('name', 'resid ' + str(r))
        cords = np.zeros((3, 3))
        i = 0
        for nam in names[5:8]:
            temp = mol.get('coords', 'resid ' + str(r) +
                           ' and name ' + nam)
            if temp.ndim > 1:
                cords[i, :] = np.mean(temp, 0)
            else:
                cords[i, :] = temp
            i += 1
        vec += np.cross(cords[0, :]-cords[1, :], cords[0, :]-cords[2, :])
    vec /= n
    return vec


def histog(mol, bins, plot=True, code='PDB'):
    """
        This function return a set of plots about the distribution
        of hydrophobic and charged residues through a set of z-axis segments
    """
    chcord = mol.get('coords', 'charged')
    hydcord = mol.get('coords', 'hydrophobic')
    chz = chcord[:, 2]
    hydz = hydcord[:, 2]
    histch = np.histogram(chz, bins=bins)
    histhyd = np.histogram(hydz, bins=bins)
    histtot = np.histogram(co[:, 2], bins=bins)
    # histdiff is the "density" of hydrophobic residues
    histdiff = histhyd[0] / histtot[0]
    if not plot:
        return histdiff
    center = (bins[:-1] + bins[1:]) / 2
    fig = plt.figure()
    plt.plot(center, histdiff)
    plt.title('Normalized hydrophobic frequency')
    plt.xlabel('Z($\AA$)')
    fig.savefig('plots/norm_hyd_{:s}.png'.format(code),
                bbox_inches='tight', dpi=300)
    fig.show()
    fig = plt.figure()
    plt.plot(center, histtot[0])
    plt.plot(center, histhyd[0])
    plt.plot(center, histch[0])
    plt.title('Frequency of different residues')
    plt.legend(['Total frequency', 'Hydrophobic frequency',
                'Charged frequency'], loc='best')
    plt.xlabel('Z($\AA$)')
    plt.savefig('plots/comphist_{:s}.png'.format(code),
                bbox_inches='tight', dpi=300)
    fig.show()
    return histdiff


def find_z(values, threshold):
    """
    To predict the start and end position of the intramembrane atoms
        on the molecule
    """
    streak = 0
    start = 0
    end = 0
    index = 0
    for value in windowed > threshold:
        if value:
            if streak == 0:
                start = index
            streak += 1
        else:
            if streak > 3:
                end = index
                break
            streak = 0
        index += 1
    return (start, end)


def plot_vector(vec):
    """Draws a vector as a line with an arrow, currently not used """
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    orig = np.array([0, 0, 0])
    ax.quiver(np.array(orig[0]), np.array(orig[1]), np.array(orig[2]),
              np.array(vec[0]), np.array(vec[1]), np.array(vec[2]*5),
              length=0.1, color="Tomato", pivot="tail")
    ax.view_init(elev=18, azim=30)
    plt.show()

if __name__ == '__main__':
    try:
        pdb = sys.argv[1]
    except:
        ValueError("Pass the name of the pdb to analyze as an argument to the \
                   program(the PDB file should be in the same folder")
    if pdb.find('.pdb') == -1:
        mol = ht.Molecule(pdb)
    else:
        mol = ht.Molecule("pdb/{:s}".format(pdb))
    filename = mol.viewname.split('.')
    PDBcode = filename[0]
    # Plots is boolean that controls wether to print several plots
    plots = True
    # Block to refine rotation, it only performs one iteration, because it does
    # not seem to converge
    vec_mean = create_planes(mol)
    alpha = np.arccos(vec_mean[2]/np.linalg.norm(vec_mean))
    axis = np.cross(vec_mean, np.array([0, 0, 1]))
    rot = ht.rotationMatrix(axis, alpha)
    mol.rotateBy(rot)
    co = mol.get('coords')
    c = np.mean(co, axis=0)
    vec_res = np.dot(rot, vec_mean)
    # This vector is calculated to ensure that proper rotation (parallel to z
    # axis)
    assert np.sum(vec_res[0:1]) < 1e-16, "Incorrect rotaion, vector not \
                parallel to z axis"
    zmin = min(co[:, 2])
    zmax = max(co[:, 2])
    pmin = np.array([min(co[:, 0])-10, min(co[:, 1])-10])
    pmax = np.array([max(co[:, 0])+10, max(co[:, 1])+10])

    zrange = zmax - zmin
    step = 2
    nsteps = int(zrange/step)
    bins = np.linspace(zmin, zmax, num=nsteps)
    center = (bins[:-1] + bins[1:]) / 2
    histdiff = histog(mol, bins, plot=plots, code=PDBcode)
    # Apply window means to soften the curve
    windowed = histdiff[0:-2]+histdiff[1:-1]+histdiff[2:]
    windowed /= 3
    threshold = np.mean(windowed)

    if plots:
        fig = plt.figure()
        plt.plot(center[1:-1], windowed)
        plt.plot((center[1], center[-2]), (threshold, threshold), 'k-')
        plt.title('Normalized hydroph. freq. windowed')
        plt.xlabel('Z($\AA$)')
        plt.legend(['Windowed frequency', 'Mean frequency'])
        fig.savefig("plots/windowed_{:s}.png".format(PDBcode),
                    bbox_inches='tight', dpi=300)
        fig.show()

    (start, end) = find_z(windowed, threshold)
    # Index is off by 1 because we have used a window previously
    z1 = center[start+1]
    z2 = center[end+1]
    mem1 = dummy_leaflet(pmin, pmax, z1)
    mem2 = dummy_leaflet(pmin, pmax, z2, name='N')
    mol.append(mem1)
    mol.append(mem2)
    mol.view(style='Licorice')
    mol.write('pdb/'+PDBcode+'_mem.pdb')
    sys.exit()
