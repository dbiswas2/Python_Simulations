
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as solvmod
import steps.utilities.meshio as  meshio
import numpy as np
import pylab
# import pyplot
import time

start_time = time.time()


# Defining the Constituents:
def gen_model():


    mdl = smodel.Model()

    U = smodel.Spec('U', mdl)
    V = smodel.Spec('V', mdl)

    vsysA = smodel.Volsys('vsysA', mdl)
    vsysB = smodel.Volsys('vsysB', mdl)

    dcst_U_A = 0.1e-9
    dcst_U_B = 0.05e-9
    dcst_V_A = 0.05e-9
    dcst_V_B = 0.1e-9

    diff_U_A = smodel.Diff('diff_U_A', vsysA, U, dcst_U_A)
    diff_U_B = smodel.Diff('diff_U_B', vsysB, U, dcst_U_B)
    diff_V_A = smodel.Diff('diff_V_A', vsysA, V, dcst_V_A)
    diff_V_B = smodel.Diff('diff_V_B', vsysB, V, dcst_V_B)


    return mdl

# Defining the Geometry:
def gen_geom():

    #n = 5
    #print n

    mesh = meshio.loadMesh('testMeshCylinder1')[0]

    ntets = mesh.countTets()

    tets_compA = []
    tets_compB = []
    tris_compA = set()
    tris_compB = set()

    z_max = mesh.getBoundMax()[2]
    z_min = mesh.getBoundMin()[2]
    x_max = mesh.getBoundMax()[0]
    x_min = mesh.getBoundMin()[0]
    y_max = mesh.getBoundMax()[1]
    y_min = mesh.getBoundMin()[1]

    print z_max
    print z_min
    print x_max
    print x_min
    print y_max
    print y_min

    z_mid = z_min + (z_max - z_min) / 2.0

    for t in range(ntets):

        barycz = mesh.getTetBarycenter(t)[2]

        tris = mesh.getTetTriNeighb(t)

        if (barycz < z_mid):
            tets_compA.append(t)
            tris_compA.add(tris[0])
            tris_compA.add(tris[1])
            tris_compA.add(tris[2])
            tris_compA.add(tris[3])
        else:
            tets_compB.append(t)
            tris_compB.add(tris[0])
            tris_compB.add(tris[1])
            tris_compB.add(tris[2])
            tris_compB.add(tris[3])

    compA = sgeom.TmComp('compA', mesh, tets_compA)
    compB = sgeom.TmComp('compB', mesh, tets_compB)

    compA.addVolsys('vsysA')
    compB.addVolsys('vsysB')

    tris_DB = tris_compA.intersection(tris_compB)
    tris_DB = list(tris_DB)

    diffusion_boundary = sgeom.DiffBoundary('diffusion_boundary', mesh, tris_DB)

    return mesh, tets_compA, tets_compB


# print gen_geom()

mdl = gen_model()

mesh, tets_compA, tets_compB = gen_geom()

rng = srng.create('mt19937', 256)
rng.initialize(654)

sim = solvmod.Tetexact(mdl, mesh, rng)
sim.reset()

tpnts = np.arange(0.0, 0.101, 0.001)
ntpnts = tpnts.shape[0]

ntets = mesh.countTets()
print 'ntets', ntets


resU = np.zeros((ntpnts, ntets))
resV = np.zeros((ntpnts, ntets))

tetU = mesh.findTetByPoint([0, 0, -4.99e-10])
tetV = mesh.findTetByPoint([0, 0, 4.99e-10])

# print 'tetU', tetU
# print 'tetV', tetV

sim.setTetCount(tetU , 'U', 1000)
sim.setTetCount(tetV, 'V', 500)


sim.setDiffBoundaryDiffusionActive('diffusion_boundary', 'U',True)
sim.setDiffBoundaryDiffusionActive('diffusion_boundary', 'V',True)


print ntpnts

for i in range(ntpnts):
    sim.run(tpnts[i])
    print i

    for k in range(ntets):
        resU[i,k] = sim.getTetCount(k, 'U')
        resV[i,k] = sim.getTetCount(k, 'V')



# Plotting the Data

def plot_binned(t_idx, bin_n):

        if (t_idx > tpnts.size):
            print "Time index is out of range."
            return


        z_tets = np.zeros(ntets)

        zbound_min = mesh.getBoundMin()[2]


        for i in range(ntets):
            baryc = mesh.getTetBarycenter(i)
            z = baryc[2] - zbound_min

            z_tets[i] = z * 1.0e6

        z_max = z_tets.max()
        z_min = z_tets.min()


        z_seg = (z_max - z_min) / bin_n
        bin_mins = np.zeros(bin_n + 1)
        z_tets_binned = np.zeros(bin_n)
        bin_vols = np.zeros(bin_n)


        z = z_min
        for b in range(bin_n + 1):
            bin_mins[b] = z
            if (b != bin_n): z_tets_binned[b] = z + z_seg / 2.0
            z += z_seg
        bin_counts = [None] * bin_n
        for i in range(bin_n): bin_counts[i] = []
        for i in range((resU[t_idx].size)):
            i_z = z_tets[i]
            for b in xrange(bin_n):
                if (i_z >= bin_mins[b] and i_z < bin_mins[b + 1]):
                    bin_counts[b].append(resU[t_idx][i])
                    bin_vols[b] += sim.getTetVol(i)
                    break

        # Converting concentration in arbitrary units
        bin_concs = np.zeros(bin_n)
        for c in range(bin_n):
            for d in range(bin_counts[c].__len__()):
                bin_concs[c] += bin_counts[c][d]
            bin_concs[c] /= (bin_vols[c] * 1.0e24)

        t = tpnts[t_idx]

        # Plot the data
        pylab.scatter(z_tets_binned, bin_concs, label='X', color='black')


        z = z_min
        for b in range(bin_n + 1):
            bin_mins[b] = z
            if (b != bin_n): z_tets_binned[b] = z + z_seg / 2.0
            z += z_seg
        bin_counts = [None] * bin_n
        for i in range(bin_n): bin_counts[i] = []
        for i in range((resV[t_idx].size)):
            i_z = z_tets[i]
            for b in xrange(bin_n):
                if (i_z >= bin_mins[b] and i_z < bin_mins[b + 1]):
                    bin_counts[b].append(resV[t_idx][i])
                    break
        bin_concs = np.zeros(bin_n)
        for c in range(bin_n):
            for d in range(bin_counts[c].__len__()):
                bin_concs[c] += bin_counts[c][d]
            bin_concs[c] /= (bin_vols[c] * 1.0e24)

        pylab.scatter(z_tets_binned, bin_concs, label='Y', color='red')


        pylab.xlabel('Z axis', fontsize=16)
        pylab.ylabel('Bin concentration (N/m^3)', fontsize=16)
        pylab.ylim(0)
        pylab.xlim(0, 0.10)
        pylab.legend(numpoints=1)
        pylab.show()
        return
print("--- %s seconds ---" % (time.time() - start_time))
print plot_binned(100, 50)
