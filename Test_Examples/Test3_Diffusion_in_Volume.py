
# simulating the diffusion in volumes

import math
import numpy
import pylab
import random
import time

import steps.model as smodel

import steps.solver as solvmod
import steps.geom as stetmesh
import steps.rng as srng
import steps.utilities.meshio as smeshio

# Number of iterations to run
NITER = 100

# Data collection time increment
DT = 0.001

# Simulation end time
INT = 0.101

# Number of molecules to be injected into the centre
NINJECT = 10000

# Number of tetrahedral elements to sample data from
SAMPLE = 2000

# Diffusion constant (m^2/s)
DCST = 20.00e-12

# Array to hold tetrahedron indices (integers)
tetidxs = numpy.zeros(SAMPLE, dtype = 'int')

# Array to hold tetrahedron radial distances (floats)
tetrads = numpy.zeros(SAMPLE)


# Model Specification

def gen_model():

    mdl = smodel.Model()
    A = smodel.Spec('A', mdl)
    vsys1 = smodel.Volsys('cytosolv', mdl)
    vsys2 = smodel.Volsys('nucleusv', mdl)
    diff_A_1 = smodel.Diff('diff_A_1', vsys1, A, DCST)
    diff_A_2 = smodel.Diff('diff_A_2', vsys2, A, DCST)

    return mdl


def gen_geom():
    mesh = smeshio.loadMesh('testMeshSphere')[0]

    # Total no. of Tetrahedron in the mesh
    ntets = mesh.countTets()

    tets_comp1 = []
    tets_comp2 = []

    tris_comp1 = set()
    tris_comp2 = set()

    # Crestion of a compartment object
    comp = stetmesh.TmComp('cyto',mesh, range(ntets))
    comp.addVolsys('cytosolv')

    # Central tetrahedron index
    ctetidx = mesh.findTetByPoint([0, 0, 0])
    tetidxs[0] = ctetidx

    # Central tetrahedron's four neighbours
    neighbidcs = mesh.getTetTetNeighb(ctetidx)
    tetidxs[1], tetidxs[2], tetidxs[3], tetidxs[4] = neighbidcs

    stored = 5

    max = mesh.getBoundMax()
    min = mesh.getBoundMin()

    while(stored<SAMPLE):
        rnx = random.random()
        rny = random.random()
        rnz = random.random()

        xcrd = min[0] + (max[0] - min[0]) * rnx
        ycrd = min[1] + (max[1] - min[1]) * rny
        zcrd = min[2] + (max[2] - min[2]) * rnz

        tidx = mesh.findTetByPoint([xcrd, ycrd, zcrd])

        if tidx == -1:
            continue

        if tidx not in tetidxs:
            tetidxs[stored] = tidx
            stored += 1

        cbaryc = mesh.getTetBarycenter(ctetidx)

        for i in range(SAMPLE):
            baryc  = mesh.getTetBarycenter(tetidxs[i])
            r = math.sqrt(math.pow((baryc[0]-cbaryc[0]),2)+math.pow((baryc[1]-cbaryc[1]),2)+math.pow((baryc[2]-cbaryc[2]),2))

            tetrads[i] = r*1.0e6

    R = tetrads.getBoundMax()
    print tetrads.getBoundMax()

    for j in range(ntets):
        baryc = mesh.getTetBarycenter(j)
        tris = mesh.getTetTriNeighb(j)

        r = math.sqrt(math.pow((baryc[0] - cbaryc[0]), 2) + math.pow((baryc[1] - cbaryc[1]), 2) + math.pow((baryc[2] - cbaryc[2]), 2))

        if r > 0.01*R:
            tets_comp1.append(j)
            tris_comp1.add(tris[0])
            tris_comp1.add(tris[1])
            tris_comp1.add(tris[2])
            tris_comp1.add(tris[3])

        else:
            tets_comp2.append(j)
            tris_comp2.add(tris[0])
            tris_comp2.add(tris[1])
            tris_comp2.add(tris[2])
            tris_comp2.add(tris[3])

    comp1 = stetmesh.TmComp('comp1', mesh, tets_comp1)
    comp2 = stetmesh.TmComp('comp2', mesh, tets_comp2)

    comp1.addVolsys('vsys1')
    comp2.addVolsys('vsys2')

    tris_DB = tris_comp1.intersection(tris_comp2)
    tris_DB = list(tris_DB)
    diffb = stetmesh.DiffBoundary('diffb', mesh, tris_DB)

    return mesh

if __name__ == "__main__":
    gen_geom()
    gen_model()