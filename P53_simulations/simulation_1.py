
import steps.model as smodel
import steps.geom as sgeom
import steps.rng as srng
import steps.solver as solvmod
import steps.utilities.meshio as  smeshio
import numpy as np
import pylab
import math
# import time
#
# start_time = time.time()


def Lambda(X):
    return math.exp(X / (1 + X))


def gen_model():


    mdl = smodel.Model()

    A_c = smodel.Spec('A_c', mdl)
    B_c = smodel.Spec('B_c', mdl)
    C_c = smodel.Spec('C_c', mdl)

    A_n = smodel.Spec('A_n', mdl)
    B_n = smodel.Spec('B_n', mdl)
    C_n = smodel.Spec('C_n', mdl)



    vsysX = smodel.Volsys('vsysX', mdl)
    vsysY = smodel.Volsys('vsysY', mdl)

    #
    # lambda B_c: math.exp(B_c / (1 + B_c))
    # lambda B_n: math.exp(B_n / (1 + B_n))

    R_Ac_f = smodel.SReac('R_Ac_f', vsysX, lhs = [C_c], rhs = [A_c], k=Lambda(B_c))
    R_An_f = smodel.SReac('R_An_f', vsysY, lhs = [C_n], rhs=  [0.1, A_n], k=Lambda(B_n))

    R_Bc_f = smodel.SReac('R_Bn_f', vsysX, lhs = [B_c], rhs = [0], k = 0.1)
    R_Bn_f = smodel.SReac('R_Bn_f', vsysX, lhs = [B_n], rhs = [0], k = 0.1)

    R_Cc_f = smodel.SReac('R_Bn_f', vsysX, lhs=[C_c], rhs=[0], k = 0.1)
    R_Cn_f = smodel.SReac('R_Bn_f', vsysX, lhs=[C_n], rhs=[0], k = 0.1)

    dcst_Ac_X = 0.1e-9
    dcst_An_Y = 0.1e-9
    dcst_Bc_X = 0.1e-9
    dcst_Bn_Y = 0.1e-9



    diff_Ac_X = smodel.Diff('diff_Ac_x', vsysX, A_c, dcst_Ac_X)
    diff_An_Y = smodel.Diff('diff_An_y', vsysY, A_n, dcst_An_Y)
    diff_Bc_X = smodel.Diff('diff_Bc_x', vsysX, B_c, dcst_Bc_X)
    diff_Bn_Y = smodel.Diff('diff_Bn_y', vsysX, B_n, dcst_Bn_Y)


    return mdl


def gen_geom():
    mesh = smeshio.loadMesh('testMeshSphere')[0]

    # Total no. of Tetrahedron in the mesh
    ntets = mesh.countTets()

    tets_comp1 = []