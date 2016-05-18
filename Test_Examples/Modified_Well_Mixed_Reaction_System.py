import steps.model as smodel
import steps.geom as swm
import steps.rng as srng
import steps.solver as ssolver
import numpy
import pylab

# Model: 3A + B <=> 2C

mdl = smodel.Model()

molA = smodel.Spec('molA', mdl)
molB = smodel.Spec('molB', mdl)
molC = smodel.Spec('molC', mdl)

mol2A = smodel.Spec('mol2A', mdl)
mol3A = smodel.Spec('mol3A', mdl)

volsys = smodel.Volsys('vsys', mdl)

kreac_I1 = smodel.Reac('kreac_I1', volsys, lhs=[molA, molA], rhs=[mol2A], kcst = 0.1e6)
kreac_I2 = smodel.Reac('kreac_I2', volsys, lhs=[mol2A, molA], rhs=[mol3A], kcst = 0.1e6)

kreac_f = smodel.Reac('kreac_f', volsys, lhs=[mol3A, molB], rhs=[molC], kcst = 0.3e6)

print kreac_I1.getOrder()
print kreac_I2.getOrder()

print kreac_f.getOrder()

kreac_b = smodel.Reac('kreac_b', volsys, lhs=[molC], rhs=[molA, molB])
kreac_b.kcst = 0.1e-2

wmgeom = swm.Geom()
comp = swm.Comp('comp', wmgeom)
comp.addVolsys('vsys')
comp.setVol(1.6667e-21)

r = srng.create('mt19937', 256)
r.initialize(23412)

sim = ssolver.Wmdirect(mdl, wmgeom, r)

sim.reset()
sim.setCompConc('comp', 'molA', 31.4e-6)
sim.setCompConc('comp', 'molB', 22.3e-6)

tpnt = numpy.arange(0.0, 2.001, 0.001)
res = numpy.zeros([2001, 5])

for t in range(0, 2001):
    sim.run(tpnt[t])
    res[t, 0] = sim.getCompCount('comp', 'molA')
    res[t, 1] = sim.getCompCount('comp', 'molB')
    res[t, 2] = 0.5*sim.getCompCount('comp', 'molC')
    res[t, 3] = sim.getCompCount('comp', 'mol2A')
    res[t, 4] = sim.getCompCount('comp', 'mol3A')


# Plot number of molecules of 'molA' over the time range:
pylab.plot(tpnt, res[:,0], label = 'A')
# Plot number of molecules of 'molB' over the time range:
pylab.plot(tpnt, res[:,1], label = 'B')
# Plot number of molecules of 'molC' over the time range:
pylab.plot(tpnt, res[:,2], label = 'C')
pylab.plot(tpnt, res[:,3], label = '2A')
pylab.plot(tpnt, res[:,4], label = '3A')

pylab.xlabel('Time (sec)')
pylab.ylabel('#molecules')
pylab.legend()
pylab.show()

# # Plot number of molecules of 'molA' over the time range:
# pylab.plot(tpnt, res[:,0], label = 'A')
# # Plot number of molecules of 'mol2A' over the time range:
# pylab.plot(tpnt, res[:,3], label = '2A')
#
# pylab.xlabel('Time (sec)')
# pylab.ylabel('#molecules')
# pylab.legend()
# pylab.show()

NITER = 100
res = numpy.zeros([NITER, 2001, 5])
tpnt = numpy.arange(0, 2.001, 0.001)

for i in range(0, NITER):
    sim.reset()
    sim.setCompConc('comp', 'molA', 31.4e-6)
    sim.setCompConc('comp', 'molB', 22.3e-6)

    for t in range(0, 2001):
        sim.run(tpnt[t])
        res[i, t, 0] = sim.getCompCount('comp', 'molA')
        res[i, t, 1] = sim.getCompCount('comp', 'molB')
        res[i, t, 2] = sim.getCompCount('comp', 'molC')
        res[i, t, 3] = sim.getCompCount('comp', 'mol2A')
        res[i, t, 4] = sim.getCompCount('comp', 'mol3A')

res_mean = numpy.mean(res, 0)

pylab.plot(tpnt, res_mean[:,0], label = 'A')
pylab.plot(tpnt, res_mean[:,1], label = 'B')
pylab.plot(tpnt, 0.5*res_mean[:,2], label = 'C')
pylab.plot(tpnt, res_mean[:,3], label = '2A')
pylab.plot(tpnt, res_mean[:,4], label = '3A')

pylab.xlabel('Time (sec)')
pylab.ylabel('# molecules')
pylab.legend()
pylab.show()


for i in range(0, NITER):
    sim.reset()
    sim.setCompConc('comp', 'molA', 31.4e-6)
    sim.setCompConc('comp', 'molB', 22.3e-6)

    for t in range(0,1001):
        sim.run(tpnt[t])
        res[i, t, 0] = sim.getCompConc('comp', 'molA')
        res[i, t, 1] = sim.getCompConc('comp', 'molB')
        res[i, t, 2] = 0.5*sim.getCompConc('comp', 'molC')
        res[i, t, 3] = sim.getCompConc('comp', 'mol2A')
        res[i, t, 4] = sim.getCompConc('comp', 'mol3A')

    sim.setCompCount('comp', 'molA', sim.getCompCount('comp', 'molA') + 20)

    for t in range(1001, 2001):
        sim.run(tpnt[t])
        res[i, t, 0] = sim.getCompConc('comp', 'molA')
        res[i, t, 1] = sim.getCompConc('comp', 'molB')
        res[i, t, 2] = 0.5*sim.getCompConc('comp', 'molC')
        res[i, t, 3] = sim.getCompConc('comp', 'mol2A')
        res[i, t, 4] = sim.getCompConc('comp', 'mol3A')

res_mean = numpy.mean(res, 0)

pylab.plot(tpnt, res_mean[:,0], label = 'A')
pylab.plot(tpnt, res_mean[:,1], label = 'B')
pylab.plot(tpnt, res_mean[:,2], label = 'C')
pylab.plot(tpnt, res_mean[:,3], label = '2A')
pylab.plot(tpnt, res_mean[:,4], label = '3A')


pylab.xlabel('Time (sec)')
pylab.ylabel('# molecules')
pylab.legend()
pylab.show()