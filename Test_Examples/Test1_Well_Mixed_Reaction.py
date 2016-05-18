
# Well Mixed Reaction A <=> B

import steps.model as smodel
import steps.geom as swm
import steps.rng as srng
import steps.solver as ssolver
import numpy
import pylab


mdl = smodel.Model()

molA = smodel.Spec('molA', mdl)
molB = smodel.Spec('molB', mdl)
molC = smodel.Spec('molC', mdl)

volsys = smodel.Volsys('vsys',mdl)

kreac_f = smodel.Reac('kreac_f', volsys, lhs = [molA, molB], rhs = [molC], kcst = 0.3e6)
kreac_b = smodel.Reac('kreac_b', volsys, lhs = [molC], rhs = [molA, molB], kcst = 0.7)

#kreac_f.setKcst = 0.3e6
#kreac_b.setKcst = 0.7

wmgeom = swm.Geom()

comp = swm.Comp('comp', wmgeom)
comp.addVolsys('vsys')
comp.setVol(1.6667e-21)

r = srng.create('mt19937', 256)
r.initialize(23412)



sim = ssolver.Wmdirect(mdl, wmgeom, r)
sim.reset()
sim.setCompConc('comp','molA', 31.4e-6)
sim.setCompConc('comp','molB', 22.3e-6)

tpnt = numpy.arange(0, 2.001, 0.001)
res =  numpy.zeros([2001, 3])

for t in range(0, 2001):
    sim.run(tpnt[t])
    res[t, 0] = sim.getCompCount('comp', 'molA')
    res[t, 1] = sim.getCompCount('comp', 'molB')
    res[t, 2] = sim.getCompCount('comp', 'molC')

pylab.plot(tpnt, res[:,0], label = 'A')
pylab.plot(tpnt, res[:,1], label = 'B')
pylab.plot(tpnt, res[:,2], label = 'C')

pylab.xlabel('Time (sec)')
pylab.ylabel('# molecules')
pylab.legend()
pylab.show()

NITER = 100
res = numpy.zeros([NITER, 2001, 3])
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

res_mean = numpy.mean(res, 0)

pylab.plot(tpnt, res_mean[:,0], label = 'A')
pylab.plot(tpnt, res_mean[:,1], label = 'B')
pylab.plot(tpnt, res_mean[:,2], label = 'C')

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
        res[i, t, 2] = sim.getCompConc('comp', 'molC')

    sim.setCompCount('comp', 'molA', sim.getCompCount('comp', 'molA') + 10)

    for t in range(1001, 2001):
        sim.run(tpnt[t])
        res[i, t, 0] = sim.getCompConc('comp', 'molA')
        res[i, t, 1] = sim.getCompConc('comp', 'molB')
        res[i, t, 2] = sim.getCompConc('comp', 'molC')

res_mean = numpy.mean(res, 0)

pylab.plot(tpnt, res_mean[:,0], label = 'A')
pylab.plot(tpnt, res_mean[:,1], label = 'B')
pylab.plot(tpnt, res_mean[:,2], label = 'C')

pylab.xlabel('Time (sec)')
pylab.ylabel('# molecules')
pylab.legend()
pylab.show()

# import matplotlib.pyplot as plt
# import matplotlib.font_manager as font_manager
#
# # Set the font dictionaries (for plot title and axis titles)
# title_font = {'fontname':'Arial', 'size':'16', 'color':'black', 'weight':'normal',
#               'verticalalignment':'bottom'} # Bottom vertical alignment for more space
# axis_font = {'fontname':'Arial', 'size':'14'}
#
# # Set the font properties (for use in legend)
# font_path = 'C:\Windows\Fonts\Arial.ttf'
# font_prop = font_manager.FontProperties(fname=font_path, size=14)
#
# ax = plt.subplot() # Defines ax variable by creating an empty plot
#
# # Set the tick labels font
# for label in (ax.get_xticklabels() + ax.get_yticklabels()):
#     label.set_fontname('Arial')
#     label.set_fontsize(13)
#
# x = numpy.linspace(0, 10)
# y = x + numpy.random.normal(x) # Just simulates some data
#
# plt.plot(x, y, 'b+', label='Data points')
# plt.xlabel("x axis", **axis_font)
# plt.ylabel("y axis", **axis_font)
# plt.title("Misc graph", **title_font)
# plt.legend(loc='lower right', prop=font_prop, numpoints=1)
# plt.text(0, 0, "Misc text", **title_font)
# plt.show()