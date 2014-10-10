import svr
M=[[1, 1, 0], [0, 1, 0], [0, 0, 1]]; print svr.remap([[0, 0, 1.]], M)

import domain
from mpi4py import MPI
import numpy

d2d = domain.Grid2D(
        gridx=numpy.linspace(0, 100., 3, endpoint=True),
        gridy=numpy.linspace(0, 100., 3, endpoint=True),
        periodic=True,
        bleeding=1)

print d2d.gridx
print d2d.gridy
pos = numpy.empty((10,2))
pos[:, 0] = numpy.linspace(0, 100, len(pos))
pos[:, 1] = 40

layout = d2d.decompose(pos)

pos2 = layout.exchange(pos)
with domain.Rotator(MPI.COMM_WORLD):
    print MPI.COMM_WORLD.rank, pos2


import painter
values = numpy.ones((len(pos), 1))
sml = numpy.ones(len(pos)) * 10
image = numpy.zeros((100, 100, 1))
painter.paint(pos, sml, values, image)
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
fig = Figure()
ax = fig.add_subplot(111)
ax.imshow(image[..., 0])
canvas = FigureCanvasAgg(fig)
fig.savefig('test.png')
