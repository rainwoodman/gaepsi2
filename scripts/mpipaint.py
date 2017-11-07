import sys
import os.path

from pypm import domain

d = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.append(d)

from bigfile import BigFile

import svr
import painter
import camera

import numpy
import os
from color import CoolWarm, N, NL

from sys import argv
from mpi4py import MPI

class UniqueBarrier(object):
    def __init__(self, comm):
        self.comm = comm

    def check(self, tag):
        self.comm.Barrier()
        newtag = None
        if self.comm.rank == 0:
            newtag = tag
        newtag = self.comm.bcast(newtag)
        if newtag != tag:
            raise Exception("rank = %d, tag is %s, on root tag is %s" %
                    (self.comm.rank, tag, newtag))
        self.comm.Barrier()
        if self.comm.rank == 0:
            print(tag)

config = argv[1]
input = argv[2]
output = argv[3]

bigfile = BigFile(input)

world = MPI.COMM_WORLD

if world.rank == 0:
    attrs = bigfile.open("Header").attrs
    HEADER = dict([(i, attrs[i]) for i in attrs])
else:
    HEADER = None
HEADER = world.bcast(HEADER)
BoxSize = HEADER['BoxSize']

# set up some defaults
SMLFACTOR = 1.0
TilePadding = 256

M = [
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1]]

def ieye2Tmass(ie, ye, mass, Xh=0.76, out=None):
    fac = 121.14740013761634  # GADGET unit Protonmass / Bolztman const
    GAMMA = 5 / 3.
    mu = 1.0 / (ye * Xh + (1 - Xh) * 0.25 + Xh)
    if out != None:
      out[:] = ie[:] * mu * (GAMMA - 1) * fac * mass
      return out
    else:
      return ie * mu * (GAMMA - 1) * fac * mass

PIXEL_WIDTH = 2200

CHUNK_SIZE = 1024 * 128 
SMALL_IMAGE = PIXEL_WIDTH <= 10000
POS_BLOCK = "0/Position"
SML_BLOCK = "0/SmoothingLength"
DATA_BLOCK =  ["0/Mass", 
        "0/StarFormationRate", 
        (ieye2Tmass, "0/InternalEnergy", "0/ElectronAbundance", "0/Mass")]
CENTER = [0.5, 0.5]
VIEW_SIZE = [1.0, 1.0]
# pull in M
execfile(config)

if M is not None:
    FULL_SIZE = svr.remap_query_size(M)
else:
    FULL_SIZE = numpy.ones(3)
FULL_SIZE *= BoxSize

NEAR = 0
FAR = FULL_SIZE[2]
EXTENT = (-FULL_SIZE[0] * 0.5 * VIEW_SIZE[0], FULL_SIZE[0] * 0.5 * VIEW_SIZE[0],
        -FULL_SIZE[1] * 0.5 * VIEW_SIZE[1], FULL_SIZE[1] * 0.5 * VIEW_SIZE[1])
CAMERA = (FULL_SIZE[0] * CENTER[0], FULL_SIZE[1] * CENTER[1], 0)
FOCUS = (FULL_SIZE[0] * CENTER[0], FULL_SIZE[1] * CENTER[1], FULL_SIZE[2])
UP = (0, 1, 0)

PIXEL_HEIGHT = int(FULL_SIZE[0] * VIEW_SIZE[0] / (FULL_SIZE[1] *VIEW_SIZE[1]) * PIXEL_WIDTH)

# do it again to make sure we do not override other confs
execfile(config)
def buildgrid(PX, NX):
    ntile = PX // TilePadding + (PX % TilePadding != 0)
    gridx = numpy.arange(NX + 1) * ntile // NX * TilePadding
    return gridx
def process_chunk(image,
       BoxSize,
       pos, sml, data, d2d=None):
    # apply a data mask
    mask = (data!= 0).any(axis=-1)

    pos = pos[mask]
    sml = sml[mask]
    data = data[mask]

    # to uniform coordinate [0-1]
    svr.wrap(pos, BoxSize, output=pos)

    # to survey coordinate
    if M is not None:
        svr.remap(pos, M, output=pos)

    # to simulation coordinate 
    pos *= BoxSize

    # simulate a camera:

    model = camera.lookat(pos=CAMERA, target=FOCUS, up=UP)
    pers = camera.ortho(near=NEAR, far=FAR, extent=EXTENT)
    matrix = numpy.dot(pers, model)

    pos = camera.apply(matrix, pos)

    # apply near field far field cut
    mask = pos[:, 2] >= -1
    if mask is not None:
        mask &= pos[:, 2] <= 1.0
    else:
        mask = pos[:, 2] <= 1.0
    if mask is not None:
        pos = pos[mask]
        sml = sml[mask]
        data = data[mask]

    pos += 1.0
    if d2d is None or d2d.comm.rank == 0:
        print('boxsize', BoxSize)
    # to device coordinate:
    pos[:, 0] *= PIXEL_HEIGHT / 2.
    pos[:, 1] *= PIXEL_WIDTH / 2.

    sml[:] *= PIXEL_WIDTH / (FULL_SIZE[1] * VIEW_SIZE[1])
    assert len(sml) == len(pos)
    if d2d is not None:

        layout = d2d.decompose(pos, smoothing=sml)
        pos = layout.exchange(pos)
        sml = layout.exchange(sml)
        data = layout.exchange(data)
        
        # to rank device coordinate
        pos[:, 0] -= d2d.mystart[0]
        pos[:, 1] -= d2d.mystart[1]

    sml[sml > 100] = 100
    if d2d is None or d2d.comm.rank == 0:
        print(sml.max())
    image[...] += painter.paint(pos, sml, data, image.shape[1:], np=0)

    if d2d is not None:
        # get some stats
        try:
            posstats = (pos.min(axis=0), pos.max(axis=0), d2d.mystart, d2d.myend)
        except:
            posstats = None
        try:
            datastats = ('data sum', data.sum(axis=0))
        except:
            datastats = None
        try:
            imgstats = ('imgsum', image.sum(axis=(0, 1)))
        except:
            imgstats = None
    
        comm = d2d.comm
        with domain.Rotator(comm):
            print(comm.rank, '----')
            print('rank', 'pos', comm.rank, pos.shape, image.shape, data.shape)
            print(posstats)
            print(datastats)
            print(imgstats)
    else:
        print('d2d is None')

def main(comm):
    posblock = bigfile.open(POS_BLOCK)
    smlblock = bigfile.open(SML_BLOCK)
    
    NumPartTotal = smlblock.size
    
    NX = int(comm.size ** 0.5)
    if NX < 1: NX = 1
    while comm.size % NX : 
        NX = NX - 1
    NY = comm.size // NX
     
    if not SMALL_IMAGE:
        d2d = domain.GridND(
                grid=[buildgrid(PIXEL_HEIGHT, NX),
                      buildgrid(PIXEL_WIDTH, NY)],
                comm=comm,
                periodic=True)

        myoffset = (d2d.mystart[0], d2d.mystart[1])
        mysize = (d2d.myend[0] - d2d.mystart[0], 
                d2d.myend[1] - d2d.mystart[1])
    else:
        d2d = None
        mysize = (PIXEL_HEIGHT, PIXEL_WIDTH)
        
    image = numpy.zeros((len(DATA_BLOCK),) + mysize, dtype=('f4'))

    datablock = [
        tuple([item[0]] 
      + [bigfile.open(b) for b in item[1:]])
        if isinstance(item, tuple) else \
        bigfile.open(item) 
      for item in DATA_BLOCK ]

    FULL_CHUNK_SIZE = comm.size * CHUNK_SIZE
    for i in range(0, NumPartTotal, FULL_CHUNK_SIZE):
        cstart, cend, stop = slice(i, i + FULL_CHUNK_SIZE).indices(NumPartTotal)
        start = comm.rank * (cend - cstart) // comm.size
        end = (comm.rank + 1) * (cend - cstart) // comm.size
        start += cstart 
        end += cstart

        pos = posblock[start:end]
        sml = numpy.float32(smlblock[start:end]) * SMLFACTOR

        data = [
            item[0] (*[numpy.float32(b[start:end]) for b in item[1:]])
            if isinstance(item, tuple) else \
            numpy.float32(item[start:end])
          for item in datablock ]
        data = numpy.vstack(data).T

        process_chunk(image, BoxSize,
                pos, sml, data, d2d=d2d)

        comm.Barrier()
        if comm.rank == 0:
            print('chunk', i, '/', NumPartTotal, cstart, cend)
        
    try:
        os.makedirs(output)
    except OSError:
        pass
    if not SMALL_IMAGE:
        if image.size > 0:
            numpy.save(os.path.join(output, 'imagetile-%05d-%05d.npy' % myoffset), image)
    else:
        gimage = image.copy()
        image[...] = 0
        comm.Reduce(gimage, image, MPI.SUM)
        if comm.rank == 0:
            numpy.save(os.path.join(output, 'imagetile.npy'), image)


main(world)
