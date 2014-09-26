
from bigfilepy import BigFile
import svr
import painter
import numpy
from color import CoolWarm, N, NL

from sys import argv

# This version assumes a small image

M = [
    [2, 2, -3],
    [7, 3, 6],
    [2, 1, 1]]

FULL_SIZE = svr.remap_query_size(M)
print FULL_SIZE
from mpi4py import MPI
import domain
bigfile = BigFile(argv[1])
#bigfile = BigFile('/physics2/yfeng1/BWSim/TEST/TEST-fof4/PART_027')
b = bigfile.open("0/Position")
world = MPI.COMM_WORLD

PIXEL_WIDTH = 2200
PIXEL_HEIGHT = int(FULL_SIZE[0] / FULL_SIZE[1] * PIXEL_WIDTH)
CHUNK_SIZE = 1024 * 128

def readheader(comm):
    if comm.rank == 0:
        attrs = bigfile.open("header").attrs
        HEADER = dict([(i, attrs[i]) for i in attrs])
    else:
        HEADER = None
    return comm.bcast(HEADER)

def ieye2T(ie, ye, Xh=0.76, out=None):
    fac = 121.14740013761634  # GADGET unit Protonmass / Bolztman const
    GAMMA = 5 / 3.
    mu = 1.0 / (ye * Xh + (1 - Xh) * 0.25 + Xh)
    if out != None:
      out[:] = ie[:] * mu * (GAMMA - 1) * fac
      return out
    else:
      return ie * mu * (GAMMA - 1) * fac

def process_chunk(comm, image,
       BoxSize,
       posblock, massblock, smlblock, ieblock, yeblock, 
        start, end):
    T = ieye2T(ieblock[start:end], yeblock[start:end])
    pos = posblock[start:end]

    # to uniform coordinate
    svr.wrap(pos, BoxSize, output=pos)

    # to survey coordinate
    svr.remap(pos, M, output=pos)

    # to device coordinate:
    pos[:, 0] *= PIXEL_HEIGHT / FULL_SIZE[0]
    pos[:, 1] *= PIXEL_WIDTH / FULL_SIZE[1]


    sml = numpy.float32(smlblock[start:end])
    mass = numpy.float32(massblock[start:end])
    T = T * mass
    data = numpy.vstack((mass, T)).T
    sml[:] *= 1.0 / BoxSize * PIXEL_WIDTH / FULL_SIZE[1]

    # to rank device coordinate
    painter.paint(pos, sml, data, image)

    with domain.Rotator(comm):
        print comm.rank, '----'
        print 'rank', 'pos', comm.rank, pos.shape
        try:
            print pos.min(axis=0), pos.max(axis=0)
        except:
            pass
        print 'max 0', image[:, :, 0].max()
        print 'max 1', image[:, :, 1].max()


def main(comm):
    HEADER = readheader(comm)

    NumPartTotal = HEADER['TotNumPart'][0]
    BoxSize = HEADER['BoxSize']
    
    image = numpy.zeros((PIXEL_HEIGHT, PIXEL_WIDTH), dtype=('f4', 2))


    posblock = bigfile.open("0/Position")
    massblock = bigfile.open("0/Mass")
    smlblock = bigfile.open("0/SmoothingLength")
    ieblock = bigfile.open("0/InternalEnergy")
    yeblock = bigfile.open("0/ElectronAbundance")

    FULL_CHUNK_SIZE = comm.size * CHUNK_SIZE
    for i in range(0, NumPartTotal, FULL_CHUNK_SIZE):
        cstart, cend, stop = slice(i, i + FULL_CHUNK_SIZE).indices(NumPartTotal)
        start = comm.rank * (cend - cstart) // comm.size
        end = (comm.rank + 1) * (cend - cstart) // comm.size
        start += cstart 
        end += cstart

        process_chunk(comm, image, 
                BoxSize,
                posblock, massblock, smlblock, ieblock, yeblock, 
                start, end)
        comm.barrier()
        if comm.rank == 0:
            print 'chunk', i, '/', NumPartTotal, cstart, cend
    cpy = image.copy()
    image[...] = 0

    comm.Reduce(cpy, image, MPI.SUM)

    if comm.rank == 0:
        numpy.save('test.npy', image)
        makefig()

def makefig():
    gimage = numpy.load('test.npy')
    meanT = gimage[:, :, 1] / gimage[:, : , 0]
    density = gimage[:, :, 0]
#    print numpy.histogram(NL(density, range=-3), range=(-0.1, 1.2))
    rgb = CoolWarm(NL(meanT), NL(density))
    print rgb.shape, rgb.dtype
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    figure = Figure((PIXEL_WIDTH / 100., PIXEL_HEIGHT / 100.), dpi=100)
    ax = figure.add_axes([0, 0, 1, 1])
    ax.imshow(rgb)
    canvas = FigureCanvasAgg(figure)
    figure.savefig('test.png')

#makefig()
main(MPI.COMM_WORLD)

