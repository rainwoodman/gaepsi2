
from bigfilepy import BigFile
import svr
import domain
import painter
import numpy
from color import CoolWarm, N, NL

from sys import argv
from mpi4py import MPI

#bigfile = BigFile(argv[1])
bigfile = BigFile('/physics2/yfeng1/BWSim/TEST/TEST-fof4/PART_027')
b = bigfile.open("0/Position")
world = MPI.COMM_WORLD

TilePadding = 256

M = [
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1]]

FULL_SIZE = svr.remap_query_size(M)

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
def buildgrid(PX, NX):
    ntile = PX // TilePadding + (PX % TilePadding != 0)
    gridx = numpy.arange(NX + 1) * ntile // NX * TilePadding
    return gridx

def ieye2T(ie, ye, Xh=0.76, out=None):
    fac = 121.14740013761634  # GADGET unit Protonmass / Bolztman const
    GAMMA = 5 / 3.
    mu = 1.0 / (ye * Xh + (1 - Xh) * 0.25 + Xh)
    if out != None:
      out[:] = ie[:] * mu * (GAMMA - 1) * fac
      return out
    else:
      return ie * mu * (GAMMA - 1) * fac

def process_chunk(comm, d2d, image,
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


    layout = d2d.decompose(pos[:, :2])

    pos = layout.exchange(pos)
    sml = layout.exchange(numpy.float32(smlblock[start:end]))
    mass = layout.exchange(numpy.float32(massblock[start:end]))
    T = layout.exchange(T)
    T = T * mass
    data = numpy.vstack((mass, T)).T
    sml[:] *= 1.0 / BoxSize * PIXEL_WIDTH / FULL_SIZE[1]

    if comm.rank == 0:
        print d2d.gridx
        print d2d.gridy

    # to rank device coordinate
    pos[:, 0] -= d2d.start[0]
    pos[:, 1] -= d2d.start[1]

    painter.paint(pos, sml, data, image)

    with domain.Rotator(comm):
        print comm.rank, '----'
        print 'px start', d2d.start
        print 'px end', d2d.end
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
    
    NX = int(comm.size ** 0.5)
    while comm.size % NX : 
        NX = NX - 1
    NY = comm.size // NX
    
    d2d = domain.Grid2D(
            gridx=buildgrid(PIXEL_HEIGHT, NX),
            gridy=buildgrid(PIXEL_WIDTH, NY),
            comm=comm,
            bleeding=30,
            periodic=False)

    mysize = (d2d.end[0] - d2d.start[0], 
            d2d.end[1] - d2d.start[1])
    image = numpy.zeros(mysize, dtype=('f4', 2))


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

        process_chunk(comm, d2d, image, 
                BoxSize,
                posblock, massblock, smlblock, ieblock, yeblock, 
                start, end)
        comm.barrier()
        if comm.rank == 0:
            print 'chunk', i, '/', NumPartTotal, cstart, cend

    tiles = image.reshape(
            image.shape[0] // TilePadding,
            TilePadding,
            image.shape[1] // TilePadding,
            TilePadding,
            image.shape[2]).transpose((1,3,0,2,4))

    images = comm.gather(image)
    starts = comm.gather(d2d.start)
    ends = comm.gather(d2d.end)
    if comm.rank == 0:
        gimage = numpy.zeros(
                (
                    d2d.gridx[-1],
                    d2d.gridy[-1], 
                    image.shape[2]))
        for start, end, image in zip(starts, ends, images):
            gimage[start[0]:end[0],
                    start[1]:end[1],
                    :] = image

        gimage[PIXEL_HEIGHT:, :, :] = 0
        gimage[:, PIXEL_WIDTH:, :] = 0
        numpy.save('test.npy', gimage)
        makefig()

def makefig():
    gimage = numpy.load('test.npy')
    meanT = gimage[:, :, 1] / gimage[:, : , 0]
    density = gimage[:, :, 0]
    print numpy.histogram(NL(density, range=-3), range=(-0.1, 1.2))
    rgb = CoolWarm(NL(meanT), NL(density, range=-3))
    print rgb.shape, rgb.dtype
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    figure = Figure((PIXEL_WIDTH / 100., PIXEL_HEIGHT / 100.), dpi=100)
    ax = figure.add_axes([0, 0, 1, 1])
    ax.imshow(rgb)
    canvas = FigureCanvasAgg(figure)
    figure.savefig('test.png')

def tilename(i, j, imax, jmax):
    """ i is y, j is x
        randy: 
             0 1
             2 3 """
    t = 1
    r = []
    while t < imax or t < jmax:
        n = 0
        if (i & t) != 0:
            n += 2
        if (j & t) != 0:
            n += 1
        t = t * 2
        r.append('0%s' % n)
    return '/'.join(r[::-1])

class TileStore(object):
    def __init__(self, dtype, nx, ny):
        self.dtype = dtype

    def save(self, image, tilex, tiley):
        # build the linear index
        # cut out the first piece for filename
        # cut out the first piece for filename
        filename, remainding = cutoff(index)
        with file(filename) as f:
            f.seek(index * imagesize)
            image.tofile(f)
makefig()
#main(MPI.COMM_WORLD)

