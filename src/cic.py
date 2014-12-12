# cloud in cell painting
#
import numpy

#mesh = numpy.zeros(shape=(Nmesh, ) * Ndim, 
#        dtype=dtype, order='C')
def cic(pos, mesh, boxsize=1.0, weights=1.0, periodic=False):
    """ CIC approximation from points to Nmesh,
        each point has a weight given by weights.
        This does not give density.
        pos is supposed to be row vectors. aka for 3d input
        pos.shape is (?, 3).

        the mesh is an array of size (Nmesh, Nmesh, Nmesh)

        pos[:, i] is mesh.shape[i]
        thus z is the fast moving index
    """
    pos = numpy.array(pos)
    chunksize = 1024 * 16 * 4
    BoxSize = 1.0 * boxsize
    Ndim = pos.shape[-1]
    Np = pos.shape[0]
    Nmesh = mesh.shape[0]

    neighbours = ((numpy.arange(2 ** Ndim)[:, None] >> \
            numpy.arange(Ndim)[None, :]) & 1)
    for start in range(0, Np, chunksize):
        chunk = slice(start, start+chunksize)
        if numpy.isscalar(weights):
          wchunk = weights
        else:
          wchunk = weights[chunk]
        if periodic:
            gridpos = numpy.remainder(pos[chunk], BoxSize) * (Nmesh / BoxSize)
            mode = 'wrap'
            intpos = numpy.intp(gridpos)
        else:
            gridpos = pos[chunk] * (Nmesh / BoxSize)
            mode = 'raise'
            intpos = numpy.intp(numpy.floor(gridpos))
        for i, neighbour in enumerate(neighbours):
            neighbour = neighbour[None, :]
            targetpos = intpos + neighbour

            targetindex = numpy.ravel_multi_index(
                    targetpos.T, mesh.shape, mode=mode)

            kernel = (1.0 - numpy.abs(gridpos - targetpos)).prod(axis=-1)
            add = wchunk * kernel
            u, label = numpy.unique(targetindex, return_inverse=True)
            mesh.flat[u] += numpy.bincount(label, add, minlength=len(u))
    return mesh

