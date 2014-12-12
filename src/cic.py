# cloud in cell painting
#
import numpy

def cic(pos, mesh, boxsize=1.0, weights=1.0, mode="raise"):
    """ CIC approximation from points to Nmesh,
        each point has a weight given by weights.
        This does not give density.
        pos is supposed to be row vectors. aka for 3d input
        pos.shape is (?, 3).

        the mesh is an array of size (Nmesh, Nmesh, Nmesh)

        pos[:, i] is mesh.shape[i]
        thus z is the fast moving index

        mode can be :
            "raise" : raise exceptions if a particle is painted
             outside the mesh
            "wrap"  : wrap with periodic boundry
            "ignore": ignore particle contribution outside of the mesh

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
        if mode == 'wrap':
            gridpos = numpy.remainder(pos[chunk], BoxSize) * (Nmesh / BoxSize)
            rmi_mode = 'wrap'
            intpos = numpy.intp(gridpos)
        elif mode == 'raise':
            gridpos = pos[chunk] * (Nmesh / BoxSize)
            rmi_mode = 'raise'
            intpos = numpy.intp(numpy.floor(gridpos))
        elif mode == 'ignore':
            gridpos = pos[chunk] * (Nmesh / BoxSize)
            rmi_mode = 'raise'
            intpos = numpy.intp(numpy.floor(gridpos))

        for i, neighbour in enumerate(neighbours):
            neighbour = neighbour[None, :]
            targetpos = intpos + neighbour

            kernel = (1.0 - numpy.abs(gridpos - targetpos)).prod(axis=-1)
            add = wchunk * kernel

            if mode == 'ignore':
                # filter out those outside of the mesh
                mask = (targetpos >= 0).all(axis=-1)
                for d in range(Ndim):
                    mask &= (targetpos[..., d] < mesh.shape[d])
                targetpos = targetpos[mask]
                add = add[mask]

            if len(targetpos) > 0:
                targetindex = numpy.ravel_multi_index(
                        targetpos.T, mesh.shape, mode=rmi_mode)
                u, label = numpy.unique(targetindex, return_inverse=True)
                mesh.flat[u] += numpy.bincount(label, add, minlength=len(u))

    return mesh

