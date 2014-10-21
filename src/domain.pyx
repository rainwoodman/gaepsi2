"""
   Domain Decomposition in Gaepsi

    currently we have a GridND decomposition algorithm.

"""
from mpi4py import MPI
cimport numpy
import numpy
import cython
cimport cython

def sqrtint(x):
    y = int(x ** 0.5) - 1
    while y > 1:
        if x % y == 0: break
        y = y - 1
    return y

cdef numpy.ndarray _digitize(data, bins):
    if len(data) == 0:
        return numpy.empty((0), dtype='intp')
    else:
        return numpy.digitize(data, bins)

cdef numpy.ndarray A(b):
    """ we use A to convert a memoryslice /view to array"""
    return numpy.asarray(b)

class Rotator(object):
    def __init__(self, comm):
        self.comm = comm
    def __enter__(self):
        self.comm.Barrier()
        for i in range(self.comm.rank):
            self.comm.Barrier()
    def __exit__(self, type, value, tb):
        for i in range(self.comm.rank, self.comm.size):
            self.comm.Barrier()
        self.comm.Barrier()

cdef class Layout(object):
    """ A global all to all communication layout 
        
    """
    cdef readonly numpy.ndarray sendcounts
    cdef readonly numpy.ndarray sendoffsets 
    cdef readonly numpy.ndarray recvoffsets 
    cdef readonly numpy.ndarray recvcounts
    cdef readonly object comm
    cdef readonly numpy.ndarray indices
    cdef readonly object newlength
    cdef readonly object oldlength
    def __init__(self, comm, sendcounts, indices, recvcounts=None):
        """
        sendcounts is the number of items to send
        indices is the indices of the items in the data array.
        """

        self.comm = comm
        assert self.comm.size == sendcounts.shape[0]

        self.sendcounts = numpy.array(sendcounts, order='C')
        self.recvcounts = numpy.empty_like(self.sendcounts, order='C')

        self.sendoffsets = numpy.zeros_like(self.sendcounts, order='C')
        self.recvoffsets = numpy.zeros_like(self.recvcounts, order='C')

        if recvcounts is None:
            # calculate the recv counts array
            # ! Alltoall
            self.comm.Barrier()
            self.comm.Alltoall(self.sendcounts, self.recvcounts)
            self.comm.Barrier()
        else:
            self.recvcounts = recvcounts
        self.sendoffsets[1:] = self.sendcounts.cumsum()[:-1]
        self.recvoffsets[1:] = self.recvcounts.cumsum()[:-1]

        self.oldlength = self.sendcounts.sum()
        self.newlength = self.recvcounts.sum()

        self.indices = indices

    def exchange(self, data):
        """ exchange the data globally according to the layout
            data shall be of the same length of the input position
            that builds the layout
        """
        #build buffer
        # Watch out: 
        # take produces C-contiguous array, 
        # friendly to alltoallv.
        # fancy indexing does not always return C_contiguous
        # array (2 days to realize this!)
        
        cdef numpy.ndarray buffer = data.take(self.indices, axis=0)

        newshape = list(data.shape)
        newshape[0] = self.newlength

        # build a dtype for communication
        # this is to avoid 2GB limit from bytes.
        duplicity = numpy.product(numpy.array(data.shape[1:], 'intp')) 
        itemsize = duplicity * data.dtype.itemsize
        dt = MPI.BYTE.Create_contiguous(itemsize)
        dt.Commit()

        cdef numpy.ndarray recvbuffer = numpy.empty(newshape, dtype=data.dtype, order='C')
        self.comm.Barrier()

        # now fire
        rt = self.comm.Alltoallv((buffer, (self.sendcounts, self.sendoffsets), dt), 
                            (recvbuffer, (self.recvcounts, self.recvoffsets), dt))
        dt.Free()
        self.comm.Barrier()
        return recvbuffer


class GridND(object):
    """
        ND domain decomposition on a uniform grid
    """
    def __init__(self, 
            grid,
            comm=MPI.COMM_WORLD,
            periodic=True):
        """ 
            grid is a list of  grid edges. 
            grid[0] or pos[:, 0], etc.
        
            grid[i][-1] are the boxsizes
            the ranks are set up into a mesh of 
                len(grid[0]) - 1, ...
        """
        self.dims = [len(g) - 1 for g in grid]
        self.grid = numpy.asarray(grid)
        self.periodic = periodic
        self.comm = comm
        assert comm.size == numpy.product(self.dims)
        rank = numpy.unravel_index(self.comm.rank, self.dims)

        self.myrank = numpy.array(rank)
        self.mystart = numpy.array([g[r] for g, r in zip(grid, rank)])
        self.myend = numpy.array([g[r + 1] for g, r in zip(grid, rank)])

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.overflowcheck(False)
    @cython.nonecheck(False)
    def _fill(self, int mode,
            int Ndim,
            int Npoint,
            short int [:, ::1] sil,
            short int [:, ::1] sir,
            int periodic,
            int [::1] counts
            ):
        # first count the sizes
        cdef int patchsize = 0
        cdef int p[32]
        cdef int strides[32]
        cdef int dims[32]
        cdef numpy.int32_t [::1] offset = None
        cdef int i
        cdef int j
        cdef int jj
        cdef int k
        cdef int t
        cdef int [::1] indices = None
        cdef int Nrank
        if mode == 1:
            Nrank = self.comm.size
            offset = numpy.empty(Nrank + 1, dtype='int32', order='C')
            offset[0] = 0
            for i in range(1, Nrank + 1):
                offset[i] = offset[i - 1] + counts[i - 1]
            totalsize = offset[Nrank]
            indices = numpy.empty(totalsize, dtype='int32', order='C')

        for j in range(Ndim):
            dims[j] = self.dims[j]

        strides[Ndim -1 ] = 1
        for j in range(Ndim - 2, -1, -1):
            strides[j] = strides[j + 1] * dims[j + 1]
        cdef numpy.intp_t target
        for i in range(Npoint):
            patchsize = 1
            for j in range(Ndim):
                patchsize *= sir[j, i] - sil[j, i]
                p[j] = sil[j, i]
            for k in range(patchsize):
                target = 0
                for j in range(Ndim):
                    t = p[j]
                    if periodic:
                        while t >= dims[j]:
                            t -= dims[j]
                        while t < 0:
                            t += dims[j]
                    target = target + t * strides[j]
                if mode == 0:
                    counts[target] += 1
                else:
                    indices[offset[target]] = i
                    offset[target] += 1
                p[Ndim - 1] += 1
                # advance
                for jj in range(Ndim-1, 0, -1):
                    if p[jj] == sir[jj, i]:
                        p[jj] = sil[jj, i]
                        p[jj - 1] += 1
                    else:
                        break
        return indices

    def decompose(self, pos, bleeding=0):
        """ decompose the domain according to pos,

            returns a Layout object that can be used
            to exchange data
        """

        # we can't deal with too many points per rank, by  MPI
        assert len(pos) < 1024 * 1024 * 1024 * 2
        posT = numpy.asarray(pos).T
        cdef short int [:, ::1] sil
        cdef short int [:, ::1] sir
        cdef int periodic = self.periodic
        cdef int Npoint
        cdef int Ndim
        cdef short int [:, ::1] width
        cdef numpy.intp_t i
        cdef int j
        cdef int [::1] counts
        cdef int [::1] indices

        Npoint = len(pos)
        Ndim = len(self.dims)
        counts = numpy.zeros(self.comm.size, dtype='int32')

        if Npoint != 0:
            sil = numpy.empty((Ndim, Npoint), dtype='i2', order='C')
            sir = numpy.empty((Ndim, Npoint), dtype='i2', order='C')
            for j in range(Ndim):
                dim = self.dims[j]
                if periodic:
                    tmp = numpy.remainder(posT[j], self.grid[j][-1])
                else:
                    tmp = posT[j]
                A(sil)[j, :] = _digitize(posT[j] - bleeding, self.grid[j]) - 1
                A(sir)[j, :] = _digitize(posT[j] + bleeding, self.grid[j])
                if periodic:
                    for i in range(Npoint):
                        if sir[j, i] < sil[j, i]:
                            sir[j, i] = sir[j, i] + dim
                else:
                    numpy.clip(sil[j], 0, dim, out=A(sil[j]))
                    numpy.clip(sir[j], 0, dim, out=A(sir[j]))
        else:
            sil = numpy.empty((Ndim, 1), dtype='i2', order='C')
            sir = numpy.empty((Ndim, 1), dtype='i2', order='C')

        self._fill(0, Ndim, Npoint, sil, sir, periodic, counts)

        # now lets build the indices array.
        indices = self._fill(1, Ndim, Npoint, sil, sir, periodic, counts)

        # create the layout object
        layout = Layout(
                comm=self.comm,
                sendcounts=A(counts),
                indices=A(indices))

        return layout

