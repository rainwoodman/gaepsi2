from mpi4py import MPI
cimport numpy
import numpy

def sqrtint(x):
    y = int(x ** 0.5) - 1
    while y > 1:
        if x % y == 0: break
        y = y - 1
    return y

cdef numpy.ndarray _digitize(data, bins, period):
    if period:
        data = numpy.remainder(data, bins[-1], output=data)
    if len(data) == 0:
        return numpy.empty((0), dtype='intp')
    else:
        return numpy.digitize(data, bins)

cdef numpy.ndarray A(b):
    return numpy.asarray(b)

class Rotator(object):
    def __init__(self, comm):
        self.comm = comm
    def __enter__(self):
        for i in range(self.comm.rank):
            self.comm.barrier()
    def __exit__(self, type, value, tb):
        for i in range(self.comm.rank, self.comm.size):
            self.comm.barrier()

class Layout(object):
    def __init__(self, comm, sendcounts, indices):
        self.comm = comm
        assert self.comm.size == sendcounts.shape[0]

        self.sendcounts = sendcounts
        self.recvcounts = numpy.empty_like(self.sendcounts)

        self.sendoffsets = numpy.zeros_like(self.sendcounts)
        self.recvoffsets = numpy.zeros_like(self.recvcounts)

        self.comm.Alltoall(self.sendcounts, self.recvcounts)

        self.sendoffsets[1:] = self.sendcounts.cumsum()[:-1]
        self.recvoffsets[1:] = self.recvcounts.cumsum()[:-1]

        self.oldlength = self.sendcounts.sum()
        self.newlength = self.recvcounts.sum()

        self.indices = indices
    def exchange(self, data):
        #build buffer
        buffer = data.take(self.indices, axis=0)
        sendbuffers = numpy.split(buffer, self.sendoffsets[1:])

        newshape = list(data.shape)
        newshape[0] = self.newlength

        duplicity = numpy.product(data.shape[1:]) 

        itemsize = duplicity * data.dtype.itemsize

        dt = MPI.BYTE.Create_contiguous(itemsize)

        dt.Commit()
        recvbuffer = numpy.empty(newshape, dtype=buffer.dtype)
        recvbuffers = numpy.split(recvbuffer, self.recvcounts)
        
        self.comm.Alltoallv((buffer, (self.sendcounts, self.sendoffsets), dt), 
                            (recvbuffer, (self.recvcounts, self.recvoffsets), dt))
        dt.Free()
        return recvbuffer


class Grid2D(object):
    """
        2D domain decomposition on a uniform grid
    """
    def __init__(self, 
            gridx,
            gridy,
            comm=MPI.COMM_WORLD,
            periodic=True):
        """ gridy is the fast changing dimension (aka this is not the image
            coordinate system
        """
        self.dims = (len(gridx) - 1, len(gridy) - 1)
        self.gridx = gridx
        self.gridy = gridy
        self.periodic = periodic
        self.comm = comm
        #self.comm2D = comm.Create_cart(self.dims, periods=[periodic for i in self.dims])

        cdef int ny = self.dims[1]
        rank = [ self.comm.rank // ny, self.comm.rank % ny]
        self.start = numpy.array([self.gridx[rank[0]], self.gridy[rank[1]]])
        self.end = numpy.array([self.gridx[rank[0] + 1], self.gridy[rank[1] + 1]])


    def decompose(self, pos, bleeding=0):
        """ decompose the domain according to pos,

            returns a Layout object that can be used
            to exchange data
        """

        x, y = numpy.asarray(pos).T

        target = numpy.empty((len(x), 4), 'int16')
        cdef short int [:] xl
        cdef short int [:] yl
        cdef short int [:] xr
        cdef short int [:] yr

        cdef numpy.intp_t i

        cdef numpy.intp_t u, v
        cdef int nx = self.dims[0]
        cdef int ny = self.dims[1]

        cdef numpy.intp_t [:, ::1] counts2d
        cdef numpy.intp_t [:, ::1] offsets2d 
        cdef numpy.intp_t [::1] counts

        cdef numpy.intp_t [::1] indices
        cdef numpy.intp_t [:, ::1] ptrs
        cdef int periodic = self.periodic

        xl, xr, yl, yr = target.T

        A(xl)[...] = _digitize(x - bleeding, self.gridx, self.periodic) - 1
        A(xr)[...] = _digitize(x + bleeding, self.gridx, self.periodic)
        A(yl)[...] = _digitize(y - bleeding, self.gridy, self.periodic) - 1
        A(yr)[...] = _digitize(y + bleeding, self.gridy, self.periodic)

        countsobj = numpy.zeros(self.dims, 'intp')
        counts2d = countsobj
        counts = countsobj.reshape(-1)

        if periodic:
            for i in range(xl.shape[0]):
                if xr[i] < xl[i]:
                    xr[i] += nx
                if yr[i] < yl[i]:
                    yr[i] += ny
        else:
            A(xl).clip(0, nx, out=A(xl))
            A(yl).clip(0, ny, out=A(yl))
            A(xr).clip(0, nx, out=A(xr))
            A(yr).clip(0, ny, out=A(yr))

        for i in range(xl.shape[0]):
            for u in range(xl[i], xr[i]):
                for v in range(yl[i], yr[i]):
                    u = u % nx
                    v = v % ny
                    counts2d[u, v] = counts2d[u, v] + 1

        cumsum = countsobj.cumsum()

        total = cumsum[-1]
        offsetsobj = numpy.zeros_like(cumsum)
        offsetsobj[1:] = cumsum[:-1]

        offsets2d = offsetsobj.reshape(self.dims)
        ptrs = offsets2d.copy()

        indices = numpy.empty(total, dtype='intp')
        for i in range(xl.shape[0]):
            for u in range(xl[i], xr[i]):
                for v in range(yl[i], yr[i]):
                    u = u % nx
                    v = v % ny
                    indices[ptrs[u, v]] = i
                    ptrs[u, v] = ptrs[u, v] + 1
        layout = Layout(
                comm=self.comm,
                sendcounts=A(counts2d).reshape(-1),
                indices=A(indices))

        return layout
