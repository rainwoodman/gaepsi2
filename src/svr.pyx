import numpy
cimport numpy

cdef extern from 'c/transform.c':
    pass
cdef extern from 'c/svremap.c':
    ctypedef struct SVRemap "SVRemap":
        double size[3]
        pass
    int svremap_init(SVRemap  * r, int * remap)
    double svremap_apply(SVRemap * r, double x[3], double y[3], int I[3])

def wrap(pos, BoxSize, offset=None, output=None):
    """ wrapping pos in a box of BoxSize into (0, 1)
        assuming the pos is of shape [..., Nd]
    """
    norm = numpy.empty(pos.shape[-1])
    norm[:] = BoxSize
    norm[:] = 1.0 / norm[:]
    off = numpy.empty(pos.shape[-1])
    off[:] = offset
    if output is None:
        output = numpy.empty_like(pos)
    if offset and offset != 0.0 :
        output = numpy.substract(pos, offset, output=output)
    numpy.multiply(pos, norm, output=output)
    return numpy.remainder(output, 1.0, output=output)

def remap_query_size(M):
    cdef SVRemap r 
    cdef int[:, ::1] matrix
    cdef double[::1] size
    matrix = numpy.array(M, dtype='int32')
    svremap_init(&r, &matrix[0, 0])
    return numpy.array(<double [:3]> r.size).copy()

def remap(pos, M, output=None):
    """ the output if provided shall be of shape [, 3]
        pos shall be of shape [, 3], and between [0, 1] in all dimensions.
        M is a 3x3 integer matrix
        currently no sanity checks are performed.
    """
    if output is None: 
        output = numpy.empty_like(pos, dtype='f8')
    else:
        assert output.dtype == numpy.dtype('f8')

    assert output.shape[-1] == 3

    output[...] = pos

    cdef int[:, ::1] matrix
    cdef int I[3] 
    cdef double[:, ::1] x = output
    cdef double[:, ::1] y = output
    cdef double tmp[3]
    cdef SVRemap r 

    matrix = numpy.array(M, dtype='int32')
    I[0] = 0
    I[1] = 0
    I[2] = 0
    svremap_init(&r, &matrix[0, 0])
    cdef numpy.intp_t i
    for i in range(x.shape[0]):
        svremap_apply(&r, &x[i, 0], &tmp[0], I)
        y[i, 0] = tmp[0]
        y[i, 1] = tmp[1]
        y[i, 2] = tmp[2]
    return output

