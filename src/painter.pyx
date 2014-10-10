import numpy
cimport numpy
import cython
cimport cython
from libc.stdlib cimport malloc, free, abort

cdef extern from 'c/spline.c':
    ctypedef double (*gsph_spline_kernel)(double r)
    gsph_spline_kernel gsph_spline_query(char * name)

cdef extern from 'c/sphrasterize.c':
    ctypedef void (*gsph_painter_write)(void * data, int x, int y, double * values)

    char * SPLINE_2D_PROJ_CUBIC "SPLINE_2D_PROJ_CUBIC"
    ctypedef struct GSPHPainter "GSPHPainter":
        int size[2]
        gsph_painter_write write
        gsph_spline_kernel sphkernel
        void * data
        int value
        
    void gsph_painter_init(GSPHPainter * painter, int size[2], 
        gsph_painter_write write, 
        gsph_spline_kernel sphkernel, 
        void * data, 
        int nvalue)
    void gsph_painter_rasterize(GSPHPainter * painter,
        double pos[2], double sml, double * mvalue)

ctypedef fused floatingpos:
    cython.int [:, :]
    cython.long [:, :]
    cython.float [:, :]
    cython.double [:, :]

ctypedef fused floatingdata:
    cython.int [:, :]
    cython.float [:, :]
    cython.double [:, :]

ctypedef fused floatingsml:
    cython.float [:]
    cython.double [:]

ctypedef fused floatingimage:
    cython.float [:, :, :]
    cython.double [:, :, :]

cdef numpy.ndarray A(obj):
    return numpy.asarray(obj)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.overflowcheck(False)
@cython.nonecheck(False)
cdef void _write(floatingimage * image, int x, int y, double * values):
    cdef int i
    # this is not thread safe
    for i in range(image.shape[2]):
#        if values[i] != values[i]:
#            abort()
        image[0][y, x, i] += values[i]

def paint(pos, sml, data, image):
    """ paint on 
          image shape(height, width, nvalue)
          data shape(:, nvalue)
          pos[0] : 0 .. height
          pos[1] : 0 .. width
    """
    return _paint(A(pos), A(sml), A(data), A(image))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.overflowcheck(False)
@cython.nonecheck(False)
def _paint(floatingpos pos, 
        floatingsml sml,
        floatingdata data,
        floatingimage image):
    cdef GSPHPainter painter
    cdef numpy.intp_t i
    cdef int j
    cdef int size[2]
    size[0] = image.shape[0]
    size[1] = image.shape[1]
    cdef int nvalue = image.shape[2]
    assert data.shape[1] == image.shape[2]
    cdef gsph_spline_kernel sphkernel
    sphkernel = gsph_spline_query(SPLINE_2D_PROJ_CUBIC)

    # watch out bigger than memoryslice object 
    # hack
    cdef floatingimage * msl = <floatingimage*>malloc(4096) 
    msl[0] = image
    gsph_painter_init(&painter, size, 
        <gsph_painter_write>_write[floatingimage], 
        sphkernel, 
        <floatingimage *>msl, image.shape[2])

    cdef double * mvalue = <double*>malloc(nvalue * sizeof(double))
    cdef double dpos[2]
    cdef double dsml
    for i in range(0, pos.shape[0]):
        dpos[0] = pos[i, 0]
        dpos[1] = pos[i, 1]
        dsml = sml[i]
        for j in range(nvalue):
            mvalue[j] = data[i, j]
        gsph_painter_rasterize(&painter, dpos, dsml, mvalue)
    free(mvalue)
    free(msl)
