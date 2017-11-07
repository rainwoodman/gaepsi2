import numpy
cimport numpy
import cython
cimport cython
from libc.stdlib cimport malloc, free, abort

cdef extern from 'c/kernel.c':
    ctypedef (void*) GSPHKernel
    GSPHKernel gaussian_kernel

cdef extern from 'c/sphrasterize.c':
    ctypedef struct GSPHImage "GSPHImage":
        int size[3]
        numpy.intp_t strides[3]
        void * data
        
    void gsph_image_init(GSPHImage * image, 
        char * dtype,
        int size[3],
        numpy.intp_t strides[3],
        void * data)
    void gsph_rasterize(GSPHImage * painter, GSPHKernel kernel,
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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.overflowcheck(False)
@cython.nonecheck(False)
def paint(floatingpos pos, 
        floatingsml sml,
        floatingdata data,
        floatingimage image,
        object mask):
    """
        paint to image. Image must be of shape (nc, ny, nx)
        
    """
    cdef:
        GSPHImage img
        numpy.intp_t i
        int usemask
        numpy.ndarray[numpy.uint8_t, ndim=1, cast=True] maskbool
        int j
        int size[3]
        numpy.intp_t strides[3]
        GSPHKernel kernel


    assert data.shape[0] == image.shape[0]

    if mask is None:
        usemask = False
        maskbool = numpy.array([False])
    else:
        usemask = True
        maskbool = mask

        
    # translate the shape and strides 
    size[0] = image.shape[1]
    size[1] = image.shape[2]
    size[2] = image.shape[0] # number of channels
    strides[0] = image.strides[1]
    strides[1] = image.strides[2]
    strides[2] = image.strides[0]

    if image.itemsize == 8:
        gsph_image_init(&img, 'f8', size, strides, &image[0, 0, 0]);
    else:
        gsph_image_init(&img, 'f4', size, strides, &image[0, 0, 0]);

    kernel = gaussian_kernel

    cdef double * mvalue = <double*>malloc(size[2] * sizeof(double))
    cdef double dpos[2]
    cdef double dsml
    for i in range(0, pos.shape[0]):
        if usemask and not maskbool[i]: continue
        dpos[0] = pos[i, 0]
        dpos[1] = pos[i, 1]
        dsml = sml[i]
        for j in range(size[2]):
            mvalue[j] = data[j, i]
        gsph_rasterize(&img, kernel, dpos, dsml, mvalue)

    free(mvalue)
