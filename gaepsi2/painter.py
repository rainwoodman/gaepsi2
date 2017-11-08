from . import _painter
import sharedmem
import numpy

def paint(pos, sml, data, shape, mask=None, np=0):
    """ Use SPH kernel to splat particles to an image.

        Parameters
        ----------

        pos : array_like
          (..., >2) position of particles. Only two first
          column is used. In device coordinate

        data : array_like
          (Nc, ...) or (...). Weight to use for painting.
          Nc channels will be produces on the device.
          if the array is 1d, Nc = 1
        
        sml : array_like
          smoothing length (half of effective size).
          In device coordinate; only correct in isotropic
          cameras

        shape : list, tuple
          (w[0], w[1]) the size of the device.
          shall enclose pos[..., 0] and pos[..., 1]

        mask : array_like, boolean
          If provided, elements with False will not be painted.

        np : int
          number of multiprocessing. 0 for single-processing.
          None for all available cores.

        Returns
        -------
        image: array_like
           (Nc, shape[0], shape[1])

        Notes
        -----
        Remember to transpose for imshow and pmesh to correct put x horizontaly.
    """

    if len(numpy.shape(data)) == 1:
        data = [data]

    with sharedmem.MapReduce(np=np) as pool:
        if pool.np > 0: nbuf = pool.np
        else: nbuf = 1
        buf = sharedmem.empty((nbuf, len(data)) + shape, dtype='f4')
        buf[:] = 0
        chunksize = 1024 * 8

        def work(i):
            sl = slice(i, i + chunksize)
            datas = [d[sl] for d in data]
            if mask is not None: masks = mask[sl]
            else: masks = None
            _painter.paint(pos[sl], sml[sl], numpy.array(datas),
                    buf[pool.local.rank], masks)

        pool.map(work, range(0, len(pos), chunksize))
    return numpy.sum(buf, axis=0)


