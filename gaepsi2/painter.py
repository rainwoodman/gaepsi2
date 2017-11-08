from . import _painter
import sharedmem
import numpy

def paint(pos, sml, data, shape, mask=None, np=0):
    """ paint on 
          paint (pos, sml, data, image)

          data is a list for quantities to paint per channel
          pos[0] : 0 .. height
          pos[1] : 0 .. width
          so remember to transpose using imshow.

        returns (nchan, shape[0], shape[1])
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


