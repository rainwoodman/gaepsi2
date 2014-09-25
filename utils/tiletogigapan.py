import bigfilepy
import numpy

from itertools import product

def tilename(i, j, imax, jmax):
    t = 1
    r = []
    n = 0
    while t < imax or t < jmax:
        if (i & t) != 0:
            n += 1
        if (j & t) != 0:
            n += 2
        t = t * 2
        r.append('0%s' % n)
    return '/'.join(r)

print tilename(3, 4, 6, 16)
