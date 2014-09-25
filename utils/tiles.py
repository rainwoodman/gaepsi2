def tilename(i, j, imax, jmax):
    """ i is y, j is x
        randy: 
             0 1
             2 3 """
    t = 1
    r = []
    while t < imax or t < jmax:
        n = 0
        if (i & t) != 0:
            n += 2
        if (j & t) != 0:
            n += 1
        t = t * 2
        r.append('0%s' % n)
    return '/'.join(r[::-1])

