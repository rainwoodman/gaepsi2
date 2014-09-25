import bigfilepy
import numpy
from matplotlib.figure import Figure
from itertools import product
from tiles import tilename

f = bigfilepy.BigFile('image')
f1 = bigfilepy.BigFile('image2')

w = f.open('W')
v = f.open('V')
h = f.open('header')
tp = h.attrs['TilePadding']
size = h.attrs['ImageSize']
h.attrs['Level'] = 0
ntile = h.attrs['NTile']

h1 = f1.create('header')
h1.attrs['TilePadding'] = tp
h1.attrs['Level'] = h.attrs['Level'] + 1
ntile1 = ntile // 2
h1.attrs['NTile'] = ntile1

w1 = f1.create('W', ('f8', tp * tp), ntile1.prod())
v1 = f1.create('V', ('f8', tp * tp), ntile1.prod())

print ntile1

def readtile(f, h, w, v, i, j):
    ntile = h.attrs['NTile']
    tp = h.attrs['TilePadding']
    if i >= ntile[0] or j >= ntile[1]:
        return numpy.zeros((2, tp, tp), 'f8')
    return (w.read(i * ntile[1] + j, 1)[0].reshape(tp, tp),
            v.read(i * ntile[1] + j, 1)[0].reshape(tp, tp))

def combine(i, j, f, h, w, v):
    i1 = 2 * i
    j1 = 2 * j
    i2 = 2 * i + 1
    j2 = 2 * j + 1
    
    wbig = numpy.zeros((2 * tp, 2 * tp))
    vbig = numpy.zeros((2 * tp, 2 * tp))

    wbig[:tp, :tp], vbig[:tp, :tp] = readtile(f, h, w, v, i1, j1)
    wbig[tp:, :tp], vbig[tp:, :tp] = readtile(f, h, w, v, i2, j1)
    wbig[:tp, tp:], vbig[:tp, tp:] = readtile(f, h, w, v, i1, j2)
    wbig[tp:, tp:], vbig[tp:, tp:] = readtile(f, h, w, v, i2, j2)

    print 'reading', tilename(i1, j1, *ntile)
    print 'reading', tilename(i1, j2, *ntile)
    print 'reading', tilename(i2, j1, *ntile)
    print 'reading', tilename(i2, j2, *ntile)

    wbigr = wbig.reshape(tp, 2, tp, 2)
    vbigr = vbig.reshape(tp, 2, tp, 2)
    return wbigr.sum(axis=(-1, -3)).ravel(), vbigr.sum(axis=(-1, -3)).ravel()

for j, i in product(range(ntile1[1]), range(ntile1[0])):
    wbuf, vbuf = combine(i, j, f, h, w, v)
    print 'writing', tilename(i, j, *ntile1)
    w1.write(i * ntile1[1] + j, wbuf.reshape(1, -1))
    v1.write(i * ntile1[1] + j, vbuf.reshape(1, -1))

