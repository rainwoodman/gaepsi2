import bigfile
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import LogNorm, Normalize
from matplotlib import cm

f = bigfile.BigFile('image2')
h = f.open('header')
tp = h.attrs['TilePadding']
w = f.open('W')
v = f.open('V')

ntile = h.attrs['NTile']
tiles = w.read(0, -1).reshape(ntile[0], ntile[1], tp, tp)
imgw = tiles.transpose((0,2,1,3)).copy().reshape(ntile[0] * tp, ntile[1] * tp)
tiles = v.read(0, -1).reshape(ntile[0], ntile[1], tp, tp)
imgv = tiles.transpose((0,2,1,3)).copy().reshape(ntile[0] * tp, ntile[1] * tp)

print(imgw.sum(dtype='f8'))

imgv /= imgw

fig = Figure(figsize=(ntile[1], ntile[0]), dpi=tp)

ax = fig.add_axes([0, 0, 1, 1])


import numpy
vmax = numpy.nanmax(imgw)
imgw *= (300. / vmax)
numpy.arcsinh(imgw, out=imgw)
imgw /= 5.
brightness = imgw

vmean = numpy.nanmean(imgv)
vstd = numpy.nanstd(imgv)
color = cm.coolwarm(Normalize(vmin=vmean-vstd * 2, vmax=vmean+vstd * 2)(imgv))

color[..., :3] *= brightness[..., None]

color = color.reshape(-1, 4)

bad = (color[..., :3] > 1.0).any(axis=-1)
print(bad.sum(), color.shape)
color[..., :3][bad] /= color[..., :3][bad].max(axis=-1)[..., None]

color = color.reshape(imgw.shape[0], imgw.shape[1], 4)

print numpy.histogram(brightness.ravel())
print numpy.histogram(color[..., :3].max(axis=-1).ravel())
ax.imshow(color)
#ax.imshow(img)

canvas = FigureCanvasAgg(fig)
fig.savefig('image2.png')
