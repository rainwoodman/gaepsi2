
import gaepsi

import bigfilepy

file = bigfilepy.BigFile('....')
pos = file.open('').read(...)

pos -= ...
pos /= ...
pos %= 1.0

mass = file.open('').read(...)
entropy = file.open('').read(...)
potential = file.open('').read(...)

gaepsi.box.normalize(pos, BoxSize, offset=0.0, out=pos)
gaepsi.svr.remap(pos, M, output=pos)

domain = gaepsi.domain.Grid2D(bleeding=....)

layout = domain.decompose(pos[:, :2])
pos = domain.exchange(layout, pos)
mass = domain.exchange(layout, mass)
entropy = domain.exchange(layout, mass)
potential = domain.exchange(layout, mass)

camera = gaepsi.camera()
camera.view(....)

camera.render(pos, mass, (entropy, potential), output='....')

