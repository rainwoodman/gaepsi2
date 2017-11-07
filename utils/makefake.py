import bigfile
import numpy

f = bigfile.BigFile("fakedata", create=True)
header = f.create("header")

pos = numpy.indices((16, 16, 1), dtype='f8').reshape(3, -1).T
pos += [0.5, 0.5, 8.0];
print(pos)
mass = len(pos) * 2 - (numpy.arange(len(pos)) + 1)

sml = numpy.linspace(0.1, 0.3, len(pos))

entropy = numpy.array([10.0] * len(pos))

print(mass.sum(dtype='f8'))
header.attrs['TotNumPart'] = [len(pos), 0, 0, 0, 0, 0]
header.attrs['BoxSize'] = 16.0

blk = f.create("0/Position", ('f8', 3), len(pos))
blk.write(0, pos)

blk = f.create("0/Mass", 'f8', len(mass))
blk.write(0, mass)

blk = f.create("0/SmoothingLength", 'f8', len(pos))
blk.write(0, sml)

blk = f.create("0/Entropy", 'f8', len(entropy))
blk.write(0, sml)
