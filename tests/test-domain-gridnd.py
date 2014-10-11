import sys
import os.path

d = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.append(d)
import numpy
import domain

fakecomm = lambda : None

fakecomm.size = 9
fakecomm.Alltoall = lambda a, b: None
grid = [
        [0, 3, 6, 9] for dir in [0, 1]
        ]
pos = numpy.array(list(numpy.ndindex((10, 10))))

fakecomm.rank = 0

def inspect(layout):
    art = numpy.zeros((10, 10, 2), 'c1')
    for t, chunk in enumerate(numpy.split(pos[layout.indices], 
        layout.sendcounts.cumsum()[:fakecomm.size - 1], axis=0)):
        art.fill('x')
        for i, j in numpy.ndindex(10, 10):
            rank = (i //3 ) * 3 + j // 3
            if i // 3 >= 3: continue
            if j // 3 >= 3: continue
            art[i, j, 0] = '%d' % rank

        for p in pos:
            art[p[0], p[1], 1] = ' '
        print 'sending to', t, 'counts', layout.sendcounts[t]
        for p in chunk:
            art[p[0], p[1], 1] = '%d' % t
        print art.view(dtype='S2').reshape(10, 10)

def test1():
    """ no bleeding, no periodic """
    dcop = domain.GridND(grid, 
            comm=fakecomm,
            periodic=False)
    layout = dcop.decompose(pos, bleeding=0)
    assert (layout.sendcounts == 9).all()

def test2():
    """ no bleeding, periodic """
    dcop = domain.GridND(grid, 
            comm=fakecomm,
            periodic=True)
    layout = dcop.decompose(pos, bleeding=0)
    #inspect(layout)
    assert (layout.sendcounts == [16, 12, 12, 12, 9, 9, 12, 9, 9]).all()

def test3():
    """ with bleeding, no periodic"""
    dcop = domain.GridND(grid, 
            comm=fakecomm,
            periodic=False)
    layout = dcop.decompose(pos, bleeding=1)
#    inspect(layout)
    assert (layout.sendcounts == [16, 20, 20, 20, 25, 25, 20, 25, 25]).all()

def test4():
    """ with bleeding, periodic"""
    dcop = domain.GridND(grid, 
            comm=fakecomm,
            periodic=True)
    layout = dcop.decompose(pos, bleeding=1)
    #inspect(layout)
    assert (layout.sendcounts == [36, 30, 36, 30, 25, 30, 36, 30, 36]).all()


test1()
test2()
test3()
test4()
