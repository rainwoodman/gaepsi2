import numpy

def test_svr():
    from gaepsi2 import svr
    M=[[1, 1, 0], [0, 1, 0], [0, 0, 1]];
    print(svr.remap([[0, 0, 1.]], M))

def test_painter():
    from gaepsi2 import painter
    pos = numpy.array([
            numpy.arange(0, 1),
            numpy.arange(0, 1),
            ]).T
    print(pos)
    sml = numpy.ones(len(pos)) * 4
    data = numpy.ones(len(pos))
    print(painter.paint(pos, sml, [data], (5, 5)))
    print(painter.paint(pos + 0.5, sml, [data], (5, 5)))


def test_camera():
    from gaepsi2 import camera

    pos = numpy.random.uniform(size=(1000, 3)) * 20.
    pos -= 10.
    proj = camera.ortho(0, 20, (-10, 10, -10, 10))
    mv = camera.lookat((0, 0, -10), (0, 0, 0), (0, 1, 0))
    matrix = camera.matrix(proj, mv)
    p2d = camera.apply(matrix, pos)
