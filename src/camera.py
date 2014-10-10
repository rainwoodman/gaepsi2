"""
   Camera support in Gaepsi

   a camera is represented by a 4x4 transformation matrix like OpenGL.

   If it is done correctly, we use
   
   pos[:, :4] dot  matrix for the transformation.

   lookat: build the model/view matrix
   ortho:  the orthogonal perspective matrix

   the full transformation is perspective dot lookat

   apply: apply the transformation to a 3 d position vector,

   the devide coordinate is [-1, 1]. 
"""
import numpy

def ortho(near, far, extent):
    """ set up the zoom by extent=(left, right, top, bottom """
    l, r, b, t = extent
    ortho = numpy.zeros((4,4))
    ortho[0, 0] = 2.0 / (r - l)
    ortho[1, 1] = 2.0 / (t - b)
    ortho[2, 2] = -2.0 / (far - near)
    ortho[3, 3] = 1
    ortho[0, 3] = - (1. * r + l) / (r - l)
    ortho[1, 3] = - (1. * t + b) / (t - b)
    ortho[2, 3] = - (1. * far + near) / (far - near)
    return ortho

def lookat(pos, target, up):
    """ the full transformation matrix is 
        perspective dot lookat
    """
    pos = numpy.asarray(pos)
    target = numpy.asarray(target)
    up = numpy.asarray(up)

    dir = target - pos
    dir = dir / (dir **2).sum() ** 0.5
    side = numpy.cross(dir, up)
    side[:] = side / (side **2).sum() ** 0.5
    up = numpy.cross(side, dir)
    up = up / (up **2).sum() ** 0.5

    m1 = numpy.zeros((4,4))
    m1[0, 0:3] = side
    m1[1, 0:3] = up
    m1[2, 0:3] = -dir
    m1[3, 3] = 1

    tran = numpy.eye(4)
    tran[0:3, 3] = -pos
    m2 = numpy.dot(m1, tran)
    return m2

def apply(matrix, pos, out=None):
    """ 
        apply the transform matrix to pos
        matrix = projection dot view 
        returns new pos in device coordinate

        (-1, 1) x (-1, 1) x (-1, 1)
    """
    if out is None:
        out = numpy.empty_like(pos)
    chunksize = 1024 * 32
    for i in range(0, len(pos), chunksize):
        tmppos = pos[i:chunksize+i]
        tmpout = out[i:chunksize+i]
        tmp = numpy.empty((len(tmppos), 4), dtype='f8')
        tmp[:, 3] = 1.0
        tmp[:, :3] = tmppos
        tmp = numpy.dot(tmp, matrix.T)
        tmpout[:, :] = tmp[:, :3]
    return out

