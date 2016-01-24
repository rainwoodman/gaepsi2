"""
   Camera support in Gaepsi

   a camera is represented by a 4x4 transformation matrix like OpenGL.

   If it is done correctly, we use
   
   pos[:, :4] dot  matrix for the transformation.

   lookat: build the model/view matrix
   ortho:  the orthogonal perspective matrix

   the full transformation is ortho dot lookat

   apply: apply the transformation to a 3 d position vector,

   the devide coordinate is [-1, 1]. 
"""
import numpy
import sharedmem
class perspectivematrix(numpy.ndarray):
    pass

class modelviewmatrix(numpy.ndarray):
    pass

def ortho(near, far, extent):
    """ set up the zoom by extent=(left, right, top, bottom) """
    l, r, b, t = extent
    ortho = numpy.zeros((4,4))
    ortho[0, 0] = 2.0 / (r - l)
    ortho[1, 1] = 2.0 / (t - b)
    ortho[2, 2] = -2.0 / (far - near)
    ortho[3, 3] = 1
    ortho[0, 3] = - (1. * r + l) / (r - l)
    ortho[1, 3] = - (1. * t + b) / (t - b)
    ortho[2, 3] = - (1. * far + near) / (far - near)
    ortho = ortho.view(type=perspectivematrix)
    ortho.extent = extent
    ortho.near = near
    ortho.far = far
    ortho.scale = numpy.array([ortho[0, 0], ortho[1, 1]])
    return ortho

def lookat(pos, target, up):
    """ the full transformation matrix is 
        modelviewmatrix dot lookat
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
    m2 = m2.view(type=modelviewmatrix)
    return m2

def apply(matrix, pos, np=None):
    """ 
        apply the transform matrix to pos
        matrix = projection dot modelview 
        returns new pos in device coordinate

        (-1, 1) x (-1, 1) x (-1, 1)
    """
    shmout = sharedmem.empty_like(pos)
    chunksize = 1024 * 32

    def work(i):
        tmppos = pos[i:chunksize+i]
        tmpout = shmout[i:chunksize+i]
        tmp = numpy.empty((len(tmppos), 4), dtype='f8')
        tmp[..., 3] = 1.0
        tmp[..., :3] = tmppos
        tmp = numpy.dot(tmp, matrix.T)
        tmpout[..., :] = tmp[..., :3]

    with sharedmem.MapReduce(np=np) as pool:
        pool.map(work, range(0, len(pos), chunksize))

    return shmout

def todevice(pos2d, extent, np=None):
    """
        convert to device coordinate
    """
    l, r, b, t = extent

    chunksize = 1024 * 32
    out = sharedmem.empty_like(pos2d)
    def work(i):
        tmp = (pos2d[i:i+chunksize] + 1.0)
        tmp *= 0.5
        tmp[..., 0] *= (r - l)
        tmp[..., 0] += l
        tmp[..., 1] *= (t - b)
        tmp[..., 1] += b
        out[i:i+chunksize] = tmp
    with sharedmem.MapReduce(np=np) as pool:
        pool.map(work, range(0, len(pos2d), chunksize))
    return out

def test():
    pos = numpy.random.uniform(size=(1000, 3)) * 20.
    pos -= 10.
    proj = ortho(0, 20, (-10, 10, -10, 10))
    mcam = lookat((0, 0, -10), (0, 0, 0), (0, 1, 0))
    p2d = apply(proj.dot(mcam), pos)
    print p2d.max(axis=0)
    print p2d.min(axis=0)

if __name__ == '__main__':
    test()
