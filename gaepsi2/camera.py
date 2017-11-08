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

def fov2extent(fov, aspect, D):
    l = - numpy.tan(fov * 0.5) * aspect * D
    r = - l
    b = - numpy.tan(fov * 0.5) * D
    t = - b
    return (l, r, b, t)

def extent2fov(extent, D):
    l, r, b, t = extent
    aspect = (l - r) /(b - t)
    fov = numpy.arctan2((t - b), D) * 2
    return fov, aspect

def persp(near, far, fov, aspect):
    """ set up the perspect by fov, and aspect. fov in radians."""
    
    # distance cancels out
    l, r, b, t = fov2extent(fov, aspect, 1)
    
    persp = numpy.zeros((4,4))
 
    persp[0, 0] = 2. / (r - l)
    persp[1, 1] = 2. / (t - b)
    persp[2, 2] = - (1. *(far + near)) / (far - near)
    persp[2, 3] = - (2. * far * near) / (far - near)
    persp[0, 2] = (r + l) / (r - l)
    persp[1, 2] = (t + b) / (t - b)
    persp[3, 2] = -1
    persp[3, 3] = 0

    persp = persp.view(type=perspectivematrix)
    persp.near = near
    persp.far = far
    persp.scale = numpy.array([persp[0, 0], persp[1, 1]])

    return persp

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
    m2.up = up
    m2.side = side
    return m2

def apply(matrix, pos, np=None):
    """ 
        apply the transform matrix to pos
        matrix = projection dot modelview 
        returns new pos in clip coordinate:

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
        n = tmp[..., 3].copy()

        tmpout[..., :] = tmp[..., :3] / n[..., None]
        
    with sharedmem.MapReduce(np=np) as pool:
        pool.map(work, range(0, len(pos), chunksize))

    return shmout

def clip(xc):
    return ((xc >= -1) & (xc <= 1)).all(axis=-1)

def todevice(xc, extent, np=None):
    """
        convert to device coordinate
    """
    if len(extent) == 2:
        r, t = extent
        l, b = 0, 0
    else:
        l, r, b, t = extent

    chunksize = 1024 * 32
    out = sharedmem.empty((len(xc), 2))
    def work(i):
        tmp = (xc[i:i+chunksize] + 1.0)
        tmp *= 0.5
        tmp[..., 1] *= (t - b)
        tmp[..., 1] += b
        tmp[..., 0] *= (r - l)
        tmp[..., 0] += l
        out[i:i+chunksize] = tmp[:, :2]
    with sharedmem.MapReduce(np=np) as pool:
        pool.map(work, range(0, len(xc), chunksize))

    return out

def test():
    pos = numpy.random.uniform(size=(1000, 3)) * 20.
    pos -= 10.
    proj = ortho(0, 20, (-10, 10, -10, 10))
    mcam = lookat((0, 0, -10), (0, 0, 0), (0, 1, 0))
    p2d = apply(proj.dot(mcam), pos)
    #print p2d.max(axis=0)
    #print p2d.min(axis=0)

if __name__ == '__main__':
    test()
