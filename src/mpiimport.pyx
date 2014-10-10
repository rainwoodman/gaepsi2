include "libmpi.pxd"
import cPickle

bytescomm = 0

cdef class Comm(object):
    cdef MPI_Comm comm
    cdef readonly int rank
    cdef readonly int size

    def barrier(self):
        MPI_Barrier(self.comm)

    def bcast(self, obj, root=0):
        global bytescomm
        cdef int n
        cdef bytes buf

        if self.rank == root:
            buf = cPickle.dumps(obj, 2)
            n = len(buf)
        MPI_Bcast(&n, 1, MPI_INT, root, self.comm)
        bytescomm = bytescomm + n
        if self.rank != root:
            buf = bytes(' ' * n)

        MPI_Bcast(<char*>buf, n, MPI_BYTE, root, self.comm)
        return cPickle.loads(buf)
cdef bind(MPI_Comm comm):
    self = Comm()
    self.comm = comm
    MPI_Comm_rank(self.comm, &self.rank)
    MPI_Comm_size(self.comm, &self.size)
    return self

MPI_Init(NULL, NULL)
COMM_WORLD = bind(MPI_COMM_WORLD)

import imp
import sys
import posix

__all__ = ['install', 'COMM_WORLD']

_tmpdir = '/tmp'
_tmpfiles = []
d = {
        imp.PY_SOURCE: "source",
        imp.PY_COMPILED: "compiled",
        imp.PKG_DIRECTORY: "pkg",
        imp.C_BUILTIN: "builtin",
        imp.C_EXTENSION: "extension",
        imp.PY_FROZEN: "frozen"}
cdef class Profiler:
    cdef readonly double time
    cdef readonly object title
    cdef double now
    cdef int count
    def __init__(self, name):
        self.title = name
        self.time = 0
        self.count = 0
    def start(self):
        self.now = MPI_Wtime()
    def end(self):
        self.time += MPI_Wtime() - self.now
        self.count = self.count + 1
    def __str__(self):
        return '%s: %g (%d)' % (self.title, self.time, self.count)

tio = Profiler('IO')
tload = Profiler('LOAD')
tloadlocal = Profiler('LOADDirect')
tfind = Profiler('FIND')
tcomm = Profiler('COMM')
tloadfile = Profiler('LOADFile')
tall = Profiler('ALL')

def tempnam(dir, prefix, suffix):
    l = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    s = posix.urandom(16)
    u = ''.join([l[ord(a) % len(l)] for a in s])
    return dir + '/' + prefix + u + suffix

def mkstemp(dir='', suffix='', prefix='', mode='w+', fmode=0600):
    i = 0
    while i < 100:
        fn = tempnam(dir, prefix, suffix)
        try:
            fd = posix.open(fn, posix.O_CREAT | posix.O_EXCL, fmode)
        except OSError:
            i = i + 1
            continue
        f = open(fn, mode)
        posix.close(fd)
        return f
    raise OSError("failed to create a tempfile");

def loadcextensionfromstring(fullname, string, pathname, description):
#    try:
        tio.start()
        with mkstemp(dir=_tmpdir, prefix=fullname.split('.')[-1] + '-', suffix=description[0]) as file:
            file.write(string)
            _tmpfiles.append(file.name)
            name = file.name
        tio.end()

        if _verbose:
            print 'module', fullname, 'using', name, 'via mpi'

        #with open(name, mode=description[1]) as f2:
        tio.start()
        f2 = open(name, mode=description[1])
        tio.end()
        #print file, pathname, description
        tload.start()
        mod = imp.load_module(fullname, f2, pathname, description)
        tload.end()
        if description[-1] == imp.C_EXTENSION:
            mod.filehandle = f2
        else:
            f2.close()
        #print mod
        tio.start()
        posix.unlink(name)
        tio.end()
        return mod 
#    except Exception as e:
#        print 'exception', e

if hasattr(sys, 'exitfunc'):
    oldexitfunc = sys.exitfunc
else:
    oldexitfunc = lambda : None

def cleanup():
    MPI_Finalize()
    return
    global _tmpfiles
    for f in _tmpfiles:
    #    print 'removing', f
        posix.unlink(f)
    _tmpfiles = []
    oldexitfunc()

sys.exitfunc = cleanup

class Loader(object):
    def __init__(self, file, pathname, description):
        self.file = file
        self.pathname = pathname
        self.description = description
    def load_module(self, fullname):
        if not _disable and self.file:
            if self.description[-1] == imp.PY_SOURCE:
                mod = sys.modules.setdefault(fullname,imp.new_module(fullname))
                mod.__file__ = "<%s>" % self.__class__.__name__
                mod.__package__ = fullname.rpartition('.')[0]
                if _verbose:
                    print 'module', fullname, 'using string'
                code = compile(self.file, self.pathname, 'exec', 0, 1)
                exec code in mod.__dict__
#                mod = loadcextensionfromstring(fullname, self.file, self.pathname, self.description) 
            elif self.description[-1] == imp.C_EXTENSION:
                #print "loading extension"
                mod = loadcextensionfromstring(fullname, self.file, self.pathname, self.description) 
            else:
                if _verbose:
                    print 'module', fullname, 'using', self.file
                tio.start()
                self.file = open(self.file, self.description[1])
                tio.end()
                tloadfile.start()
                mod = imp.load_module(fullname, self.file, self.pathname, self.description)
                tloadfile.end()
        else:
            tloadlocal.start()
            mod = imp.load_module(fullname, self.file, self.pathname, self.description)
            tloadlocal.end()
        mod.__loader__ = self
        return mod

class Finder(object):
    def __init__(self, comm):
        self.comm = comm
        self.rank = comm.rank
    def find_module(self, fullname, path=None):
        file, pathname, description = None, None, None
        if _disable:
            tfind.start()
            try:
                file, pathname, description = imp.find_module(fullname, path)
            except ImportError as e:
                file = e
                pass
            tfind.end()
            if self.rank == 0:
                #print fullname, file, pathname, 'disable'
                pass
        else:
            if self.rank == 0:
                tfind.start()
                try:
                    file, pathname, description = imp.find_module(fullname, path)
                except ImportError as e:
                    file = e
                tfind.end()
                #print fullname, file, pathname
                if not isinstance(file, Exception):
                    tio.start()
                    if file:
                        if description[-1] == imp.PY_SOURCE:
                            #print 'finding python module', file.name
                            s = file.read()
                            file.close()
                            file = s
                        elif description[-1] == imp.C_EXTENSION:
                            #print 'finding extension', file.name
                            s = file.read()
                            file.close()
                            file = s
                        else:
                            #print 'finding file by name', d[description[-1]]
                            file = file.name
                    tio.end()
            tcomm.start()
            file, pathname, description = self.comm.bcast((file, pathname, description))
            tcomm.end()

        if isinstance(file, Exception):
            return None
        return Loader(file, pathname, description)

def install(comm=COMM_WORLD, tmpdir='/tmp', verbose=False, disable=False):
    tall.start()
    global _tmpdir
    global _verbose
    global _disable
    _verbose = verbose or int(posix.environ.get('PYTHON_MPIIMPORT_VERBOSE', 0)) == 1
    _disable = disable
    _tmpdir = tmpdir
    sys.meta_path.append(Finder(comm))

    if sys.flags.no_site:
        import site
        site.main0()

        import sysconfig
        import _sysconfigdata
        import re

        if comm.rank == 0:
            site.main1()
        sys.path = comm.bcast(sys.path)
        site.main2()
