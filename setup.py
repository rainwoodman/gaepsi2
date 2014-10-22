from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy
def myext(*args):
    return Extension(*args, include_dirs=["./", numpy.get_include()])
extensions = [
        myext("gaepsi.svr", ["src/svr.pyx"]),
        myext("gaepsi._domain", ["src/_domain.pyx"]),
        myext("gaepsi.painter", ["src/painter.pyx"]),
        myext("mpiimport.mpiimport", ["mpiimport/mpiimport.pyx"]),
        ]

setup(
    name="bigfilepy", version="0.1",
    author="Yu Feng",
    description="python binding of BigFile, a peta scale IO format",
    package_dir = {'gaepsi': 'src'},
    install_requires=['cython', 'numpy'],
    packages= ['gaepsi'],
    requires=['numpy'],
    ext_modules = cythonize(extensions)
)

