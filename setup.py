from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy
def myext(*args):
    return Extension(*args, include_dirs=["./", numpy.get_include()])
extensions = [
        myext("gaepsi2.svr", ["src/svr.pyx"]),
        myext("gaepsi2._painter", ["src/_painter.pyx"]),
        ]

setup(
    name="gaepsi2", version="0.1",
    author="Yu Feng",
    description="gaepsi2",
    package_dir = {'gaepsi2': 'src'},
    install_requires=['cython', 'numpy'],
    packages= ['gaepsi2'],
    requires=['numpy'],
    ext_modules = cythonize(extensions)
)

