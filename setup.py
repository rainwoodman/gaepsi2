from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy
def myext(*args):
    return Extension(*args, include_dirs=["./", numpy.get_include()])
extensions = [
        myext("gaepsi2.svr", ["gaepsi2/svr.pyx"]),
        myext("gaepsi2._painter", ["gaepsi2/_painter.pyx"]),
        ]

setup(
    name="gaepsi2", version="0.1",
    author="Yu Feng",
    description="gaepsi2",
    install_requires=['cython', 'numpy'],
    packages= ['gaepsi2'],
    requires=['numpy'],
    ext_modules = cythonize(extensions)
)

