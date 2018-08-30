from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy
def myext(*args):
    return Extension(*args, include_dirs=["./", numpy.get_include()],
            extra_compile_args=['-O3'],
            extra_link_args=['-O3']
            )

extensions = [
        myext("gaepsi2.svr", ["gaepsi2/svr.pyx"]),
        myext("gaepsi2._painter", ["gaepsi2/_painter.pyx"])
        ]


def find_version(path):
    import re
    # path shall be a plain ascii text file.
    s = open(path, 'rt').read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              s, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Version not found")

setup(
    name="gaepsi2", version=find_version('gaepsi2/version.py'),
    author="Yu Feng",
    author_email="rainwoodman@gmail.com",
    url="http://github.com/rainwoodman/gaepsi2",
    description="gaepsi2. SPH visualization.",
    install_requires=['cython', 'numpy', 'sharedmem'],
    zip_safe = False,
    packages= ['gaepsi2', 'gaepsi2.tests'],
    ext_modules = cythonize(extensions),
    license='BSD-2-Clause',
)

