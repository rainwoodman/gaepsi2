import sys
import os.path

d = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.append(d)

import mpiimport; mpiimport.install(tmpdir='/tmp', verbose=False,
    disable=False)

import site

print site.__mpisite__

