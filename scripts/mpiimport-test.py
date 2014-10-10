import sys
import os.path
d = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.append(d)

import site
print __name__
print site.__mpisite__

