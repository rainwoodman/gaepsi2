
# usage python -S python-mpi.py yourscript.py args

# this will handle exceptions gracefully with a call to mpiimport.abort()

import mpiimport; mpiimport.install(tmpdir='/tmp', verbose=False, disable=False)
import sys; main = sys.argv[1]; sys.argv = sys.argv[1:]
import traceback
# save stdout in case it gets tampered
stdout = sys.stdout
try:
    execfile(main)
except Exception as e:
    stdout.write(traceback.format_exc())
    stdout.flush()
    mpiimport.abort()
