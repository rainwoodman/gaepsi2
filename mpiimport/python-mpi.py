
# usage python -S python-mpi.py yourscript.py args

# this will handle exceptions gracefully with a call to mpiimport.abort()

class OptionError(Exception):
    pass

class getopt(object):
    def help(self):
        return "Usage: " + self.argv[0] + " " \
                + \
                " ".join(
                       [ "-" +  o
                        for o in self.optset - self.optset_wa] ) \
                + \
                " ".join(
                        ["-" + o + " OPTARG"
                        for o in self.optset_wa]) \
                + \
                " .... "

    def __init__(self, argv, optstring):
        self.optstring = optstring
        self.argv = argv

        i = 0
        optlist_wa = []
        optlist = []
        last = None

        while i < len(optstring):
            if optstring[i] == ':':
                optlist_wa.append(last)
                last = None
                i = i + 1
                continue
            #optlist_wa.append(optstring[i])
            optlist.append(optstring[i])
            last = optstring[i]
            i = i + 1

        self.optset_wa = set(optlist_wa) 
        self.optset = set(optlist) 

    def __iter__(self):
        optind = 1 
        self.remain = self.argv[1:]
        while optind < len(self.argv):
            opt = self.argv[optind]
            if opt[0] == '-' and len(opt) == 2:
                optopt = opt[1]
                if optopt not in self.optset:
                    raise OptionError(optopt)
                optind = optind + 1
                self.remain = self.argv[optind:]
            else:
                break
            optarg = None
            if optopt in self.optset_wa:
                if optind >= len(self.argv):
                    raise OptionError(optopt)
                optarg = self.argv[optind]
                optind = optind + 1
                self.remain = self.argv[optind:]
            yield (optopt, optarg, optind)

import mpiimport; 
import sys; 

# save stdout in case it gets tampered
stdout = sys.stdout

#parse args first
command = None
verbose = False
opt = getopt(sys.argv, "vc:")
try:
    for optopt, optarg, optind in iter(opt):
        if optopt == 'c':
            command = optarg
        if optopt == 'v':
            verbose = True
except OptionError as e:
    if mpiimport.COMM_WORLD.rank == 0:
        print opt.help()
    raise SystemExit
sys.argv = opt.remain

# install the hook
mpiimport.install(tmpdir='/tmp', verbose=verbose, disable=False)

import traceback
try:
    if command:
        exec(command)
    else:
        main = opt.remain[0]
        execfile(main)

except Exception as e:
    stdout.write(traceback.format_exc())
    stdout.flush()
    mpiimport.abort()
