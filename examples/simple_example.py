import sys
from HyresRebuilder import Rebuilder

hyres_pdb = sys.argv[1]
at_pdb = sys.argv[2]

Rebuilder(hyres_pdb, at_pdb)

