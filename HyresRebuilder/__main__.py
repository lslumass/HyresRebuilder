import sys
from .Rebuilder import rebuild


if __name__ == '__main__':
  inp = sys.argv[1]
  out = sys.argv[2]
  rebuild(inp, out)
