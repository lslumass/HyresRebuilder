import sys
from .Rebuilder import rebuild

inp = sys.argv[1]
out = sys.argv[2]

if __name__ == '__main__':
  rebuild(inp, out)
