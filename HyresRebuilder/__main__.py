import sys
import argparse
from .Rebuilder import rebuild

parser = argparse.ArgumentParser(
  prog='HyresRebuilder',
  description='Rebuild atomistic model from HyRes model.'
)
parser.add_argument('input')
parser.add_argument('output')

if __name__ == '__main__':
  rebuild(parser.input, parser.output)
