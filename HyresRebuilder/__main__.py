import sys
import argparse
from .Rebuilder import rebuild
from .RNA_Rebuilder import RNArebuild


def main():
    parser = argparse.ArgumentParser(description="A command-line tool for my package")
    parser.add_argument('command', choices=['hyresrebuilder', 'RNArebuilder'], help="The function to execute")
    args = parser.parse_args()

    if args.command == 'hyresrebuilder':
        rebuild()
    elif args.command == 'RNArebuilder':
        RNArebuild()

if __name__ == '__main__':
  main()
