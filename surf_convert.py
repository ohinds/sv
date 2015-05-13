#!/usr/bin/env python

import sys

from sv_io import read_surf, write_surf

def main(argv):
    if len(argv) != 3:
        print "usage: %s <source surf> <dest surf>" % argv[0]
        return 0

    points, faces = read_surf(argv[1])
    write_surf(argv[2], points, faces)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
