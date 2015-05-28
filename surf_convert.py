#!/usr/bin/env python

import argparse
import nibabel
import numpy
import sys

from sv_io import read_surf, write_surf

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('src_surf')
    parser.add_argument('dst_surf')
    parser.add_argument('--transform-volume', '-v',
                        help="Volume to find affine transformation")
    parser.add_argument('--transform', '-t',
                        help="Transformation to apply. Choices are:\n"
                        "  vox2ras [default],"
                        "  inv-vox2ras,"
                        "  vox2ras-tkr,"
                        "  inv-vox2ras-tkr,"
                        "  vox2mipav,"
                        "  inv-vox2mipav"
                        )

    return parser.parse_args()

def main(argv):
    options = parse_args()

    points, faces = read_surf(options.src_surf)

    if options.transform_volume is not None:
        try:
            vol = nibabel.load(options.transform_volume)

            new_points = None
            trans = numpy.eye(4)
            if options.transform is None or options.transform.endswith("vox2ras"):
                trans = numpy.matrix(vol.get_affine())

            elif options.transform.endswith("vox2ras-tkr"):
                trans = numpy.matrix(numpy.zeros((4, 4)))
                vox_shapes = vol.get_shape()
                vox_sizes = vol.get_header().get_zooms()
                trans[0, 0] = -vox_sizes[0]
                trans[1, 2] = vox_sizes[2]
                trans[2, 1] = -vox_sizes[1]

                # TODO: This seems to be in agreement with the
                # freesurfer coordinate system info reported at
                # https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems,
                # but does not agree with the --vox2ras-tkr reporting
                # of mri_info.
                trans[0, 3] = vox_shapes[0] / 2.
                trans[1, 3] = -vox_shapes[2] / 2.
                trans[2, 3] = vox_shapes[1] / 2.
                trans[3, 3] = 1

            elif options.transform.endswith("vox2mipav"):
                vox_shapes = vol.get_shape()
                vox_sizes = vol.get_header().get_zooms()

                new_points = numpy.empty((points.shape))
                for i, p in enumerate(points):
                    new_points[i, 0] = vox_sizes[0] * p[0]
                    new_points[i, 1] = vox_sizes[1] * ((vox_shapes[1] - 1) - p[1])
                    new_points[i, 2] = vox_sizes[2] * ((vox_shapes[2] - 1) - p[2])

            else:
                raise(ValueError,
                      "ERROR: unknown transform type %s" % options.transform)

            if options.transform.startswith("inv"):
                trans = trans.I

            if new_points is None:
                print "Applying transform:"
                print trans
                new_points = numpy.empty((points.shape))
                for i, p in enumerate(points):
                    pl = list(p)
                    pl.append(1)
                    np = trans * numpy.array([pl]).T
                    new_points[i] = numpy.array(np)[0:3, 0]

            points = new_points
        except:
            print "ERROR: could not load %s" % options.transform_volume
            raise

    write_surf(options.dst_surf, points, faces)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
