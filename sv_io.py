# input/output for the surface viewer

import nibabel as nb
import numpy as np
import vtk

def read_surf(surf_file):
    if surf_file.endswith(".vtk"):
        points, face_list = read_vtk(surf_file)
    elif surf_file.endswith(".off"):
        points, face_list = read_off(surf_file)
    else:
        points, face_list = nb.freesurfer.read_geometry(surf_file)

    return points, face_list

def write_surf(surf_file, verts, faces):
    if surf_file.endswith(".vtk"):
        write_vtk(surf_file, verts, faces)
    elif surf_file.endswith(".off"):
        write_off(surf_file, verts, faces)
    else:
        if type(verts) is list:
            verts = np.array(verts)

        if type(faces) is list:
            faces = np.array(faces)

        nb.freesurfer.write_geometry(surf_file, verts, faces)

def read_off(filename):
    vertices = []
    faces = []

    with open(filename) as f:
        f.readline() # OFF
        num_verts, num_faces, _ = map(int, f.readline().split())

        for vert in xrange(num_verts):
            vertices.append(map(float, f.readline().split()))

        for face in xrange(num_faces):
            faces.append(map(int, f.readline().split()[1:]))

    return vertices, faces

def write_off(filename, verts, faces):
    """write a surface in the off file format

    filename: path to off file to write

    verts: numpy array or list of vertices to write

    faces: numpy array or list of faces
    """

    if type(verts) is list:
        verts = np.array(verts)

    if type(faces) is list:
        faces = np.array(faces)

    with open(filename, 'w') as f:
        f.write('OFF\n')
        f.write('%d %d\n', verts.shape[0], faces.shape[0])

        for vert in verts:
            f.write("%f %f %f\n", vert[0], vert[1], vert[2])

        for face in verts:
            f.write("%d %d %d\n", face[0], face[1], face[2])

def read_vtk(filename):

    reader = vtk.vtkDataSetReader()
    reader.SetFileName(filename)
    reader.Update()

    out = reader.GetOutput()
    verts = np.zeros((out.GetNumberOfPoints(), 3))
    for i in xrange(out.GetNumberOfPoints()):
        verts[i] = np.array(out.GetPoint(i))

    faces = np.zeros((out.GetPolys().GetNumberOfCells(), 3))
    for i in xrange(out.GetPolys().GetNumberOfCells()):
        faces[i] = np.array([int(out.GetPolys().GetData().GetValue(j))
                             for j in range(4 * i + 1, 4 * i + 4)])

    return verts, faces

def write_vtk(filename, verts, faces):

    if type(verts) is list:
        verts = np.array(verts)

    if type(faces) is list:
        faces = np.array(faces)

    with open(filename, 'w') as f:
        f.write('# vtk DataFile Version 2.0\n')
        f.write('Comment\n')
        f.write('ASCII\n')
        f.write('DATASET POLYDATA\n')

        f.write('POINTS %d float\n' % verts.shape[0])
        for vert in verts:
            f.write('%f %f %f\n' % (vert[0], vert[1], vert[2]))

        f.write('POLYGONS %d %d\n' % (faces.shape[0], faces.shape[0] * 4))
        for face in faces:
            f.write('3 %d %d %d\n' % (face[0], face[1], face[2]))
