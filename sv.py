#!/usr/bin/env python

import numpy as np
import sys
import vtk

from mindboggle.utils.io_vtk import read_faces_points, read_scalars

class SurfaceInteractor(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, interactor, callback):
        self.interactor = interactor
        self.callback = callback
        self.AddObserver("KeyPressEvent", self.key_pressed)

    def key_pressed(self, obj, event):
        key = self.interactor.GetKeyCode()

        if key == ",":
            self.callback(False)
        elif key == ".":
            self.callback(True)

class SurfaceViewer(object):

    def __init__(self, surf_file, overlay_files):

        # surface
        face_list, points, npoints = read_faces_points(surf_file)

        vertices = vtk.vtkPoints()
        for pt in points:
            vertices.InsertNextPoint(pt)

        faces = vtk.vtkCellArray()
        for face in face_list:
            faces.InsertNextCell(3)
            for pt in face:
                faces.InsertCellPoint(pt)

        self.poly_data = vtk.vtkPolyData()
        self.poly_data.SetPoints(vertices)
        self.poly_data.SetPolys(faces)

        # overlays
        self.scalars = []
        self.current_overlay = -1
        for overlay_file in overlays:
            data = read_scalars(overlay_file)
            scalar = vtk.vtkUnsignedCharArray()
            scalar.SetNumberOfComponents(3)

            ran = (np.amin(data[0]), np.amax(data[0]))
            mid = (ran[0] + ran[1]) / 2.
            for val in data[0]:
                if val < mid:
                    gr = 255.0 * (val - ran[0]) / (mid - ran[0])
                else:
                    gr = 255.0 * (ran[1] - val) / (ran[1] - mid)
                scalar.InsertNextTuple3(
                    255 * (val - ran[0]) / (ran[1] - ran[0]),
                    gr,
                    255 * (ran[1] - val) / (ran[1] - ran[0]))

            self.scalars.append(scalar)

    def change_overlay(self, move_up):
        if len(self.scalars) < 1:
            return

        if move_up:
            self.current_overlay += 1
            if self.current_overlay >= len(self.scalars):
                self.current_overlay = 0

        else:
            self.current_overlay -= 1
            if self.current_overlay < 0:
                self.current_overlay = len(self.scalars) - 1

        self.poly_data.GetPointData().SetScalars(
            self.scalars[self.current_overlay])
        self.poly_data.Modified()
        self.ren_win.Render()

    def show(self):
        ren = vtk.vtkRenderer()
        self.ren_win = vtk.vtkRenderWindow()
        self.ren_win.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(self.ren_win)

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInput(self.poly_data)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        iren.SetInteractorStyle(SurfaceInteractor(iren, self.change_overlay))
        ren.AddActor(actor)

        ren.ResetCamera()
        iren.Initialize()
        self.change_overlay(True)
        iren.Start()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "usage %s <surf> [<overlay> ...]"
        sys.exit(1)

    surf_file = sys.argv[1]
    overlays = sys.argv[2:] if len(sys.argv) > 2 else []
    SurfaceViewer(surf_file, overlays).show()
    sys.exit(0)
