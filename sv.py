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

    def __init__(self, *args):

        self.surface = None
        self.surface_overlays = []
        self.points = []
        self.point_clouds = []

        if len(args) == 2:
            surf_file = args[0]
            overlay_files = args[1]

            face_list, points, npoints = read_faces_points(surf_file)
            overlays = []
            for overlay_file in overlay_files:
                overlays.append(read_scalars(overlay_file))

        elif len(args) == 3:
            points = args[0]
            face_list = args[1]
            overlays = args[2]

        else:
            raise ValueError("invalid arguments")

        vertices = vtk.vtkPoints()
        for pt in points:
            vertices.InsertNextPoint(pt)

        faces = vtk.vtkCellArray()
        for face in face_list:
            faces.InsertNextCell(3)
            for pt in face:
                faces.InsertCellPoint(pt)

        self.set_surface(vertices, faces)

        # overlays
        self.scalars = []
        self.current_overlay = -1
        for data in overlays:
            if len(data[0]) > 0:
                self.add_surface_overlay(data[0])

    def set_surface(self, points, faces):
        assert self.surface is None

        poly_data = vtk.vtkPolyData()
        poly_data.SetPoints(points)
        poly_data.SetPolys(faces)
        self.points.append(points)
        self.surface = poly_data

    def add_surface_overlay(self, data):
        scalar = vtk.vtkUnsignedCharArray()
        scalar.SetNumberOfComponents(3)

        ran = (np.amin(data), np.amax(data))
        mid = (ran[0] + ran[1]) / 2.
        for val in data:
            if val < mid:
                gr = 255.0 * (val - ran[0]) / (mid - ran[0])
            else:
                gr = 255.0 * (ran[1] - val) / (ran[1] - mid)
            scalar.InsertNextTuple3(
                255 * (val - ran[0]) / (ran[1] - ran[0]),
                gr,
                255 * (ran[1] - val) / (ran[1] - ran[0]))

        self.surface_overlays.append(scalar)

    def add_point_cloud(self, points, indices=None):
        point_data = vtk.vtkPoints()
        cell_data = vtk.vtkCellArray()
        scalar_data = vtk.vtkUnsignedCharArray()
        scalar_data.SetNumberOfComponents(3)

        for idx in indices if indices is not None else xrange(len(points)):
            cell_data.InsertNextCell(1)
            cell_data.InsertCellPoint(point_data.InsertNextPoint(points[idx]))
            scalar_data.InsertNextTuple3(255, 255, 0)

        point_poly = vtk.vtkPolyData()
        point_poly.SetPoints(point_data)
        point_poly.SetVerts(cell_data)
        point_poly.GetPointData().SetScalars(scalar_data)
        self.point_clouds.append(point_poly)

    def add_point_cloud_from_surf_indices(self, indices):
        cell_data = vtk.vtkCellArray()
        scalar_data = vtk.vtkUnsignedCharArray()
        scalar_data.SetNumberOfComponents(3)

        for idx in indices:
            cell_data.InsertNextCell(1)
            cell_data.InsertCellPoint(idx)
            scalar_data.InsertNextTuple3(255, 255, 0)

        point_poly = vtk.vtkPolyData()
        point_poly.SetPoints(self.points[0])
        point_poly.SetVerts(cell_data)
        point_poly.GetCellData().SetScalars(scalar_data)
        self.point_clouds.append(point_poly)


    def change_overlay(self, move_up):
        if len(self.surface_overlays) < 1:
            return

        if move_up:
            self.current_overlay += 1
            if self.current_overlay >= len(self.surface_overlays):
                self.current_overlay = 0

        else:
            self.current_overlay -= 1
            if self.current_overlay < 0:
                self.current_overlay = len(self.surface_overlays) - 1

        self.surface.GetPointData().SetScalars(
            self.surface_overlays[self.current_overlay])
        self.surface.Modified()
        self.ren_win.Render()

    def show(self):
        ren = vtk.vtkRenderer()
        self.ren_win = vtk.vtkRenderWindow()
        self.ren_win.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(self.ren_win)
        iren.SetInteractorStyle(SurfaceInteractor(iren, self.change_overlay))

        surf_mapper = vtk.vtkDataSetMapper()
        surf_mapper.SetInput(self.surface)
        surf_actor = vtk.vtkActor()
        surf_actor.SetMapper(surf_mapper)
        ren.AddActor(surf_actor)

        for cloud in self.point_clouds:
            pt_mapper = vtk.vtkDataSetMapper()
            pt_mapper.SetInput(cloud)
            pt_mapper.SetScalarVisibility(1)
            pt_actor = vtk.vtkActor()
            pt_actor.SetMapper(pt_mapper)
            ren.AddActor(pt_actor)

        ren.ResetCamera()
        iren.Initialize()
        self.change_overlay(True)
        iren.Start()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "usage %s <surf> [<overlay> ...]"
        sys.exit(1)

    surf_file = sys.argv[1]
    overlays = sys.argv[2:] if len(sys.argv) > 2 else [sys.argv[1]]
    SurfaceViewer(surf_file, overlays).show()
    sys.exit(0)
