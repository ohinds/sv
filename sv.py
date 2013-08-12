#!/usr/bin/env python

"""This is for ohinds only. Don't complain if you aren't ohinds.
"""

import nibabel as nb
import numpy as np
import sys
import vtk

from mindboggle.utils.io_vtk import read_faces_points, read_scalars

class SurfaceInteractor(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, interactor, key_callback, mouse_callback):
        self.interactor = interactor
        self.picker = vtk.vtkPointPicker()
        self.key_callback = key_callback
        self.mouse_callback = mouse_callback
        self.AddObserver("KeyPressEvent", self.key_pressed)
        self.AddObserver("MouseMoveEvent", self.mouse_moved)

    def key_pressed(self, obj, event):
        key = self.interactor.GetKeyCode()

        if key == ",":
            self.key_callback(False)
        elif key == ".":
            self.key_callback(True)

    def mouse_moved(self, obj, event):
        coords = self.interactor.GetEventPosition()

        if self.picker.Pick(
            coords[0], coords[1], 0,
            self.interactor.GetRenderWindow().GetRenderers().GetFirstRenderer()):
            self.mouse_callback(self.picker.GetPointId())

        self.OnMouseMove()


class SurfaceViewer(object):

    def __init__(self, *args):

        self.surface = None
        self.surface_overlays = []
        self.surface_overlay_data = []
        self.last_overlay_value = None
        self.points = []
        self.point_clouds = []

        if len(args) == 2:
            surf_file = args[0]
            overlay_files = args[1]

            if surf_file.endswith(".vtk"):
                face_list, points, npoints = read_faces_points(surf_file)
            else:
                (points, face_list) = nb.freesurfer.read_geometry(surf_file)
                npoints = len(points)

            overlays = []
            for overlay_file in overlay_files:
                if overlay_file.endswith(".vtk"):
                    overlays.append(np.array(read_scalars(overlay_file)[0]))
                else:
                    labels, _, _ = nb.freesurfer.read_annot(overlay_file)
                    overlays.append(labels)


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
            if data.shape[0] > 0:
                self.add_surface_overlay(data)

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
        self.surface_overlay_data.append(data)

    def add_point_cloud(self, points, indices=None, rgb=(255, 255, 0)):
        point_data = vtk.vtkPoints()
        cell_data = vtk.vtkCellArray()
        scalar_data = vtk.vtkUnsignedCharArray()
        scalar_data.SetNumberOfComponents(3)

        for idx in indices if indices is not None else xrange(len(points)):
            cell_data.InsertNextCell(1)
            cell_data.InsertCellPoint(point_data.InsertNextPoint(points[idx]))
            scalar_data.InsertNextTuple3(*rgb)

        point_poly = vtk.vtkPolyData()
        point_poly.SetPoints(point_data)
        point_poly.SetVerts(cell_data)
        point_poly.GetPointData().SetScalars(scalar_data)
        self.point_clouds.append(point_poly)

    def add_point_cloud_from_surf_indices(self, indices, rgb=(255, 255, 0)):
        cell_data = vtk.vtkCellArray()
        scalar_data = vtk.vtkUnsignedCharArray()
        scalar_data.SetNumberOfComponents(3)

        for idx in indices:
            cell_data.InsertNextCell(1)
            cell_data.InsertCellPoint(idx)
            scalar_data.InsertNextTuple3(*rgb)

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

        print "showing overlay %d" % self.current_overlay

        self.surface.GetPointData().SetScalars(
            self.surface_overlays[self.current_overlay])
        self.surface.Modified()
        self.ren_win.Render()

    def over_point(self, point_id):
        if len(self.surface_overlays) < 1:
            return

        if (self.last_overlay_value ==
            self.surface_overlay_data[self.current_overlay][point_id]):
            return

        print "%d: %f" % (
            point_id, self.surface_overlay_data[self.current_overlay][point_id])

        self.last_overlay_value = \
            self.surface_overlay_data[self.current_overlay][point_id]


    def show(self):
        ren = vtk.vtkRenderer()
        self.ren_win = vtk.vtkRenderWindow()
        self.ren_win.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(self.ren_win)
        iren.SetInteractorStyle(
            SurfaceInteractor(iren, self.change_overlay, self.over_point))

        surf_mapper = vtk.vtkDataSetMapper()
        surf_mapper.SetInput(self.surface)
        surf_actor = vtk.vtkActor()
        # surf_actor.GetProperty().SetEdgeColor(1, 1, 1)
        # surf_actor.GetProperty().EdgeVisibilityOn()
        surf_actor.SetMapper(surf_mapper)
        ren.AddActor(surf_actor)

        for cloud in self.point_clouds:
            pt_mapper = vtk.vtkDataSetMapper()
            pt_mapper.SetInput(cloud)
            pt_mapper.SetScalarVisibility(1)
            pt_actor = vtk.vtkActor()
            pt_actor.GetProperty().SetPointSize(5)
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
