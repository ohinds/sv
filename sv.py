#!/usr/bin/env python

"""Simple visualization for vtk, off, or freesurfer surfaces
"""

import argparse
import nibabel as nb
import numpy as np
import sv_io
import sys
import vtk

class SurfaceInteractor(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, interactor, key_callback, mouse_callback):
        self.interactor = interactor
        self.picker = vtk.vtkCellPicker()
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
            self.interactor.GetRenderWindow().GetRenderers(
                ).GetFirstRenderer()):
            if self.mouse_callback is not None:
                self.mouse_callback(self.picker.GetCellId(),
                                    self.picker.GetPickPosition())

        self.OnMouseMove()


class SurfaceViewer(object):

    def __init__(self, options, points=None, face_list=None, overlays=None):
        self.options = options
        self.surface = None
        self.surface_overlays = []
        self.surface_overlay_data = []
        self.last_overlay_value = None
        self.points = []
        self.point_clouds = []
        self._no_show = False

        if points is None:
            assert face_list is None
            surf_file = options.surface[0]
            overlay_files = options.overlays

            points, face_list = sv_io.read_surf(surf_file)

            overlays = []
            for overlay_file in overlay_files:
                if overlay_file.endswith(".vtk"):
                    overlays.append(np.array(read_scalars(overlay_file)[0]))
                else:
                    try:
                        labels, _, _ = nb.freesurfer.read_annot(
                            overlay_file)
                    except:
                        try:
                            labels = nb.freesurfer.read_w_data(
                                overlay_file)
                        except:
                            labels = nb.freesurfer.read_morph_data(
                                overlay_file)
                    overlays.append(labels)

            if options.combine_overlays:
                overlays = [np.sum(np.vstack(overlays), 0)]

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

        if options.colormap:
            colormap = np.loadtxt(options.colormap)
        else:
            colormap = None

        for data in overlays:
            if data.shape[0] > 0:
                self.add_surface_overlay(data, colormap)


    def set_surface(self, points, faces):
        assert self.surface is None

        poly_data = vtk.vtkPolyData()
        poly_data.SetPoints(points)
        poly_data.SetPolys(faces)
        self.points.append(points)
        self.surface = poly_data

    def add_surface_overlay(self, data, colormap=None):
        scalar = vtk.vtkUnsignedCharArray()
        scalar.SetNumberOfComponents(3)

        if colormap is None:
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
        else:
            for val in data:
                if val < colormap[0][0]:
                    ind = 0
                elif val > colormap[-1][0]:
                    ind = colormap.shape[0] - 1
                else:
                    for ind in xrange(colormap.shape[0] - 1):
                        if val < colormap[ind + 1][0]:
                            break
                scalar.InsertNextTuple3(255 * colormap[ind][1],
                                        255 * colormap[ind][2],
                                        255 * colormap[ind][3])

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

    def over_cell(self, cell_id, position):
        if len(self.surface_overlays) < 1:
            return

        face_vertices = self.surface.GetCell(cell_id).GetPointIds()
        min_dist = float('inf')
        min_point = None
        for i in range(0, 3):
            dist = sum((np.array(self.surface.GetPoint(
                            face_vertices.GetId(i))) -
                            np.array(position))**2)
            if dist < min_dist:
                min_dist = dist
                min_point = face_vertices.GetId(i)

        point_id = min_point

        if (self.last_overlay_value ==
            self.surface_overlay_data[self.current_overlay][point_id]):
            return

        print "%d: %f" % (
            point_id, self.surface_overlay_data[self.current_overlay][point_id])

        self.last_overlay_value = \
            self.surface_overlay_data[self.current_overlay][point_id]


    def show(self):
        if self._no_show:
            return 0

        ren = vtk.vtkRenderer()
        self.ren_win = vtk.vtkRenderWindow()
        self.ren_win.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(self.ren_win)
        iren.SetInteractorStyle(
            SurfaceInteractor(iren, self.change_overlay,
                              self.over_cell if self.options.show_values
                              else None))

        surf_mapper = vtk.vtkDataSetMapper()
        surf_mapper.SetInput(self.surface)
        surf_actor = vtk.vtkActor()

        # TODO: allow interactive toggling of edges.
        # surf_actor.GetProperty().SetEdgeColor(0.5, 0.5, 0.5)
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

def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument('--show-values', '-s', action="store_true",
                        help="Mouse over vertices to show their values")
    parser.add_argument('--combine-overlays', '-c', action="store_true",
                        help="Combine the values of all overlays into one")
    parser.add_argument('--colormap', '-m',
                        help="File containing a colormap definition")
    parser.add_argument('surface', nargs=1,
                        help="VTK, OFF, or FreeSurfer surface file to display")
    parser.add_argument('overlays', nargs='*',
                        help="VTK, LAB, or FreeSurfer files to overlay")


    return parser.parse_args(args)

if __name__ == "__main__":
    SurfaceViewer(parse_args(sys.argv[1:])).show()
    sys.exit(0)
