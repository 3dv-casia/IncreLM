import numpy as np
from .vis_utils import test_point_inside_ranges, test_line_inside_ranges


def pyvista_vis_3d_lines(lines, img_hw=(600, 800), width=2, ranges=None, scale=1.0):
    """
    Visualize a 3D line map with `PyVista <https://docs.pyvista.org/version/stable/>`_

    Args:
        lines (list[:class:`limap.base.Line3d`]): The 3D line map
        width (float, optional): width of the line
    """
    import pyvista as pv
    plotter = pv.Plotter(window_size=[img_hw[1], img_hw[0]])
    for line in lines:
        if ranges is not None:
            if not test_line_inside_ranges(line, ranges):
                continue
        plotter.add_lines(line.as_array() * scale, '#000000', width=width)

    # Helpful options for capturing screenshots of different line mapping methods from the same viewpoint.
    set_camera_position = False
    save_camera_position = False

    if set_camera_position:
        # You can assign `camera_position` in terms of `camera_position.txt`. An example is as follows.
        camera_position = [(-3.7999621397657775, -15.216386879718808, 10.954815430513474),
                           (1.8069981276499998, -2.47011542035, 1.13601308155),
                           (0.2172003156476655, 0.5339318829465624, 0.8171540658013596)]
        plotter.camera_position = camera_position
        plotter.show()
    else:
        plotter.show()
        if save_camera_position:
            # Save the camera position when closing the window.
            camera_position = plotter.camera_position
            with open('camera_position.txt', 'w') as file:
                file.write(str(camera_position))


def open3d_add_points(w, points, color=[0.0, 0.0, 0.0], psize=1.0, name="pcd", ranges=None, scale=1.0):
    if np.array(points).shape[0] == 0:
        return w
    import open3d as o3d
    o3d_points, o3d_colors = [], []
    for idx in range(np.array(points).shape[0]):
        if ranges is not None:
            if not test_point_inside_ranges(points[idx], ranges):
                continue
        o3d_points.append(points[idx] * scale)
        o3d_colors.append(color)
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(np.stack(o3d_points))
    pcd.colors = o3d.utility.Vector3dVector(np.stack(o3d_colors))
    mat = o3d.visualization.rendering.MaterialRecord()
    mat.shader = "defaultUnlit"
    mat.point_size = psize
    w.add_geometry(name, pcd, mat)
    return w


def open3d_get_line_set(lines, color=[0.0, 0.0, 0.0], width=2, ranges=None, scale=1.0):
    import open3d as o3d
    o3d_points, o3d_lines, o3d_colors = [], [], []
    counter = 0
    for line in lines:
        if ranges is not None:
            if not test_line_inside_ranges(line, ranges):
                continue
        o3d_points.append(line.start * scale)
        o3d_points.append(line.end * scale)
        o3d_lines.append([2 * counter, 2 * counter + 1])
        counter += 1
        o3d_colors.append(color)
    line_set = o3d.geometry.LineSet()
    line_set.points = o3d.utility.Vector3dVector(o3d_points)
    line_set.lines = o3d.utility.Vector2iVector(o3d_lines)
    line_set.colors = o3d.utility.Vector3dVector(o3d_colors)
    return line_set


def open3d_add_line_set(w, lines, color=[0.0, 0.0, 0.0], width=2, name="lineset", ranges=None, scale=1.0):
    import open3d as o3d
    o3d_points, o3d_lines, o3d_colors = [], [], []
    counter = 0
    for line in lines:
        if ranges is not None:
            if not test_line_inside_ranges(line, ranges):
                continue
        o3d_points.append(line.start * scale)
        o3d_points.append(line.end * scale)
        o3d_lines.append([2 * counter, 2 * counter + 1])
        counter += 1
        o3d_colors.append(color)
    line_set = o3d.geometry.LineSet()
    line_set.points = o3d.utility.Vector3dVector(o3d_points)
    line_set.lines = o3d.utility.Vector2iVector(o3d_lines)
    line_set.colors = o3d.utility.Vector3dVector(o3d_colors)
    mat = o3d.visualization.rendering.MaterialRecord()
    mat.shader = "unlitLine"
    mat.line_width = width
    w.add_geometry(name, line_set, mat)
    return w


def open3d_get_cameras(imagecols, color=[1.0, 0.0, 0.0], ranges=None, scale_cam_geometry=1.0, scale=1.0):
    import open3d as o3d
    import copy
    cameras = o3d.geometry.LineSet()

    camera_lines = {}
    for cam_id in imagecols.get_cam_ids():
        cam = imagecols.cam(cam_id)
        camera_lines[cam_id] = o3d.geometry.LineSet.create_camera_visualization(
            cam.w(), cam.h(), cam.K(), np.eye(4), scale=0.005 * scale_cam_geometry * scale)
    for img_id in imagecols.get_img_ids():
        camimage = imagecols.camimage(img_id)
        if ranges is not None:
            if not test_point_inside_ranges(camimage.pose.center(), ranges):
                continue
        T = np.eye(4)
        T[:3, :3] = camimage.R()
        T[:3, 3] = camimage.T() * scale
        T = np.linalg.inv(T)
        cam = copy.deepcopy(camera_lines[camimage.cam_id]).transform(T)
        cam.paint_uniform_color(color)
        cameras += cam
    return cameras


def open3d_add_cameras(w, imagecols, color=[1.0, 0.0, 0.0], ranges=None, scale_cam_geometry=1.0, scale=1.0):
    import open3d as o3d
    import copy
    camera_lines = {}
    for cam_id in imagecols.get_cam_ids():
        cam = imagecols.cam(cam_id)
        camera_lines[cam_id] = o3d.geometry.LineSet.create_camera_visualization(
            cam.w(), cam.h(), cam.K(), np.eye(4), scale=0.005 * scale_cam_geometry * scale)
    for img_id in imagecols.get_img_ids():
        camimage = imagecols.camimage(img_id)
        if ranges is not None:
            if not test_point_inside_ranges(camimage.pose.center(), ranges):
                continue
        T = np.eye(4)
        T[:3, :3] = camimage.R()
        T[:3, 3] = camimage.T() * scale
        T = np.linalg.inv(T)
        cam = copy.deepcopy(camera_lines[camimage.cam_id]).transform(T)
        cam.paint_uniform_color(color)
        w.add_geometry(camimage.image_name(), cam)
    return w


def open3d_vis_3d_lines(lines, width=2, ranges=None, scale=1.0):
    """
    Visualize a 3D line map with `Open3D <http://www.open3d.org/>`_

    Args:
        lines (list[:class:`limap.base.Line3d`]): The 3D line map
        width (float, optional): width of the line
    """
    import open3d as o3d
    vis = o3d.visualization.Visualizer()
    vis.create_window(height=1080, width=1920)
    line_set = open3d_get_line_set(lines, width=width, ranges=ranges, scale=scale)
    vis.add_geometry(line_set)
    vis.run()
    vis.destroy_window()
