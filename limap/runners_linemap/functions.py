import numpy as np
from tqdm import tqdm
import limap.structures as _structures
import copy


def compute_2d_bipartites_from_colmap(reconstruction, imagecols, all_2d_lines, cfg=dict(), min_support_images=3):
    all_bpt2ds = {}
    cfg_bpt2d = _structures.PL_Bipartite2dConfig(cfg)
    reconstruction_copy = copy.deepcopy(reconstruction)
    colmap_cameras, colmap_images, colmap_points = reconstruction_copy[
        "cameras"], reconstruction_copy["images"], reconstruction_copy["points"]
    print("Start computing 2D bipartites...")

    filtered_point_ids = set()

    for img_id, colmap_image in tqdm(colmap_images.items()):
        n_points = colmap_image.xys.shape[0]
        indexes = np.arange(0, n_points)
        xys = colmap_image.xys
        point3D_ids = colmap_image.point3D_ids
        mask = []
        for point3D_id in point3D_ids:
            if point3D_id >= 0:
                if colmap_points[point3D_id].image_ids.shape[0] >= min_support_images:
                    mask.append(True)
                else:
                    filtered_point_ids.add(point3D_id)
                    mask.append(False)
            else:
                mask.append(False)

        # resize xys if needed
        cam_id = imagecols.camimage(img_id).cam_id
        orig_size = (colmap_cameras[cam_id].width, colmap_cameras[cam_id].height)
        cam = imagecols.cam(cam_id)
        new_size = (cam.w(), cam.h())
        if orig_size != new_size:
            xys[:, 0] = xys[:, 0] * new_size[0] / orig_size[0]
            xys[:, 1] = xys[:, 1] * new_size[1] / orig_size[1]

        # init bpt2d
        bpt2d = _structures.PL_Bipartite2d(cfg_bpt2d)
        bpt2d.init_lines(all_2d_lines[img_id])
        bpt2d.add_keypoints_with_point3D_ids(xys[mask], point3D_ids[mask], indexes[mask])
        all_bpt2ds[img_id] = bpt2d

    points = {}
    for point3d_id, p in tqdm(colmap_points.items()):
        if p.image_ids.shape[0] >= min_support_images:
            points[point3d_id] = p.xyz

    print("[before filtering] Number of SfM points = {}".format(len(colmap_points)))
    print("[after filtering]  Number of SfM points (# support images >= {}) = {}".format(
        min_support_images, len(colmap_points) - len(filtered_point_ids)))

    return all_bpt2ds, points


def filter_tracks_by_num_support_images(linetracks, th_num_support_images=4):
    new_linetracks = []
    for linetrack in linetracks:
        if linetrack.count_images() >= th_num_support_images:
            new_linetracks.append(linetrack)
    print("[before filtering] num tracks = {}".format(len(linetracks)))
    print("[after filtering] num tracks = {} (nv >= {})".format(len(new_linetracks), th_num_support_images))
    return new_linetracks
