import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Hypersim import read_raydepth, raydepth2depth

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import limap.base as _base


class HypersimDepthReader(_base.BaseDepthReader):
    def __init__(self, filename, K, img_hw):
        super(HypersimDepthReader, self).__init__(filename)
        self.K = K
        self.img_hw = img_hw

    def read(self, filename):
        raydepth = read_raydepth(filename, resize_hw=self.img_hw)
        depth = raydepth2depth(raydepth, self.K, self.img_hw)
        return depth


def read_scene_hypersim(cfg, dataset, scene_id, cam_id=0, load_depth=False, image_ids=None):
    # set scene id
    dataset.set_scene_id(scene_id)
    dataset.set_max_dim(cfg["max_image_dim"])

    # generate image indexes
    if image_ids is not None:
        index_list = image_ids
    else:
        index_list = np.arange(0, cfg["input_n_views"], cfg["input_stride"]).tolist()
        index_list = dataset.filter_index_list(index_list, cam_id=cam_id)

    # get image collections
    K = dataset.K.astype(np.float32)
    img_hw = [dataset.h, dataset.w]
    Ts, Rs = dataset.load_cameras(cam_id=cam_id)
    cameras, camimages = {}, {}
    cameras[0] = _base.Camera("SIMPLE_PINHOLE", K, cam_id=0, hw=img_hw)
    for image_id in index_list:
        pose = _base.CameraPose(Rs[image_id], Ts[image_id])
        imname = dataset.load_imname(image_id, cam_id=cam_id)
        camimage = _base.CameraImage(0, pose, image_name=imname)
        camimages[image_id] = camimage
    imagecols = _base.ImageCollection(cameras, camimages)

    if load_depth:
        # get depths
        depths = {}
        for image_id in index_list:
            depth_fname = dataset.load_raydepth_fname(image_id, cam_id=cam_id)
            depth = HypersimDepthReader(depth_fname, K, img_hw)
            depths[image_id] = depth
        return imagecols, depths
    else:
        return imagecols
