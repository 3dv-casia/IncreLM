import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Hypersim import Hypersim
from loader import read_scene_hypersim

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import limap.util.config as cfgutils
import limap.util.evaluation as limapeval
import limap.util.io as limapio
import limap.runners
import limap.optimize
import limap.pointsfm

from runners_linemap.colmap_triangulation import run_colmap_triangulation
from limap.runners_linemap import filter_tracks_by_num_support_images


def run_scene_hypersim(imagecols_output_dir, tri_cfg, hypersim_dataset, scene_id, cam_id=0):
    imagecols_gt = read_scene_hypersim(tri_cfg, hypersim_dataset, scene_id, cam_id=cam_id, load_depth=False)

    ###########################################################################
    # [1] run colmap to obtain point tracks
    ###########################################################################

    import limap.pointsfm as _psfm
    tri_cfg = limap.runners.setup(tri_cfg)

    # global_dir_save = tri_cfg["dir_save"]
    limapio.check_makedirs(imagecols_output_dir)
    limapio.save_npy(os.path.join(imagecols_output_dir, "imagecols_gt.npy"), imagecols_gt)

    ###################################################
    # If the following code had been executed, the colmap results will be saved in `colmap_path`.
    # You can comment the following code and set `colmap_path` manually to reuse colmap results.
    colmap_path = os.path.join(tri_cfg["dir_save"], "colmap_sfm")
    _psfm.run_colmap_sfm(tri_cfg["sfm"], imagecols_gt, output_path=colmap_path,
                         skip_exists=tri_cfg["skip_exists"], map_to_original_image_names=False)
    ###################################################

    imagecols, _, _ = _psfm.read_infos_colmap(tri_cfg["sfm"], colmap_path, model_path="sparse/0", image_path="images")
    limapio.save_npy(os.path.join(imagecols_output_dir, "imagecols_sfm.npy"), imagecols)

    ###########################################################################
    # [2] given colmap sfm poses, run IncreLM to obtain line tracks
    ###########################################################################

    colmap_folder = os.path.join(colmap_path, "sparse/0")
    tri_cfg["info_path"] = None
    _, _, linetracks = run_colmap_triangulation(tri_cfg, colmap_path, model_path="sparse/0", image_path="images")
    linetracks = filter_tracks_by_num_support_images(linetracks, 4)

    ###########################################################################
    # [3] run point-line joint BA to optimize poses
    ###########################################################################

    reconstruction = limap.pointsfm.PyReadCOLMAP(colmap_folder)
    pointtracks = limap.pointsfm.ReadPointTracks(reconstruction, imagecols)
    cfg_ba = limap.optimize.HybridBAConfig()
    ba_engine = limap.optimize.solve_hybrid_bundle_adjustment(cfg_ba, imagecols, pointtracks, linetracks)
    new_imagecols = ba_engine.GetOutputImagecols()

    # save optimized poses
    limapio.save_npy(os.path.join(imagecols_output_dir, "imagecols_optimized.npy"), new_imagecols)

    ###########################################################################
    # [4] evaluate absolute and relative pose error using GT poses
    ###########################################################################

    print("=" * 78)
    print("Evaluate absolute and relative pose error using GT poses")
    print("=" * 78)
    print("eval scene_id =", scene_id)
    trans_errs_orig, rot_errs_orig, acc_orig = limapeval.eval_imagecols_with_acc_percentage(imagecols, imagecols_gt)
    trans_errs_new, rot_errs_new, acc_new = limapeval.eval_imagecols_with_acc_percentage(new_imagecols, imagecols_gt)
    print(f"[median] rot_errs_orig = {np.median(rot_errs_orig)}, trans_errs_orig = {np.median(trans_errs_orig)}")
    print(f"[median] rot_errs_new = {np.median(rot_errs_new)}, trans_errs_new = {np.median(trans_errs_new)}")


def parse_config():
    import argparse
    arg_parser = argparse.ArgumentParser(
        description='refining COLMAP SfM poses using line tracks from IncreLM')
    arg_parser.add_argument('-c', '--tri_config_file', type=str,
                            default='cfgs_linemap/triangulation/hypersim.yaml', help='triangulation config file')
    arg_parser.add_argument('--tri_default_config_file', type=str,
                            default='cfgs_linemap/triangulation/default.yaml', help='triangulation default config file')
    arg_parser.add_argument('--npyfolder', type=str, default=None,
                            help='folder to load precomputed results')
    arg_parser.add_argument('--tri_output_dir', type=str,
                            default="tmp/tri", help='triangulation output dir (unoptimized line map)')
    arg_parser.add_argument('--imagecols_output_dir', type=str,
                            default="tmp/imagecols", help='imagecols output dir')

    args, unknown = arg_parser.parse_known_args()
    tri_cfg = cfgutils.load_config(
        args.tri_config_file, default_path=args.tri_default_config_file)
    shortcuts = dict()
    shortcuts['-nv'] = '--n_visible_views'
    shortcuts['-sid'] = '--scene_id'
    tri_cfg = cfgutils.update_config(tri_cfg, unknown, shortcuts)
    tri_cfg["folder_to_load"] = args.npyfolder
    if tri_cfg["folder_to_load"] is None:
        tri_cfg["folder_to_load"] = os.path.join("precomputed", "hypersim", tri_cfg["scene_id"])
    tri_cfg["output_dir"] = args.tri_output_dir
    print("args.imagecols_output_dir =", args.imagecols_output_dir)
    print("tri_cfg =", tri_cfg)
    return args.imagecols_output_dir, tri_cfg


def main():
    imagecols_output_dir, tri_cfg = parse_config()
    dataset = Hypersim(tri_cfg["data_dir"])
    run_scene_hypersim(imagecols_output_dir, tri_cfg, dataset,
                       tri_cfg["scene_id"], cam_id=tri_cfg["cam_id"])


if __name__ == '__main__':
    main()
