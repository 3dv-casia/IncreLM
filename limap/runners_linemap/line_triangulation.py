import os
from tqdm import tqdm

import time

import limap.base as _base
import limap.vplib as _vplib
import limap.pointsfm as _psfm
import limap.runners as _runners
import limap.util.io as limapio
import limap.visualize as limapvis

import _limap._linemap as _linemap
import limap.runners_linemap as _runners_linemap
import limap.optimize as _optim
import limap.structures as _structures


def line_triangulation(cfg, imagecols, neighbors=None, ranges=None, colmap_model_path=None):
    '''
    [IncreLM] Main interface of line triangulation over multi-view images.

    Args:
        cfg (dict): Configuration. Fields refer to :file:`cfgs_linemap/triangulation/default.yaml` as an example
        imagecols (:class:`limap.base.ImageCollection`): The image collection corresponding to all the images of interest
        neighbors (dict[int -> list[int]], optional): visual neighbors for each image. By default we compute neighbor information from the covisibility of COLMAP triangulation.
        ranges (pair of :class:`np.array` each of shape (3,), optional): robust 3D ranges for the scene. By default we compute range information from the COLMAP triangulation.
        colmap_model_path: undistorted colmap reconstruction path
    Returns:
        linetracks (list[:class:`limap.base.LineTrack`]): list of output 3D line tracks
        colmap_model_path: the colmap reconstruction path offering point sfm information
        linetracks_nv (list[:class:`limap.base.LineTrack`]): list of output 3D line tracks filtered by cfg["n_visible_views"]
    '''
    time_dict = dict()

    print("[LOG] Number of images: {0}".format(imagecols.NumImages()))
    cfg = _runners.setup(cfg)
    detector_name = cfg["line2d"]["detector"]["method"]
    if cfg["triangulation"]["var2d"] == -1:
        cfg["triangulation"]["var2d"] = cfg["var2d"][detector_name]

    if not imagecols.IsUndistorted():
        raise RuntimeError("error: input imagecols is undistorted.")

    # resize cameras
    assert imagecols.IsUndistorted() == True
    if cfg["max_image_dim"] != -1 and cfg["max_image_dim"] is not None:
        imagecols.set_max_image_dim(cfg["max_image_dim"])
    limapio.save_txt_imname_dict(os.path.join(cfg["dir_save"], 'image_list.txt'), imagecols.get_image_name_dict())
    limapio.save_npy(os.path.join(cfg["dir_save"], 'imagecols.npy'), imagecols.as_dict())

    ##########################################################
    # [1] sfm metainfos (neighbors, ranges)
    ##########################################################

    t1 = time.time()
    sfminfos_colmap_folder = None
    if neighbors is None:
        sfminfos_colmap_folder, neighbors, ranges = _runners.compute_sfminfos(cfg, imagecols)
    else:
        limapio.save_txt_metainfos(os.path.join(cfg["dir_save"], "metainfos.txt"), neighbors, ranges)
        neighbors = imagecols.update_neighbors(neighbors)
        for img_id, neighbor in neighbors.items():
            neighbors[img_id] = neighbors[img_id][:cfg["n_neighbors"]]
    limapio.save_txt_metainfos(os.path.join(cfg["dir_save"], "metainfos.txt"), neighbors, ranges)
    t2 = time.time()
    time_dict["point sfm"] = t2 - t1

    ##########################################################
    # [2] get 2D line segments for each image
    ##########################################################

    t1 = time.time()
    compute_descinfo = ((not cfg["load_match"]) and (not cfg["load_det"])) or cfg["line2d"]["compute_descinfo"]
    all_2d_segs, descinfo_folder = _runners.compute_2d_segs(cfg, imagecols, compute_descinfo=compute_descinfo)
    t2 = time.time()
    time_dict["detect line"] = t2 - t1

    # convert np.array to limap::Line2d
    all_2d_lines = _base.get_all_lines_2d(all_2d_segs)

    ##########################################################
    # [3] get line matches
    ##########################################################

    if not cfg["triangulation"]["match_lines_by_epipolar_IoU"]:
        t1 = time.time()
        print("descinfo_folder =", descinfo_folder)
        matches_dir = _runners.compute_matches(cfg, descinfo_folder, imagecols.get_img_ids(), neighbors)
        t2 = time.time()
        time_dict["match line"] = t2 - t1

    ##########################################################
    # [4] build 2D point-line association
    ##########################################################

    t1 = time.time()
    if cfg["triangulation"]["use_pointsfm"] or cfg["global_pl_association"]['use_pointsfm']:
        print("input colmap_model_path =", colmap_model_path)
        if colmap_model_path is None:
            if cfg["triangulation"]["pointsfm_config"]["colmap_folder"] is None:
                # check if colmap model exists from sfminfos computation
                if cfg["triangulation"]["pointsfm_config"]["reuse_sfminfos_colmap"] and sfminfos_colmap_folder is not None:
                    colmap_model_path = os.path.join(sfminfos_colmap_folder, "sparse")
                    if not _psfm.check_exists_colmap_model(colmap_model_path):
                        colmap_model_path = None
                # retriangulate
                if colmap_model_path is None:
                    colmap_output_path = os.path.join(cfg["dir_save"], "colmap_outputs_junctions")
                    input_neighbors = None
                    if cfg["triangulation"]["pointsfm_config"]["use_neighbors"]:
                        input_neighbors = neighbors
                    _psfm.run_colmap_sfm_with_known_poses(
                        cfg["sfm"], imagecols, output_path=colmap_output_path, skip_exists=cfg["skip_exists"], neighbors=input_neighbors)
                    colmap_model_path = os.path.join(colmap_output_path, "sparse")
            else:
                colmap_model_path = cfg["triangulation"]["pointsfm_config"]["colmap_folder"]
        print("used colmap_model_path = ", colmap_model_path)
        reconstruction = _psfm.PyReadCOLMAP(colmap_model_path)
        pointtracks = _psfm.ReadPointTracks(reconstruction, imagecols)
        all_bpt2ds, sfm_points = _runners_linemap.compute_2d_bipartites_from_colmap(
            reconstruction, imagecols, all_2d_lines, cfg["structures"]["bpt2d"], min_support_images=3)
    t2 = time.time()
    time_dict["build 2D point-line association"] = t2 - t1

    ##########################################################
    # [5] detect vanishing points (i.e., build 2D line-VP association)
    ##########################################################

    t1 = time.time()
    if cfg["triangulation"]["use_vp"] or cfg["global_pl_association"]["use_vp"]:
        if not cfg["load_vpdet"]:
            print("Start detecting VPs...", flush=True)
            vpdetector = _vplib.get_vp_detector(cfg["triangulation"]["vpdet_config"],
                                                n_jobs=cfg["triangulation"]["vpdet_config"]["n_jobs"])
            vpresults = vpdetector.detect_vp_all_images(all_2d_lines, imagecols.get_map_camviews())
            vp_save_dir = os.path.join(cfg["dir_save"], "vp_results", cfg["line2d"]["detector"]["method"])
            limapio.check_makedirs(vp_save_dir)
            limapio.write_pkl(vpresults, os.path.join(vp_save_dir, "vp_results.pkl"))
        else:
            print("Start loading VPs...", flush=True)
            vp_load_dir = os.path.join(cfg["dir_load"], "vp_results", cfg["line2d"]["detector"]["method"])
            vpresults = limapio.read_pkl(os.path.join(vp_load_dir, "vp_results.pkl"))
    t2 = time.time()
    time_dict["build 2D line-VP association (i.e. detect VP)"] = t2 - t1

    if cfg["triangulation"]["only_execute_pre_processing"]:
        print("Exit the program, because only the pre-processing is required.")
        return None, None

    ##########################################################
    # [6] initialize line mapper
    ##########################################################

    t1 = time.time()
    line_mapper = _linemap.LineMapper(cfg["triangulation"], imagecols, all_2d_lines, ranges)
    if cfg["triangulation"]["use_pointsfm"]:
        line_mapper.SetPLBipartite2d(all_bpt2ds)
        line_mapper.SetSfMPoints(sfm_points)
        line_mapper.SetPointTracks(pointtracks)
    if cfg["triangulation"]["use_vp"]:
        line_mapper.SetVPResults(vpresults)
    line_mapper.Initialize()
    t2 = time.time()
    time_dict["initialize line mapper"] = t2 - t1

    ##########################################################
    # [7] load or match line matches
    ##########################################################

    if not cfg["triangulation"]["match_lines_by_epipolar_IoU"]:
        t1 = time.time()
        print('Start loading line matches...', flush=True)
        for img_id in tqdm(imagecols.get_img_ids()):
            matches = limapio.read_npy(os.path.join(matches_dir, "matches_{0}.npy".format(img_id))).item()
            line_mapper.LoadMatches(img_id, matches)
        t2 = time.time()
        time_dict["load line matches"] = t2 - t1
    else:
        t1 = time.time()
        print('Start matching lines by epipolar IoU...', flush=True)
        line_mapper.MatchLines2dByEpipolarIoU(
            neighbors, cfg["triangulation"]["max_num_one_to_many_matches"], cfg["triangulation"]["th_IoU"])
        t2 = time.time()
        time_dict["match lines"] = t2 - t1
    line_mapper.CountMatches()

    ##########################################################
    # [8] construct tri nodes
    ##########################################################

    t_mapping_begin = time.time()

    # config hybrid ransac
    hybrid_ransac_options = _linemap.HybridLineTriangulationEstimatorOptions()
    ransac_cfg = cfg["triangulation"]["ransac"]
    hybrid_ransac_options.ransac_options.min_num_iterations_ = ransac_cfg["min_num_iterations"]
    hybrid_ransac_options.ransac_options.squared_inlier_thresholds_ = [ransac_cfg["inlier_si_dist_th"] * ransac_cfg["inlier_si_dist_th"],
                                                                       ransac_cfg["inlier_angle_th"] * ransac_cfg["inlier_angle_th"]]
    hybrid_ransac_options.ransac_options.data_type_weights_ = [
        ransac_cfg["weight_point"], ransac_cfg["weight_direction"]]
    hybrid_ransac_options.ransac_options.max_num_iterations_ = 200
    hybrid_ransac_options.cheirality_min_depth = ransac_cfg["cheirality_min_depth"]
    solver_flags = []
    if ransac_cfg["use_two_points_solver"]:
        solver_flags.append(True)
    else:
        solver_flags.append(False)
    if ransac_cfg["use_point_vp_solver"]:
        solver_flags.append(True)
    else:
        solver_flags.append(False)
    hybrid_ransac_options.solver_flags = solver_flags
    hybrid_ransac_options.random = ransac_cfg["random"]
    hybrid_ransac_options.th_angle = ransac_cfg["th_angle"]
    hybrid_ransac_options.th_perp = ransac_cfg["th_perp"]
    hybrid_ransac_options.Print()  # print options

    line_mapper.SetHybridRansacOptions(hybrid_ransac_options)

    t1 = time.time()
    print('Start constructing tri nodes (i.e. triangulating 3D line segments from two views)...', flush=True)
    line_mapper.ConstructTriNodes()
    t2 = time.time()
    time_dict["construct tri nodes"] = t2 - t1

    # statistics
    # all_hypotheses = line_mapper.GetAllHypotheses()
    # linetracks_folder = os.path.join(cfg["dir_save"], "all_hypotheses", cfg["output_folder"])
    # limapio.save_folder_linetracks_with_info(
    #     linetracks_folder, all_hypotheses, config=cfg, imagecols=imagecols, all_2d_segs=all_2d_segs)
    # assert False

    ##########################################################
    # [9] construct tri edges
    ##########################################################

    t1 = time.time()
    print('Start constructing tri edges...', flush=True)
    line_mapper.ConstructTriEdges()
    t2 = time.time()
    time_dict["construct tri edges"] = t2 - t1

    ##########################################################
    # [10] incrementally reconstruct 3D LS one-by-one with known camera poses
    ##########################################################

    print('Start incremental reconstruction with known camera poses...', flush=True)

    time_dict["find best tracks"] = 0
    time_dict["extend tracks"] = 0
    time_dict["merge tracks"] = 0
    time_dict["optimization"] = 0
    time_dict["filter tracks"] = 0

    n_iter = 0
    while True:
        print("=" * 78)
        print("incremental iter =", n_iter)
        print("=" * 78, flush=True)

        t1 = time.time()
        print('Start finding the best line track...', flush=True)
        if cfg["triangulation"]["use_fast_find_max_version"]:
            new_track_id = line_mapper.FindAndAddBestLineTrackFastFindMax()
        else:
            new_track_id = line_mapper.FindAndAddBestLineTrack()
        t2 = time.time()
        time_dict["find best tracks"] += t2 - t1

        if line_mapper.Stop():
            if cfg["triangulation"]["output_submaps"]:
                linetracks = line_mapper.GetLineTrackList()  # has `n_iter` line tracks
                linetracks_folder = os.path.join(cfg["dir_save"], "submaps", str(n_iter))
                limapio.save_folder_linetracks_with_info(
                    linetracks_folder, linetracks, config=None, imagecols=None, all_2d_segs=None)
            break

        t1 = time.time()
        print('Start extending the new line track...', flush=True)
        line_mapper.ExtendTrack(new_track_id)
        t2 = time.time()
        time_dict["extend tracks"] += t2 - t1

        if cfg["triangulation"]["output_submaps"] and ((n_iter + 1) in [100, 300, 500]):
            linetracks = line_mapper.GetLineTrackList()  # has (n_iter + 1) line tracks
            linetracks_folder = os.path.join(cfg["dir_save"], "submaps", str(n_iter))
            limapio.save_folder_linetracks_with_info(
                linetracks_folder, linetracks, config=None, imagecols=None, all_2d_segs=None)

        n_iter += 1

    ##########################################################
    # [11] post-processing
    ##########################################################

    linetrack_before_merging = line_mapper.GetLineTrackList()
    linetracks_folder = os.path.join(cfg["dir_save"], "before_merging", cfg["output_folder"])
    limapio.save_folder_linetracks_with_info(
        linetracks_folder, linetrack_before_merging, config=cfg, imagecols=imagecols, all_2d_segs=all_2d_segs)

    t1 = time.time()
    if cfg["triangulation"]["enable_merging"]:
        while True:
            num_merged_tracks = line_mapper.MergeLineTracks()
            if num_merged_tracks == 0:
                break
    t2 = time.time()
    time_dict["merge tracks"] += t2 - t1

    linetrack_after_merging = line_mapper.GetLineTrackList()
    linetracks_folder = os.path.join(cfg["dir_save"], "after_merging", cfg["output_folder"])
    limapio.save_folder_linetracks_with_info(
        linetracks_folder, linetrack_after_merging, config=cfg, imagecols=imagecols, all_2d_segs=all_2d_segs)

    t1 = time.time()
    if cfg["global_pl_association"]["enable"]:
        # the following joint optimization program is inspired by LIMAP (https://github.com/cvg/limap/blob/main/runners/pointline_association.py)
        line_mapper.FilterLineTracksByNumImages(3)
        linetracks_map = line_mapper.GetLineTrackMap()
        linetracks_list = line_mapper.GetLineTrackList()
        # optimize association
        cfg_associator = _optim.GlobalAssociatorConfig(cfg["global_pl_association"])
        cfg_associator.solver_options.max_num_iterations = cfg["global_pl_association"]["max_num_iterations"]
        # cfg_associator.solver_options.logging_type = _ceresbase.LoggingType.STDOUT
        associator = _optim.GlobalAssociator(cfg_associator)
        associator.InitImagecols(imagecols)
        associator.InitLineTracks(linetracks_map)
        if cfg["global_pl_association"]['use_pointsfm']:
            associator.InitPointTracks(pointtracks)
            associator.Init2DBipartites_PointLine(all_bpt2ds)
        # associator.ReassociateJunctions()
        if cfg["global_pl_association"]["use_vp"]:
            vptrack_constructor = _vplib.GlobalVPTrackConstructor()
            vptrack_constructor.Init(vpresults)
            vptracks = vptrack_constructor.ClusterLineTracks(linetracks_list, imagecols)
            all_bpt2ds_vp = _structures.GetAllBipartites_VPLine2d(all_2d_lines, vpresults, vptracks)
            associator.InitVPTracks(vptracks)
            associator.Init2DBipartites_VPLine(all_bpt2ds_vp)
        associator.SetUp()
        associator.Solve()
        linetracks_map = associator.GetOutputLineTracks_linemap(
            cfg["global_pl_association"]["num_outliers_aggregator"], cfg["triangulation"]["var2d"])  # std::map format
        line_mapper.SetLineTracks(linetracks_map)
    t2 = time.time()
    time_dict["optimization"] += t2 - t1

    t1 = time.time()
    if cfg["triangulation"]["enable_filtering"]:
        line_mapper.FilterLineTracksByReprojection()
        line_mapper.FilterLineTracksBySensitivity()
        line_mapper.FilterLineTracksByOverlap()
    t2 = time.time()
    time_dict["filter tracks"] += t2 - t1

    t_mapping_end = time.time()
    time_dict["total 3D line mapping time"] = t_mapping_end - t_mapping_begin

    ##########################################################
    # [12] output & visualize
    ##########################################################

    linetracks = line_mapper.GetLineTrackList()
    line_mapper.FilterLineTracksByNumImages(cfg["n_visible_views"])
    linetracks_nv = line_mapper.GetLineTrackList()

    # save filtered line tracks according to cfg["n_visible_views"]
    limapio.save_txt_linetracks(os.path.join(cfg["dir_save"], "alltracks_nv{0}.txt".format(
        cfg["n_visible_views"])), linetracks, n_visible_views=cfg["n_visible_views"])

    # save folder-to-linetracks
    # Note: none of the line tracks contained in linetracks_folder are filtered by n_visible_views
    linetracks_folder = os.path.join(cfg["dir_save"], cfg["output_folder"])
    limapio.save_folder_linetracks_with_info(
        linetracks_folder, linetracks, config=cfg, imagecols=imagecols, all_2d_segs=all_2d_segs)

    VisTrack = limapvis.Open3DTrackVisualizer(linetracks)
    VisTrack.report()
    limapio.save_obj(os.path.join(cfg["dir_save"], 'triangulated_lines_nv{0}.obj'.format(
        cfg["n_visible_views"])), VisTrack.get_lines_np(n_visible_views=cfg["n_visible_views"]))

    # visualize
    if cfg["visualize"]:
        validtracks = [track for track in linetracks if track.count_images() >= cfg["n_visible_views"]]

        def report_track(track_id):
            limapvis.visualize_line_track(
                imagecols, validtracks[track_id], prefix="track.{0}".format(track_id))
        import pdb
        pdb.set_trace()
        VisTrack.vis_reconstruction(
            imagecols, n_visible_views=cfg["n_visible_views"], width=2)
        pdb.set_trace()

    print("Report running times:")
    for k, v in time_dict.items():
        print(f"{k} : {v:.4f} seconds")

    return linetracks, colmap_model_path, linetracks_nv
