import os
import cv2
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import copy
import math
import random


def random_color():
    r = int(255 * np.random.rand())
    g = int(255 * np.random.rand())
    b = 255 - r
    return b, g, r


def draw_points(image, points, color=None, thickness=1):
    image = copy.deepcopy(image)
    for p in points:
        c = random_color() if color is None else color
        pos_x, pos_y = int(round(p[0])), int(round(p[1]))
        cv2.circle(image, (pos_x, pos_y), thickness * 2, c, -1)
    return image


def draw_segments(image, segments, color=None, thickness=1, endpoints=True):
    image = copy.deepcopy(image)
    for s in segments:
        c = random_color() if color is None else color
        p1 = (int(s[0]), int(s[1]))
        p2 = (int(s[2]), int(s[3]))
        cv2.line(image, p1, p2, c, thickness)
        if endpoints:
            cv2.circle(image, p1, thickness * 2, c, -1)
            cv2.circle(image, p2, thickness * 2, c, -1)
    return image


def draw_salient_segments(image, segments, saliency, color1=(0, 255, 0), color2=(0, 0, 255), thickness=1,
                          endpoints=True):
    assert len(segments) == len(saliency)
    image = copy.deepcopy(image)
    max_saliency, min_saliency = np.max(saliency), np.min(saliency)
    for s, s_saliency in zip(segments, saliency):
        r = (s_saliency - min_saliency) / (max_saliency - min_saliency)
        c = (int(color1[0] * r + color2[0] * (1 - r)),
             int(color1[1] * r + color2[1] * (1 - r)),
             int(color1[2] * r + color2[2] * (1 - r)))
        p1 = (int(s[0]), int(s[1]))
        p2 = (int(s[2]), int(s[3]))
        cv2.line(image, p1, p2, c, thickness)
        if endpoints:
            cv2.circle(image, p1, thickness * 2, c, -1)
            cv2.circle(image, p2, thickness * 2, c, -1)
    return image


def draw_multiscale_segments(img, segs, color=None, endpoints=True, thickness=2):
    img = copy.deepcopy(img)
    for s in segs:
        mycolor = random_color() if color is None else color
        octave, l = s[0]
        l = l * np.sqrt(2) ** octave
        cv2.line(img, (int(l[0]), int(l[1])), (int(l[2]), int(l[3])), mycolor, thickness)
        if endpoints:
            cv2.circle(img, (int(l[0]), int(l[1])), int(thickness * 1.5), mycolor, -1)
            cv2.circle(img, (int(l[2]), int(l[3])), int(thickness * 1.5), mycolor, -1)
    return img


def draw_multiscale_matches(img_left, img_right, segs_left, segs_right, matches):
    assert img_left.ndim == 2
    h, w = img_left.shape

    # store the matching results of the first and second images into a single image
    left_color_img = cv2.cvtColor(img_left, cv2.COLOR_GRAY2BGR)
    right_color_img = cv2.cvtColor(img_right, cv2.COLOR_GRAY2BGR)
    r1, g1, b1 = [], [], []  # the line colors

    for pair in range(len(matches)):
        r, g, b = random_color()
        r1.append(r)
        g1.append(g)
        b1.append(b)
        line_id_l, line_id_r = matches[pair]

        octave, l = segs_left[line_id_l][0]
        l = l * np.sqrt(2) ** octave
        cv2.line(left_color_img, (int(l[0]), int(l[1])), (int(l[2]), int(l[3])), (r1[pair], g1[pair], b1[pair]), 3)

        octave, l = segs_right[line_id_r][0]
        l = l * np.sqrt(2) ** octave
        cv2.line(right_color_img, (int(l[0]), int(l[1])), (int(l[2]), int(l[3])), (r1[pair], g1[pair], b1[pair]), 3)

    result_img = np.hstack([left_color_img, right_color_img])
    result_img_tmp = result_img.copy()
    for pair in range(27, len(matches)):
        line_id_l, line_id_r = matches[pair]
        octave_left, seg_left = segs_left[line_id_l][0]
        octave_right, seg_right = segs_right[line_id_r][0]
        seg_left = seg_left[:4] * (np.sqrt(2) ** octave_left)
        seg_right = seg_right[:4] * (np.sqrt(2) ** octave_right)

        start_ptn = (int(seg_left[0]), int(seg_left[1]))
        end_ptn = (int(seg_right[0] + w), int(seg_right[1]))
        cv2.line(result_img_tmp, start_ptn, end_ptn, (r1[pair], g1[pair], b1[pair]), 2, cv2.LINE_AA)

    result_img = cv2.addWeighted(result_img, 0.5, result_img_tmp, 0.5, 0.0)
    return result_img


def draw_singlescale_matches(img_left, img_right, matched_segs, mask=None):
    assert img_left.ndim == 2
    h, w = img_left.shape

    # store the matching results of the first and second images into a single image
    left_color_img = cv2.cvtColor(img_left, cv2.COLOR_GRAY2BGR)
    right_color_img = cv2.cvtColor(img_right, cv2.COLOR_GRAY2BGR)
    colors = []
    if mask is None:
        for l, r in matched_segs:
            color = random_color()
            colors.append(color)
            cv2.line(left_color_img, (int(l[0]), int(l[1])), (int(l[2]), int(l[3])), color, 3)
            cv2.line(right_color_img, (int(r[0]), int(r[1])), (int(r[2]), int(r[3])), color, 3)
    else:
        for correct, (l, r) in zip(mask, matched_segs):
            color = (0, 255, 0) if correct else (0, 0, 255)
            cv2.line(left_color_img, (int(l[0]), int(l[1])), (int(l[2]), int(l[3])), color, 3)
            cv2.line(right_color_img, (int(r[0]), int(r[1])), (int(r[2]), int(r[3])), color, 3)

    result_img = np.hstack([left_color_img, right_color_img])
    result_img_tmp = result_img.copy()
    if mask is None:
        for i in range(len(matched_segs)):
            l, r = matched_segs[i]
            start_ptn = (int(l[0]), int(l[1]))
            end_ptn = (int(r[0] + w), int(r[1]))
            cv2.line(result_img_tmp, start_ptn, end_ptn, colors[i], 2, cv2.LINE_AA)
    else:
        for correct, (l, r) in zip(mask, matched_segs):
            color = (0, 255, 0) if correct else (0, 0, 255)
            start_ptn = (int(l[0]), int(l[1]))
            end_ptn = (int(r[0] + w), int(r[1]))
            cv2.line(result_img_tmp, start_ptn, end_ptn, color, 2, cv2.LINE_AA)

    result_img = cv2.addWeighted(result_img, 0.5, result_img_tmp, 0.5, 0.0)
    return result_img


def crop_to_patch(img, center, patch_size=50):
    img_h, img_w = img.shape[0], img.shape[1]
    if img_h < patch_size or img_w < patch_size:
        raise ValueError("Error! Image size should be larger than the given patch size!")

    # get range
    cx, cy = center[0], center[1]
    start_x = int(cx - patch_size / 2.0)
    if start_x < 0:
        start_x = 0
    if start_x + patch_size > img_w:
        start_x = img_w - patch_size
    end_x = start_x + patch_size
    start_y = int(cy - patch_size / 2.0)
    if start_y < 0:
        start_y = 0
    if start_y + patch_size > img_h:
        start_y = img_h - patch_size
    end_y = start_y + patch_size

    # crop
    if len(img.shape) == 2:
        newimg = img[start_y:end_y, start_x:end_x]
        return newimg
    elif len(img.shape) == 3:
        newimg = img[start_y:end_y, start_x:end_x, :]
        return newimg
    else:
        raise NotImplementedError


def cat_to_bigimage(imgs, shape, pad=20):
    '''
    make a big image with 2d image collections
    all images should have the same size
    '''
    img = imgs[0][0]
    img_h, img_w = img.shape[0], img.shape[1]
    channel = img.shape[2]
    n_rows, n_cols = shape[0], shape[1]

    final_img_h = img_h * n_rows + pad * (n_rows - 1)
    final_img_w = img_w * n_cols + pad * (n_cols - 1)
    bigimg = (np.ones((final_img_h, final_img_w, channel)) * 255).astype(np.uint8)
    for row in range(n_rows):
        for col in range(n_cols):
            row_start = (img_h + pad) * row
            row_end = row_start + img_h
            col_start = (img_w + pad) * col
            col_end = col_start + img_w
            img = imgs[row][col]
            bigimg[row_start:row_end, col_start:col_end, :] = img
    return bigimg


def make_bigimage(imgs, pad=20):
    n_images = len(imgs)
    if n_images <= 5:
        return cat_to_bigimage([imgs], (1, n_images), pad=pad)

    n_cols = int(np.ceil(np.sqrt(n_images)))
    n_rows = int(np.ceil(n_images / n_cols))
    imgs_collections = []
    imgs_now = []
    for image in imgs:
        imgs_now.append(image)
        if len(imgs_now) == n_cols:
            imgs_collections.append(imgs_now)
            imgs_now = []
    blank_image = (np.ones((imgs[0].shape[0], imgs[0].shape[1], imgs[0].shape[2])) * 255).astype(np.uint8)
    while len(imgs_now) < n_cols:
        imgs_now.append(blank_image)
    imgs_collections.append(imgs_now)
    return cat_to_bigimage(imgs_collections, (n_rows, n_cols), pad=pad)


def test_point_inside_ranges(point, ranges):
    point = np.array(point)
    if ~np.all(point > ranges[0]) or ~np.all(point < ranges[1]):
        return False
    return True


def test_line_inside_ranges(line, ranges):
    if not test_point_inside_ranges(line.start, ranges):
        return False
    if not test_point_inside_ranges(line.end, ranges):
        return False
    return True


def compute_robust_range(arr, range_robust=[0.05, 0.95], k_stretch=2.0):
    N = arr.shape[0]
    start_idx = int(round((N - 1) * range_robust[0]))
    end_idx = int(round((N - 1) * range_robust[1]))
    arr_sorted = np.sort(arr)
    start = arr_sorted[start_idx]
    end = arr_sorted[end_idx]
    start_stretched = (start + end) / 2.0 - k_stretch * (end - start) / 2.0
    end_stretched = (start + end) / 2.0 + k_stretch * (end - start) / 2.0
    return start_stretched, end_stretched


def compute_robust_range_lines(lines, range_robust=[0.05, 0.95], k_stretch=2.0):
    lines_array = np.array([line.as_array() for line in lines])
    x_array = lines_array.reshape(-1, 3)[:, 0]
    y_array = lines_array.reshape(-1, 3)[:, 1]
    z_array = lines_array.reshape(-1, 3)[:, 2]

    x_start, x_end = compute_robust_range(x_array, range_robust=range_robust, k_stretch=k_stretch)
    y_start, y_end = compute_robust_range(y_array, range_robust=range_robust, k_stretch=k_stretch)
    z_start, z_end = compute_robust_range(z_array, range_robust=range_robust, k_stretch=k_stretch)
    ranges = np.array([[x_start, y_start, z_start], [x_end, y_end, z_end]])
    return ranges


def filter_ranges(lines_np, counts_np, ranges):
    new_lines_np, new_counts_np = [], []
    for idx in range(lines_np.shape[0]):
        if test_line_inside_ranges(lines_np[idx], ranges):
            new_lines_np.append(lines_np[idx])
            new_counts_np.append(counts_np[idx])
    return np.array(new_lines_np), np.array(new_counts_np)


def report_dist_reprojection(line3d, line2d, camview, prefix=None):
    import limap.base as _base
    line2d_proj = line3d.projection(camview)
    angle = _base.compute_distance_2d(line2d, line2d_proj, _base.LineDistType.ANGULAR)
    perp_dist = _base.compute_distance_2d(line2d, line2d_proj, _base.LineDistType.PERPENDICULAR_ONEWAY)
    overlap = _base.compute_distance_2d(line2d_proj, line2d, _base.LineDistType.OVERLAP)
    sensitivity = line3d.sensitivity(camview)
    if prefix is None:
        print("angle = {0:.4f}, perp = {1:.4f}, overlap = {2:.4f}, sensi = {3:.4f}".format(
            angle, perp_dist, overlap, sensitivity))
    else:
        print("{4}: angle = {0:.4f}, perp = {1:.4f}, overlap = {2:.4f}, sensi = {3:.4f}".format(
            angle, perp_dist, overlap, sensitivity, prefix))


def visualize_2d_line(fname, imagecols, all_lines_2d, img_id, line_id):
    img = imagecols.read_image(img_id)
    img = draw_segments(img, [all_lines_2d[img_id][line_id].as_array().reshape(-1)], (0, 255, 0))
    cv2.imwrite(fname, img)


def visualize_point_track_v2(imagecols, pointtrack, output_dir="tmp", prefix='pointtrack__id_xxx'):
    print("[VISUALIZE]: prefix: {0} num_supporting_points: {1}".format(prefix, len(pointtrack.image_id_list)))
    image_point_dict = dict()
    for idx, (img_id, point2d) in enumerate(zip(pointtrack.image_id_list, pointtrack.p2d_list)):
        if img_id not in image_point_dict:
            image_point_dict[img_id] = []
        image_point_dict[img_id].append(point2d)

    for img_id, point2d_list in image_point_dict.items():
        img = imagecols.read_image(img_id)
        if len(img.shape) == 2:
            img = img[:, :, None].repeat(3, 2)
        point2d_proj = imagecols.camview(img_id).projection(pointtrack.p)
        img = draw_points(img, [np.array(point2d_proj)], (255, 0, 0), thickness=1)
        for point2d in point2d_list:
            img = draw_points(img, [np.array(point2d)], (0, 0, 255), thickness=1)
        fname = os.path.join(output_dir, '{0}__num_{1}__image_{2}.png'.format(
            prefix, len(point2d_list), os.path.basename(imagecols.camimage(img_id).image_name())[:-4]))
        cv2.imwrite(fname, img)


def visualize_line_track_v2(imagecols, linetrack, output_dir="tmp", prefix='linetrack__id_xxx'):
    print("[VISUALIZE]: prefix: {0} line length: {1}, num_supporting_lines: {2}".format(
        prefix, linetrack.line.length(), len(linetrack.image_id_list)))
    image_line_dict = dict()
    for idx, (img_id, line2d) in enumerate(zip(linetrack.image_id_list, linetrack.line2d_list)):
        if img_id not in image_line_dict:
            image_line_dict[img_id] = []
        image_line_dict[img_id].append(line2d)

    for img_id, line2d_list in image_line_dict.items():
        img = imagecols.read_image(img_id)
        if len(img.shape) == 2:
            img = img[:, :, None].repeat(3, 2)
        line2d_proj = linetrack.line.projection(imagecols.camview(img_id))
        img = draw_segments(img, [line2d_proj.as_array().reshape(-1)], (255, 0, 0), thickness=1)
        for line2d in line2d_list:
            img = draw_segments(img, [line2d.as_array().reshape(-1)], (0, 0, 255), thickness=1)
        fname = os.path.join(output_dir, '{0}__num_{1}__image_{2}.png'.format(
            prefix, len(line2d_list), os.path.basename(imagecols.camimage(img_id).image_name())[:-4]))
        cv2.imwrite(fname, img)


def visualize_point_line_track(imagecols, pointtrack, linetrack, output_dir="tmp", prefix='pointtrack_id_xxx__linetrack__id_xxx'):
    shared_image_ids = list(set(pointtrack.image_id_list).intersection(set(linetrack.image_id_list)))
    print("len(shared_image_ids) = ", len(shared_image_ids))
    image_point_dict = dict()
    for idx, (img_id, point2d) in enumerate(zip(pointtrack.image_id_list, pointtrack.p2d_list)):
        if img_id not in image_point_dict:
            image_point_dict[img_id] = []
        image_point_dict[img_id].append(point2d)

    image_line_dict = dict()
    for idx, (img_id, line2d) in enumerate(zip(linetrack.image_id_list, linetrack.line2d_list)):
        if img_id not in image_line_dict:
            image_line_dict[img_id] = []
        image_line_dict[img_id].append(line2d)

    for shared_image_id in shared_image_ids:
        img = imagecols.read_image(shared_image_id)
        if len(img.shape) == 2:
            img = img[:, :, None].repeat(3, 2)
        point2d_proj = imagecols.camview(shared_image_id).projection(pointtrack.p)
        img = draw_points(img, [np.array(point2d_proj)], (125, 0, 0), thickness=2)

        for point2d in image_point_dict[shared_image_id]:
            img = draw_points(img, [np.array(point2d)], (0, 0, 125), thickness=2)

        line2d_proj = linetrack.line.projection(imagecols.camview(shared_image_id))
        img = draw_segments(img, [line2d_proj.as_array().reshape(-1)], (255, 0, 0), thickness=1)
        for line2d in image_line_dict[shared_image_id]:
            img = draw_segments(img, [line2d.as_array().reshape(-1)], (0, 0, 255), thickness=1)

        fname = os.path.join(output_dir, '{0}__sharedimage_numpoint2d_{1}__numline2d_{2}__image_{3}.png'.format(
            prefix, len(image_point_dict[shared_image_id]), len(image_line_dict[shared_image_id]), os.path.basename(imagecols.camimage(shared_image_id).image_name())[:-4]))
        cv2.imwrite(fname, img)

    for img_id, point2d_list in image_point_dict.items():
        if img_id in shared_image_ids:
            continue
        img = imagecols.read_image(img_id)
        if len(img.shape) == 2:
            img = img[:, :, None].repeat(3, 2)
        point2d_proj = imagecols.camview(img_id).projection(pointtrack.p)
        img = draw_points(img, [np.array(point2d_proj)], (125, 0, 0), thickness=2)
        for point2d in point2d_list:
            img = draw_points(img, [np.array(point2d)], (0, 0, 125), thickness=2)
        fname = os.path.join(output_dir, '{0}__unsharedimage__numpoint2d_{1}__image_{2}.png'.format(
            prefix, len(point2d_list), os.path.basename(imagecols.camimage(img_id).image_name())[:-4]))
        cv2.imwrite(fname, img)

    for img_id, line2d_list in image_line_dict.items():
        if img_id in shared_image_ids:
            continue
        img = imagecols.read_image(img_id)
        if len(img.shape) == 2:
            img = img[:, :, None].repeat(3, 2)
        line2d_proj = linetrack.line.projection(imagecols.camview(img_id))
        img = draw_segments(img, [line2d_proj.as_array().reshape(-1)], (255, 0, 0), thickness=1)
        for line2d in line2d_list:
            img = draw_segments(img, [line2d.as_array().reshape(-1)], (0, 0, 255), thickness=1)
        fname = os.path.join(output_dir, '{0}__unsharedimage__numline2d_{1}__image_{2}.png'.format(
            prefix, len(line2d_list), os.path.basename(imagecols.camimage(img_id).image_name())[:-4]))
        cv2.imwrite(fname, img)


def visualize_line_track(imagecols, linetrack, prefix='linetrack', report=False):
    print("[VISUALIZE]: line length: {0}, num_supporting_lines: {1}".format(
        linetrack.line.length(), len(linetrack.image_id_list)))
    for idx, (img_id, line2d) in enumerate(zip(linetrack.image_id_list, linetrack.line2d_list)):
        img = imagecols.read_image(img_id)
        if len(img.shape) == 2:
            img = img[:, :, None].repeat(3, 2)
        report_dist_reprojection(linetrack.line, line2d, imagecols.camview(
            img_id), prefix="Reprojecting to line {0} (img {1}, line {2})".format(idx, img_id, linetrack.line_id_list[idx]))
        line2d_proj = linetrack.line.projection(imagecols.camview(img_id))
        img = draw_segments(img, [line2d_proj.as_array().reshape(-1)], (255, 0, 0), thickness=1)
        img = draw_segments(img, [line2d.as_array().reshape(-1)], (0, 0, 255), thickness=1)
        fname = os.path.join('tmp_draw_track', '{0}.{1}.{2}.png'.format(
            prefix, idx, os.path.basename(imagecols.camimage(img_id).image_name())[:-4]))
        cv2.imwrite(fname, img)


def vis_vpresult(img, lines, vpres, vp_id=-1, show_original=False, endpoints=False):
    import seaborn as sns
    import cv2
    n_vps = vpres.count_vps()
    colors = sns.color_palette("husl", n_colors=n_vps)
    colors = (np.array(colors) * 255).astype(np.uint8).tolist()
    if n_vps == 1:
        colors = [[255, 0, 0]]
    for line_id, line in enumerate(lines):
        c = [255, 255, 255]  # default color: white
        if not vpres.HasVP(line_id):
            if not show_original:
                continue
        elif vp_id >= 0 and vpres.labels[line_id] != vp_id:
            if not show_original:
                continue
        else:
            c = colors[vpres.labels[line_id]]
        cv2.line(img, (int(line.start[0]), int(line.start[1])), (int(line.end[0]), int(line.end[1])), c, 2)
        if endpoints:
            cv2.circle(img, (int(line.start[0]), int(line.start[1])), 3, [0, 0, 0], -1)
            cv2.circle(img, (int(line.end[0]), int(line.end[1])), 3, [0, 0, 0], -1)
    return img


def features_to_RGB(*Fs, skip=1):
    """Copied from pixloc repo. [LINK] https://github.com/cvg/pixloc/blob/master/pixloc/visualization/viz_2d.py"""
    """Project a list of d-dimensional feature maps to RGB colors using PCA."""
    from sklearn.decomposition import PCA

    def normalize(x):
        return x / np.linalg.norm(x, axis=-1, keepdims=True)
    flatten = []
    shapes = []
    for F in Fs:
        c, h, w = F.shape
        F = np.rollaxis(F, 0, 3)
        F = F.reshape(-1, c)
        flatten.append(F)
        shapes.append((h, w))
    flatten = np.concatenate(flatten, axis=0)

    pca = PCA(n_components=3)
    if skip > 1:
        pca.fit(normalize(flatten[::skip]))
        flatten = normalize(pca.transform(normalize(flatten)))
    else:
        flatten = normalize(pca.fit_transform(normalize(flatten)))
    flatten = (flatten + 1) / 2

    Fs = []
    for h, w in shapes:
        F, flatten = np.split(flatten, [h * w], axis=0)
        F = F.reshape((h, w, 3))
        Fs.append(F)
    assert flatten.shape[0] == 0
    return Fs


def plot_images(imgs, titles=None, cmaps='gray', dpi=100, size=6, pad=.5):
    """Plot a set of images horizontally.
    Args:
        imgs: a list of NumPy or PyTorch images, RGB (H, W, 3) or mono (H, W).
        titles: a list of strings, as titles for each image.
        cmaps: colormaps for monochrome images.
    """
    n = len(imgs)
    if not isinstance(cmaps, (list, tuple)):
        cmaps = [cmaps] * n
    figsize = (size * n, size * 3 / 4) if size is not None else None
    fig, ax = plt.subplots(1, n, figsize=figsize, dpi=dpi)
    if n == 1:
        ax = [ax]
    for i in range(n):
        ax[i].imshow(imgs[i], cmap=plt.get_cmap(cmaps[i]))
        ax[i].get_yaxis().set_ticks([])
        ax[i].get_xaxis().set_ticks([])
        ax[i].set_axis_off()
        for spine in ax[i].spines.values():  # remove frame
            spine.set_visible(False)
        if titles:
            ax[i].set_title(titles[i])
    fig.tight_layout(pad=pad)


def plot_keypoints(kpts, colors='lime', ps=2):
    """Plot keypoints for existing images.
    Args:
        kpts: list of ndarrays of size (N, 2).
        colors: string, or list of list of tuples (one for each keypoints).
        ps: size of the keypoints as float.
    """
    if not isinstance(colors, list):
        colors = [colors] * len(kpts)
    axes = plt.gcf().axes
    for a, k, c in zip(axes, kpts, colors):
        a.scatter(k[:, 0], k[:, 1], c=c, s=ps)
