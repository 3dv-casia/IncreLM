import os
import sys
import numpy as np
import pdb

MAX_ERROR = 0.01
# colmap_output_path = os.path.expanduser('~/data/TanksTemples/colmap/training')
# input_meta_path = os.path.expanduser('~/data/TanksTemples/meta_train')


def get_imname_list(colmap_output_path, scene_id):
    image_path = os.path.join(colmap_output_path, scene_id, 'dense/images')
    flist = os.listdir(image_path)
    n_images = len(flist)
    imname_list = []
    for idx in range(n_images):
        fname = '{0:06d}.jpg'.format(idx + 1)
        # fname = os.path.join(image_path, fname)
        imname_list.append(fname)
    return imname_list


def read_positions(log_file):
    with open(log_file, 'r') as f:
        lines = f.readlines()
    n_images = int(len(lines) / 5)
    positions = []
    counter = 0
    for img_id in range(n_images):
        counter += 1
        mat = []
        for idx in range(4):
            l = lines[counter].strip('\n').split(' ')
            mat.append([float(k) for k in l])
            counter += 1
        mat = np.array(mat)
        pos = mat[:3, 3]
        positions.append(pos)
    return positions


def read_trans(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    mat = []
    for idx in range(4):
        l = lines[idx].strip('\n').split(' ')
        mat.append([float(k) for k in l])
    mat = np.array(mat)
    assert np.all(mat[3, :] == np.array([0, 0, 0, 1]))
    return mat[:3, :]


def write_geoinfo_txt(fname, imname_list, positions):
    with open(fname, 'w') as f:
        for imname, pos in zip(imname_list, positions):
            f.write('{0} {1} {2} {3}\n'.format(imname, pos[0], pos[1], pos[2]))


def parse_config():
    import argparse
    arg_parser = argparse.ArgumentParser(
        description='align a Tanks and Temples scene')
    arg_parser.add_argument('--scene_id', type=str, required=True, help='scene_id')
    arg_parser.add_argument('--input_meta_path', type=str, required=True, help='meta train path')
    arg_parser.add_argument('--colmap_output_path', type=str, required=True, help='colmap output path')
    args, unknown = arg_parser.parse_known_args()
    cfg = dict()
    cfg["scene_id"] = args.scene_id
    cfg["input_meta_path"] = args.input_meta_path
    cfg["colmap_output_path"] = args.colmap_output_path
    return cfg


def main():
    cfg = parse_config()
    scene_id = cfg["scene_id"]
    input_meta_path = cfg["input_meta_path"]
    colmap_output_path = cfg["colmap_output_path"]

    # get geo txt
    imname_list = get_imname_list(colmap_output_path, scene_id)
    log_file = os.path.join(input_meta_path, scene_id, '{0}_COLMAP_SfM.log'.format(scene_id))
    positions = read_positions(log_file)
    trans_file = os.path.join(input_meta_path, scene_id, '{0}_trans.txt'.format(scene_id))
    trans_mat = read_trans(trans_file)
    new_positions = [trans_mat[:3, :3] @ k + trans_mat[:3, 3] for k in positions]
    output_fname = os.path.join(input_meta_path, scene_id, 'geo_positions.txt')
    write_geoinfo_txt(output_fname, imname_list, new_positions)

    # colmap align
    cmd_list = []
    basepath = os.path.join(colmap_output_path, scene_id, 'dense')
    cmd = 'mkdir -p {0}'.format(os.path.join(basepath, 'aligned'))
    cmd_list.append(cmd)
    cmd = 'colmap model_aligner --input_path {0} --output_path {1} --ref_images_path {2} --robust_alignment 1 --robust_alignment_max_error {3} --transform_path {4} --ref_is_gps false'.format(
        os.path.join(basepath, 'sparse'), os.path.join(basepath, 'aligned'), os.path.join(input_meta_path, scene_id, 'geo_positions.txt'), MAX_ERROR, os.path.join(basepath, 'transform.txt'))
    cmd_list.append(cmd)
    cmd = 'colmap model_converter --input_path {0} --output_path {1} --output_type PLY'.format(
        os.path.join(basepath, 'aligned'), os.path.join(basepath, 'aligned/points.ply'))
    cmd_list.append(cmd)
    cmd = 'colmap model_converter --input_path {0} --output_path {1} --output_type TXT'.format(
        os.path.join(basepath, 'aligned'), os.path.join(basepath, 'aligned'))
    cmd_list.append(cmd)
    for cmd in cmd_list:
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    main()
