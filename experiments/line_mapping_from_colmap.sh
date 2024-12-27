colmap_path=$1
model_path=$2
image_path=$3
output_dir=$4
echo "[$0 arg] colmap_path = ${colmap_path}"
echo "[$0 arg] model_path = ${model_path}"
echo "[$0 arg] image_path = ${image_path}"
echo "[$0 arg] output_dir = ${output_dir}"

# directory tree:
# --${colmap_path}
#     --${model_path}
#         cameras.txt/bin
#         images.txt/bin
#         points3D.txt/bin
#     --${image_path}

config_file="cfgs_linemap/triangulation/default.yaml"
default_config_file="cfgs_linemap/triangulation/default.yaml"

log_path="${output_dir}/log.txt"
mkdir -p "${output_dir}"

python runners_linemap/colmap_triangulation.py --config_file ${config_file} \
    --default_config_file ${default_config_file} \
    --colmap_path ${colmap_path} \
    --model_path ${model_path} \
    --image_path ${image_path} \
    --visualize False \
    --line2d.visualize True \
    --output_dir ${output_dir} | tee ${log_path}