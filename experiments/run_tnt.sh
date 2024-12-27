tnt_output_dir=$1
tnt_meta_train_dir=$2
tnt_colmap_dir=$3
tnt_load_dir=$4
echo "[$0 arg] tnt_output_dir = ${tnt_output_dir}"
echo "[$0 arg] tnt_meta_train_dir = ${tnt_meta_train_dir}"
echo "[$0 arg] tnt_colmap_dir = ${tnt_colmap_dir}"
echo "[$0 arg] tnt_load_dir = ${tnt_load_dir}"

# detectors=("deeplsd")
detectors=("lsd" "deeplsd")
nvs=(4)

# test_scene="Barn"
test_scene="all"

for detector in ${detectors[*]}; do
    output_dir=${tnt_output_dir}/${detector}
    bash experiments/tnt/line_triangulation.sh ${tnt_colmap_dir} ${output_dir} ${detector} ${test_scene} ${tnt_load_dir}

    for nv in ${nvs[*]}; do
        bash experiments/tnt/eval.sh ${tnt_meta_train_dir} ${output_dir} ${nv} ${test_scene}
    done
done