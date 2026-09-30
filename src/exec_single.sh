#!/bin/bash

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q
#$ -pe OpenMP 16

# echo "-------------------------------------------------------------------------------------"
# echo "Number of phisical CPUs = " `grep physical.id /proc/cpuinfo | sort -u | wc -l`
# echo "Number of cores per CPU : "
# grep cpu.cores /proc/cpuinfo | sort -u
# echo "Number of logical processors = " `grep processor /proc/cpuinfo | wc -l`
# echo "System memory usage : "
# free -h
# echo "-------------------------------------------------------------------------------------"

export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

# Use ${NSLOTS} as OMP_NUM_THREADS when $HOSTNAME is "chemopgw" or "asagasumi"
# if [ "$HOSTNAME" = "chemopgw" ] || [ "$HOSTNAME" = "asagasumi" ]; then
#     export OMP_NUM_THREADS=${NSLOTS}
# else
#     export OMP_NUM_THREADS=16
# fi

if [ -n "${NSLOTS}" ]; then
    # 1. NSLOTSがあれば最優先
    export OMP_NUM_THREADS=${NSLOTS}
elif [ -n "${OMP_NUM_THREADS}" ]; then
    # 2. すでにOMP_NUM_THREADSが設定されていればそれを維持（何もしない）
    :
elif [ -n "$(nproc)" ]; then
    # 3. なければ利用可能なCPU数を取得
    export OMP_NUM_THREADS=$(nproc)
else
    # 4. それ以外（nprocが使えない環境など）はデフォルト4
    export OMP_NUM_THREADS=4
fi

# parent_dir="${PWD}"
parent_dir=$(cd $(dirname $0); pwd)

echo "Start:" `TZ=JST-9 date "+%Y/%m/%d %H:%M:%S"`


if [ $# != 1 ]; then
    echo "Usage: $0 <csv_file>"
    exit 1
elif [ ! -s "$1" ]; then
    echo "Error. No data in $1"
    exit 1
fi

chmod u+x "${parent_dir}/auto_simplex/activeshirley_gcc.exe"
chmod u+x "${parent_dir}/ga/decompositionMultipleSpectra.exe"
chmod u+x "${parent_dir}/param_simplex/active_shirley_from_params.exe"

original_file_path="$1"
original_filename=`basename "$1"`
original_dirname=`basename "${original_filename}" ".csv"`

temp_filename="_sample.csv"
cp "${original_file_path}" "${temp_filename}"

work_dir=`basename "${temp_filename}" ".csv"`
rm -rf "${original_dirname}"
mkdir -p "${work_dir}"
cp "${temp_filename}" "${work_dir}"

# Use "recipe.txt" if it exists in $PWD; otherwise, not use it.
PWD_RECIPE="${PWD}/recipe.txt"
if [ -f "${PWD_RECIPE}" ]; then
    cp "${PWD_RECIPE}" "${work_dir}"
fi

cd "${work_dir}"
    python3 ${parent_dir}/peakSeparationForXPS.py "${temp_filename}" >stdout 2>stderr
    mv "${temp_filename}" "${original_filename}"
cd ..

mv "${work_dir}" "${original_dirname}"
rm -f "${temp_filename}"

# rm -rf __pycache__

echo "End:  " `TZ=JST-9 date "+%Y/%m/%d %H:%M:%S"`
