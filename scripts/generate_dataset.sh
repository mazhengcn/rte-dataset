#!/bin/bash

# Set config and generator paths
CONFIG_DIR='rte_dataset/generator/matlab/2d-sweeping/configs'
GENERATOR_PATH='rte_dataset/generator/matlab/2d-sweeping'

# Set data save directory and destination host directory
DATA_SAVE_DIR= 'assets/raw_data'

# Generate data
GEN_CONFIG_PATH="${GENERATOR_PATH}/config.m"
if [ -d "${CONFIG_DIR}" ]; then
    for file in "${CONFIG_DIR}"/*; do
        echo "Generating data with config file: ${file}"
        cp "${file}" "${GEN_CONFIG_PATH}"

        matlab -nodisplay -r "run ${GENERATOR_PATH}/generator.m; exit"

        rm -f "${GEN_CONFIG_PATH}"
    done
else
    echo "Error: Config directory ${CONFIG_DIR} does not exist"
    exit 1
fi
